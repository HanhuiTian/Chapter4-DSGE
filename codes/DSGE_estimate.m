function DSGE_estimate(action, draws, root)
% Internal estimation driver. Use DSGE_main('estimate') as the entry point.
if nargin<1, action='estimate'; end
if nargin<2, draws=100000; end
assert(nargin>=3,'Use DSGE_main to select an isolated output workspace.');
DSGE_setup;
oldThreads=maxNumCompThreads(1);
threadCleanup=onCleanup(@() maxNumCompThreads(oldThreads));
work=fullfile(root,'results','dynare_work');
if ~isfolder(work), mkdir(work); end
addpath(fileparts(mfilename('fullpath')),fullfile(root,'data'),work);
if strcmpi(action,'estimate')
    if isfile(fullfile(root,'results','estimated_model.mat'))
        DSGE_analyze(root);
        return
    end
    if ~isfile(fullfile(root,'results','posterior_mode.mat')), DSGE_estimate('mode',draws,root); end
    if ~isfile(fullfile(root,'results','multistart_audit.mat')), DSGE_estimate('multistart',draws,root); end
    pilotFile=fullfile(root,'results','pilot_rw','metropolis','medium_dsge_mh1_blck1.mat');
    if ~isfile(pilotFile), DSGE_estimate('pilot',100000,root); end
    DSGE_mixture_estimate(root,draws);
    return
end
if strcmpi(action,'multistart')
    DSGE_multistart(root);
    return
end
if strcmpi(action,'mixture')
    DSGE_mixture_estimate(root,draws);
    return
end
if strcmpi(action,'extend') && isfile(fullfile(root,'results','mh_mixture','proposal.mat'))
    DSGE_mixture_estimate(root,draws);
    return
end
if strcmpi(action,'analyze')
    DSGE_analyze(root);
    return
end
if strcmpi(action,'mechanisms')
    DSGE_mechanisms(root);
    return
end
DSGE_setup;
originalDir=pwd;
cleanup=onCleanup(@() cd(originalDir));
D=load(fullfile(root,'data','dsge_observables.mat'));
data_beta=D.calibration.beta; data_Pi=D.calibration.Pi;
data_sC=D.calibration.s_C; data_sI=D.calibration.s_I; data_sG=D.calibration.s_G;
save(fullfile(work,'calibration_input.mat'),'data_beta','data_Pi','data_sC','data_sI','data_sG');
% Dynare's parser reserves NOBS (case insensitive). Preserve the source file
% and give the identical labor observation an unambiguous model alias.
D.N_obs=D.Nobs;
save(fullfile(work,'dsge_observables.mat'),'-struct','D','-v7');
copyfile(fullfile(root,'medium_dsge.mod'),fullfile(work,'medium_dsge.mod'));
cd(work);
if strcmpi(action,'check')
    flags='-DRUN_ESTIMATION=0';
elseif strcmpi(action,'mode')
    flags='-DRUN_ESTIMATION=1 -DMH_DRAWS=0 -DRESUME_MODE=0';
elseif ismember(lower(action),{'native_estimate','extend','pilot'})
    resume=isfile(fullfile(work,'medium_dsge','Output','medium_dsge_mode.mat'));
    extend=strcmpi(action,'extend');
    if ~extend && ~strcmpi(action,'pilot') && isfile(fullfile(root,'results','estimated_model.mat'))
        error('Saved estimates exist. Use analyze to reuse them or extend to add draws.');
    end
    if extend
        assert(isfolder(fullfile(work,'medium_dsge','metropolis')),'No existing chains to extend.');
    end
    flags=sprintf('-DRUN_ESTIMATION=1 -DMH_DRAWS=%d -DRESUME_MODE=%d -DEXTEND_CHAINS=%d',draws,resume,extend);
    if strcmpi(action,'pilot'), flags=[flags,' -DMH_CHAINS=1']; end
else
    error('Unknown action: %s',action);
end
evalin('base',['dynare medium_dsge.mod noclearall ',flags]);
snapshot=struct();
for name={'M_','oo_','options_','estim_params_','bayestopt_','dataset_','dataset_info'}
    if evalin('base',sprintf('exist(''%s'',''var'')',name{1}))
        snapshot.(name{1})=evalin('base',name{1});
    end
end
snapshot.run_action=action;
snapshot.finished_at=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
snapshot.matlab_version=version;
if strcmpi(action,'check')
    file='validation.mat';
elseif strcmpi(action,'mode')
    file='posterior_mode.mat';
elseif strcmpi(action,'pilot')
    file='pilot_estimation.mat';
else
    file='estimated_model.mat';
end
save(fullfile(root,'results',file),'-struct','snapshot','-v7.3');
fprintf('Saved %s\n',fullfile(root,'results',file));
if strcmpi(action,'pilot')
    pilotDir=fullfile(root,'results','pilot_rw','metropolis');
    if ~isfolder(pilotDir), mkdir(pilotDir); end
    copyfile(fullfile(work,'medium_dsge','metropolis','medium_dsge_mh1_blck1.mat'),pilotDir);
end
if ismember(lower(action),{'native_estimate','extend'}) && draws>0, DSGE_analyze(root); end
end
