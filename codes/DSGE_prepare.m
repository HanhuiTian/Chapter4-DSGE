function DSGE_prepare(root)
% Generate platform-local model functions and bind a working estimate copy.
DSGE_setup;
work=fullfile(root,'results','dynare_work'); if ~isfolder(work), mkdir(work); end
D=load(fullfile(root,'data','dsge_observables.mat'));
data_beta=D.calibration.beta; data_Pi=D.calibration.Pi;
data_sC=D.calibration.s_C; data_sI=D.calibration.s_I; data_sG=D.calibration.s_G;
save(fullfile(work,'calibration_input.mat'),'data_beta','data_Pi','data_sC','data_sI','data_sG');
D.N_obs=D.Nobs; save(fullfile(work,'dsge_observables.mat'),'-struct','D','-v7');
copyfile(fullfile(root,'medium_dsge.mod'),fullfile(work,'medium_dsge.mod'));
old=pwd; cleanup=onCleanup(@() cd(old)); cd(work); addpath(work);
evalin('base','dynare medium_dsge.mod noclearall -DRUN_ESTIMATION=0');
fresh=evalin('base','M_');
[~,source]=fileparts(root);
if strcmp(source,'published')
    package=fileparts(fileparts(mfilename('fullpath')));
    copyfile(fullfile(package,'estimates','posterior_mode.mat'), ...
        fullfile(work,'medium_dsge','Output','medium_dsge_mode.mat'));
end
S=load(fullfile(root,'results','estimated_model.mat'));
assert(isequal(fresh.endo_names,S.M_.endo_names) && isequal(fresh.param_names,S.M_.param_names), ...
    'Model variables/parameters changed. Estimate the revised model in a new workspace.');
if isfield(S,'chain_work'), S=rmfield(S,'chain_work'); end
S.options_.datafile='dsge_observables.mat';
S.options_.dirname='medium_dsge'; S.M_.dname='medium_dsge'; S.M_.fname='medium_dsge';
S.options_.nograph=true; S.options_.nodisplay=true;
[dr,info,S.M_.params]=resol(0,S.M_,S.options_,S.oo_.dr,S.oo_.steady_state, ...
    S.oo_.exo_steady_state,S.oo_.exo_det_steady_state);
assert(info(1)==0,'Saved posterior mean fails the solution conditions.');
residual=feval('medium_dsge.static_resid',dr.ys,zeros(1,S.M_.exo_nbr),S.M_.params);
assert(max(abs(residual))<1e-8,'Steady-state residual check failed.');
S.oo_.dr=dr; S.oo_.steady_state=dr.ys;
save(fullfile(root,'results','estimated_model.mat'),'-struct','S','-v7.3');
fprintf('Published model check: steady-state residual %.3g; BK solution successful.\n',max(abs(residual)));
end
