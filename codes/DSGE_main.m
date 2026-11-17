function result = DSGE_main(action, source, draws)
%DSGE_MAIN Chapter 4 replication entry point.
% DSGE_main                       Recompute analyses from published estimates.
% DSGE_main('figures')            Redraw figures from saved numerical results.
% DSGE_main('check')              Check the published model solution.
% DSGE_main('estimate')           Estimate afresh in output/reestimated.
% DSGE_main('estimate',[],100000)  Set the number of draws per chain.
% DSGE_main('analyze','reestimated') Analyze your own estimation results.
% DSGE_main('extend',[],100000)    Extend your own four chains.
% DSGE_main('data')               Rebuild observations from the bundled Excel.
if nargin<1 || isempty(action), action='analyze'; end
if nargin<2 || isempty(source), source='published'; end
if nargin<3 || isempty(draws), draws=100000; end
action=lower(char(action)); source=lower(char(source));
codeDir=fileparts(mfilename('fullpath')); package=fileparts(codeDir);
addpath(codeDir,fullfile(package,'data'));
assert(ismember(source,{'published','reestimated'}),'Unknown estimation source.');
assert(ismember(action,{'analyze','figures','check','estimate','extend','data'}),'Unknown action.');
if strcmp(action,'data')
    dest=fullfile(package,'output','data'); if ~isfolder(dest), mkdir(dest); end
    result=prepare_medium_dsge_data(fullfile(package,'data','us_dsge_raw.xlsx'), ...
        fullfile(dest,'dsge_observables.mat')); return
end
DSGE_setup;
oldDir=pwd; oldPath=path; oldRng=rng; oldThreads=maxNumCompThreads(1);
cleanup=onCleanup(@() restore_session(oldDir,oldPath,oldRng,oldThreads));
if ismember(action,{'estimate','extend'}), source='reestimated'; end
root=fullfile(package,'output',source);
for folder={'results','data','figures'}
    dest=fullfile(root,folder{1}); if ~isfolder(dest), mkdir(dest); end
end
copyfile(fullfile(package,'data','dsge_observables.mat'),fullfile(root,'data','dsge_observables.mat'));
copyfile(fullfile(codeDir,'medium_dsge.mod'),fullfile(root,'medium_dsge.mod'));
if strcmp(source,'published')
    % One-way copies only. No routine ever saves to estimates/ or reference_results/.
    for file={'estimated_model.mat','mcmc_diagnostics.mat','posterior_summary.csv'}
        copyfile(fullfile(package,'estimates',file{1}),fullfile(root,'results',file{1}));
    end
    if strcmp(action,'figures')
        files=dir(fullfile(package,'reference_results','*.mat'));
        for j=1:numel(files)
            copyfile(fullfile(files(j).folder,files(j).name),fullfile(root,'results',files(j).name));
        end
    end
end
switch action
    case 'figures'
        DSGE_slides_figures('all',root);
    case {'analyze','check'}
        assert(isfile(fullfile(root,'results','estimated_model.mat')), ...
            'No estimates in this workspace. Run DSGE_main(''estimate'') first.');
        DSGE_prepare(root);
        if strcmp(action,'analyze'), DSGE_analyze(root); end
    case 'estimate'
        DSGE_estimate('estimate',draws,root);
    case 'extend'
        assert(isfile(fullfile(root,'results','mh_mixture','proposal.mat')), ...
            'No reestimated chains to extend. Run estimate first.');
        DSGE_estimate('extend',draws,root);
end
result=struct('source',source,'output_directory',root,'figures_directory',fullfile(root,'figures'));
fprintf('\nResults: %s\nFigures: %s\n',root,result.figures_directory);
end

function restore_session(folder,searchPath,state,threads)
cd(folder); path(searchPath); rng(state); maxNumCompThreads(threads);
end
