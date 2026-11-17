function analysis=DSGE_analyze(root)
% Recompute all quantitative slide evidence without reestimating parameters.
DSGE_setup;
resultDir=fullfile(root,'results'); work=fullfile(resultDir,'dynare_work');
addpath(work); old=pwd; cleanup=onCleanup(@() cd(old)); cd(work);
S=load(fullfile(resultDir,'estimated_model.mat'));
D=load(fullfile(root,'data','dsge_observables.mat'));
diagFile=fullfile(resultDir,'mcmc_diagnostics.mat');
if isfield(S,'chain_work') && isfolder(S.chain_work)
    diagnostics=DSGE_mcmc_diagnostics(S.chain_work,S.bayestopt_.name,.5);
    save(diagFile,'diagnostics','-v7.3');
else
    T=load(diagFile,'diagnostics'); diagnostics=T.diagnostics;
end
writetable(diagnostics.table,fullfile(resultDir,'posterior_summary.csv'));
decompVars={'GDPobs','Cobs','Iobs','N_obs','Piobs','Robs'};
sourceFields={'GDPobs','Cobs','Iobs','Nobs','Piobs','Robs'};
S.options_.shock_decomp.parameter_set='posterior_mean';
S.options_.shock_decomp.nograph=true; S.options_.nograph=true;
[S.oo_,S.M_]=shock_decomposition(S.M_,S.oo_,S.options_,decompVars,S.bayestopt_,S.estim_params_);
S.options_.irf=40; S.options_.order=1; S.options_.noprint=true;
vars={'GDPobs','Cobs','Iobs','N_obs','wobs','Piobs','Robs','Kobs','uobs'};
[info,S.oo_,S.options_,S.M_]=stoch_simul(S.M_,S.options_,S.oo_,vars);
assert(info(1)==0);
save(fullfile(resultDir,'estimated_model.mat'),'-struct','S','-v7.3');
bands=DSGE_irf_bands(S,diagnostics,vars,2000);
save(fullfile(resultDir,'irf_posterior_bands.mat'),'bands','-v7');
shockNames={'eA','eG','eR','eb','eI','epmark','ewmark'};
irfs=bands.posterior_mean; fixed=zeros(size(irfs));
for k=1:7, for j=1:numel(vars), fixed(:,j,k)=S.oo_.irfs.([vars{j},'_',shockNames{k}])(:); end; end
vd=DSGE_variance_decomposition(S,decompVars);
comparison=DSGE_moments(S,D,decompVars,sourceFields);
assert(max(abs(comparison.model_sd-vd.standard_deviation))<1e-10);
save(fullfile(resultDir,'variance_decomposition.mat'),'vd');
save(fullfile(resultDir,'moments_comparison.mat'),'comparison');
writetable(comparison.table,fullfile(resultDir,'moments_comparison.csv'));
tab=array2table(vd.shares_percent,'VariableNames',shockNames);
tab=addvars(tab,string(decompVars(:)),vd.standard_deviation,'Before',1, ...
    'NewVariableNames',{'Variable','ModelSD'});
writetable(tab,fullfile(resultDir,'variance_decomposition.csv'));
n=numel(D.GDPobs); observed=zeros(n,6); reconstructed=observed; decomposition=zeros(n,8,6);
for j=1:6
    index=find(strcmp(S.M_.endo_names,decompVars{j}));
    pieces=squeeze(S.oo_.shock_decomposition(index,:,:))';
    decomposition(:,:,j)=pieces(:,1:8); reconstructed(:,j)=sum(pieces(:,1:8),2);
    observed(:,j)=D.(sourceFields{j});
    assert(max(abs(reconstructed(:,j)-observed(:,j)))<1e-5,'Historical decomposition does not add up.');
end
analysis=struct('irfs',irfs,'irf_variables',{vars},'shock_names',{shockNames}, ...
    'irf_summary','pointwise posterior mean','irfs_at_posterior_mean_parameters',fixed, ...
    'irf_bands',bands,'variance_decomposition',vd,'moments_comparison',comparison, ...
    'historical_contributions',decomposition,'historical_variables',{decompVars}, ...
    'observed',observed,'reconstructed',reconstructed,'dates',D.quarter_dates, ...
    'parameter_set','posterior_mean','posterior_mean',S.oo_.posterior_mean, ...
    'max_reconstruction_error',max(abs(observed-reconstructed),[],'all'));
save(fullfile(resultDir,'analysis_results.mat'),'analysis','-v7');
DSGE_mechanisms(root);
DSGE_parameter_table(root,S,diagnostics);
DSGE_slides_figures('all',root);
fprintf('Analysis complete. Historical reconstruction error %.3g.\n',analysis.max_reconstruction_error);
end
