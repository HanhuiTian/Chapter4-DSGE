function comparison=DSGE_moments(S,D,variables,sourceFields)
% Compare unconditional theoretical moments with observed sample moments.
M=S.M_; O=S.oo_; options=S.options_;
options.order=1; options.ar=1; options.hp_filter=0;
options.one_sided_hp_filter=0; options.bandpass.indicator=false;
if isempty(options.qz_criterium), options.qz_criterium=1+1e-6; end
indices=zeros(numel(variables),1); X=zeros(numel(D.GDPobs),numel(variables));
for j=1:numel(variables)
    indices(j)=find(strcmp(M.endo_names,variables{j}));
    X(:,j)=D.(sourceFields{j});
end
assert(all(isfinite(X),'all'));
[moments,stationary]=th_autocovariances(O.dr,indices,M,options,1);
assert(numel(stationary)==numel(variables));
covariance=moments{1}; modelSD=sqrt(diag(covariance));
modelCorr=covariance(:,1)./(modelSD*modelSD(1));
modelACF=diag(moments{2});
dataSD=std(X,0,1)'; correlation=corrcoef(X); dataCorr=correlation(:,1);
dataACF=zeros(numel(variables),1);
for j=1:numel(variables)
    lagged=corrcoef(X(2:end,j),X(1:end-1,j)); dataACF(j)=lagged(1,2);
end
% Independent covariance and lag covariance check from impulse responses.
irfCov=zeros(numel(variables)); irfLag=irfCov;
for k=1:M.exo_nbr
    impulse=zeros(M.exo_nbr,1); impulse(k)=sqrt(M.Sigma_e(k,k));
    path=irf(M,options,O.dr,impulse,10000,0,1,1); path=path(indices,:);
    irfCov=irfCov+path*path';
    irfLag=irfLag+path(:,2:end)*path(:,1:end-1)';
end
covError=max(abs(irfCov-covariance)./(modelSD*modelSD'),[],'all');
acfError=max(abs(diag(irfLag)./diag(irfCov)-modelACF));
assert(covError<1e-6 && acfError<1e-6,'Moment crosscheck failed.');
assert(all(modelSD>0) && all(abs([modelCorr;modelACF;dataCorr;dataACF])<=1+1e-10));
summary=table(string(variables(:)),dataSD,modelSD,dataCorr,modelCorr,dataACF,modelACF, ...
    'VariableNames',{'Variable','DataSD','ModelSD','DataCorrGDP','ModelCorrGDP','DataACF1','ModelACF1'});
comparison=struct('variables',{variables},'source_fields',{sourceFields}, ...
    'data_sd',dataSD,'model_sd',modelSD,'data_corr_gdp',dataCorr, ...
    'model_corr_gdp',modelCorr,'data_acf1',dataACF,'model_acf1',modelACF, ...
    'model_covariance',covariance,'sample_size',size(X,1),'dates',D.quarter_dates, ...
    'parameter_set','posterior_mean','model_moment_type','unconditional theoretical', ...
    'data_sd_normalization','T-1','data_acf_definition','Pearson correlation of adjacent observations', ...
    'irf_covariance_crosscheck_error',covError,'irf_acf_crosscheck_error',acfError,'table',summary);
disp(summary);
end
