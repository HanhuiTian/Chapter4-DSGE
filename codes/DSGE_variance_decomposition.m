function vd=DSGE_variance_decomposition(S,variables)
% Unconditional observation variances, at the posterior mean, without refiltering.
M=S.M_; O=S.oo_; options=S.options_;
options.order=1; options.ar=0; options.hp_filter=0;
options.one_sided_hp_filter=0; options.bandpass.indicator=false;
if isempty(options.qz_criterium), options.qz_criterium=1+1e-6; end
indices=zeros(numel(variables),1);
for j=1:numel(variables), indices(j)=find(strcmp(M.endo_names,variables{j})); end
assert(max(abs(M.Sigma_e-diag(diag(M.Sigma_e))),[],'all')<1e-12);
[moments,stationary]=th_autocovariances(O.dr,indices,M,options,0);
assert(numel(stationary)==numel(variables),'Nonstationary observation in variance decomposition.');
variance=diag(moments{1}); shares=moments{2};
assert(all(variance>0) && all(isfinite(shares),'all'));
assert(max(abs(sum(shares,2)-1))<1e-8 && min(shares,[],'all')>-1e-10);
contributions=shares.*variance;
% Independent check: sum squared impulse responses over 10000 quarters.
irfVariance=zeros(size(contributions));
for k=1:M.exo_nbr
    impulse=zeros(M.exo_nbr,1); impulse(k)=sqrt(M.Sigma_e(k,k));
    path=irf(M,options,O.dr,impulse,10000,0,1,1);
    irfVariance(:,k)=sum(path(indices,:).^2,2);
end
relativeError=max(abs(irfVariance-contributions)./max(variance,eps),[],'all');
assert(relativeError<1e-6,'IRF and Lyapunov variance calculations disagree.');
vd=struct('variables',{variables},'shock_names',{cellstr(M.exo_names)}, ...
    'variance',variance,'standard_deviation',sqrt(variance), ...
    'contributions',contributions,'shares_percent',100*shares, ...
    'parameter_set','posterior_mean','horizon','unconditional', ...
    'irf_crosscheck_relative_error',relativeError);
end
