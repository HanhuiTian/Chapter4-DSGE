function checks=DSGE_verify
% Compare a fresh published-result run with the bundled numerical reference.
package=fileparts(fileparts(mfilename('fullpath')));
root=fullfile(package,'output','published','results');
A=load(fullfile(root,'analysis_results.mat')); B=load(fullfile(package,'reference_results','analysis_results.mat'));
R=load(fullfile(root,'mechanism_comparison.mat')); T=load(fullfile(package,'reference_results','mechanism_comparison.mat'));
checks.irf=max(abs(A.analysis.irfs-B.analysis.irfs),[],'all');
checks.hpd=max(abs([A.analysis.irf_bands.low(:)-B.analysis.irf_bands.low(:); ...
    A.analysis.irf_bands.high(:)-B.analysis.irf_bands.high(:)]));
checks.history=max(abs(A.analysis.historical_contributions-B.analysis.historical_contributions),[],'all');
checks.variance=max(abs(A.analysis.variance_decomposition.shares_percent-B.analysis.variance_decomposition.shares_percent),[],'all');
checks.moments=max(abs(table2array(A.analysis.moments_comparison.table(:,2:end))- ...
    table2array(B.analysis.moments_comparison.table(:,2:end))),[],'all');
checks.mechanisms=max(abs(R.result.posterior_mean-T.result.posterior_mean),[],'all');
assert(max(struct2array(checks))<1e-7,'Published replication differs from the reference.');
assert(strcmp(R.result.shock_name,'eR') && R.result.number_of_draws==2000);
disp(checks); fprintf('PASS: all published quantitative results reproduced.\n');
end
