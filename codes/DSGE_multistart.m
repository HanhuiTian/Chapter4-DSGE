function DSGE_multistart(root)
% Refine the Monte Carlo mode and check four independent feasible prior starts.
if nargin<1, root=fileparts(mfilename('fullpath')); end
DSGE_setup;
work=fullfile(root,'results','dynare_work'); addpath(work);
previous=pwd; cleanup=onCleanup(@() cd(previous)); cd(work);
S=load(fullfile(root,'results','posterior_mode.mat'));
modeFile=fullfile('medium_dsge','Output','medium_dsge_mode.mat');
saved=load(modeFile);
options=S.options_; options.mode_compute=4; options.mode_file=modeFile;
options.mh_replic=0; options.csminwel.maxiter=500;
options.csminwel.tolerance.f=1e-7; options.silent_optimizer=true;
[dataset,datasetInfo,~,~,M,options,O,E,B,bounds]= ...
    dynare_estimation_init(options.varobs,M_name(S.M_),[],S.M_,options,S.oo_,S.estim_params_,S.bayestopt_);
args={dataset,datasetInfo,options,M,E,B,bounds,O.dr,O.steady_state,O.exo_steady_state,O.exo_det_steady_state};
rng(20261112,'twister'); prior_draw(B,options.prior_trunc);
starts=zeros(numel(saved.xparam1),5); starts(:,1)=saved.xparam1;
for k=2:5
    for tries=1:10000
        x=prior_draw; x=x(:); [value,~,valid]=dsge_likelihood(x,args{:});
        if valid==1 && isfinite(value), starts(:,k)=x; break; end
    end
    assert(any(starts(:,k)),'Could not draw feasible starting parameters.');
end
best=saved.xparam1; bestValue=dsge_likelihood(best,args{:});
solutions=zeros(size(starts)); objectives=zeros(5,1); exitflags=zeros(5,1);
for k=1:5
    [x,f,flag]=dynare_minimize_objective(@dsge_likelihood,starts(:,k),4,options, ...
        [bounds.lb bounds.ub],B.name,B,[],args{:});
    solutions(:,k)=x; objectives(k)=f; exitflags(k)=flag;
    fprintf('Mode cross-check %d/5: negative log posterior %.8f, exit %d\n',k,f,flag);
    if f<bestValue, best=x; bestValue=f; end
end
audit=struct('starts',starts,'solutions',solutions,'objectives',objectives, ...
    'exitflags',exitflags,'best',best,'best_value',bestValue);
save(fullfile(root,'results','multistart_audit.mat'),'audit');
copyfile(modeFile,fullfile(root,'results','mode_before_multistart.mat'));
xparam1=best; fval=bestValue; parameter_names=B.name;
hh=hessian(@dsge_likelihood,xparam1,options.gstep,args{:});
hh=reshape(hh,numel(best),numel(best)); hh=(hh+hh')/2;
[~,flag]=chol(hh);
if flag~=0
    hh=saved.hh;
    warning('Mode Hessian not positive definite. Keeping Monte Carlo proposal covariance.');
end
save(modeFile,'xparam1','fval','parameter_names','hh');
fprintf('Best mode saved. Negative log posterior %.8f.\n',fval);
end

function name=M_name(M), name=M.dname; end
