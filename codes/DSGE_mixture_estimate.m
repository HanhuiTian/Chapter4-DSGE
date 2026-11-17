function DSGE_mixture_estimate(root,draws)
% Frozen mixture independence MH plus random walk MH, using Dynare likelihood.
% A pilot is used only to fit the proposal. Its draws are not final posterior draws.
% Positive parameters are logged, unit-interval parameters use logit transforms.
% The Jacobian and independence-proposal density ratio enter every acceptance test.
if nargin<1, root=fileparts(mfilename('fullpath')); end
if nargin<2, draws=100000; end
DSGE_setup;
work=fullfile(root,'results','dynare_work'); addpath(work);
previous=pwd; cleanup=onCleanup(@() cd(previous)); cd(work);
S=load(fullfile(root,'results','posterior_mode.mat'));
options=S.options_; options.mode_compute=0;
options.mode_file='medium_dsge/Output/medium_dsge_mode.mat'; options.mh_replic=0;
[dataset,datasetInfo,~,~,M,options,O,E,B,bounds]= ...
    dynare_estimation_init(options.varobs,S.M_.dname,[],S.M_,options,S.oo_,S.estim_params_,S.bayestopt_);
args={dataset,datasetInfo,options,M,E,B,bounds,O.dr,O.steady_state,O.exo_steady_state,O.exo_det_steady_state};
unitMask=B.pshape==1;
mixtureWork=fullfile(root,'results','mh_mixture');
chainDir=fullfile(mixtureWork,'medium_dsge','metropolis');
if ~isfolder(chainDir), mkdir(chainDir); end
proposalFile=fullfile(mixtureWork,'proposal.mat');
if isfile(proposalFile)
    pack=load(proposalFile); proposal=pack.proposal; starts=pack.starts;
else
    pilot=load(fullfile(root,'results','pilot_rw','metropolis','medium_dsge_mh1_blck1.mat'),'x2');
    training=to_unconstrained(pilot.x2(1:10:end,:),unitMask);
    rng(20261113,'twister'); proposal=fit_proposal(training,12);
    starts=pilot.x2(round(linspace(10000,size(pilot.x2,1),4)),:);
    save(proposalFile,'proposal','starts','unitMask');
end
if license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool)
        pool=parpool('local',4);
        poolCleanup=onCleanup(@() delete(pool));
    end
    assert(pool.NumWorkers>=4,'Four MATLAB process workers are required.');
    parfor b=1:4
        mixture_chain(b,draws,starts(b,:),proposal,unitMask,args,chainDir,work);
    end
else
    for b=1:4
        mixture_chain(b,draws,starts(b,:),proposal,unitMask,args,chainDir,work);
    end
end
diagnostics=DSGE_mcmc_diagnostics(mixtureWork,B.name,.5);
disp(diagnostics.table);
fprintf('MIXTURE MH: max Rhat %.5f, min bulk ESS %.1f, min tail ESS %.1f\n', ...
    max(diagnostics.rhat),min(diagnostics.bulk_ess),min(diagnostics.tail_ess));
S.M_=set_all_parameters(diagnostics.posterior_mean,E,M);
for j=1:numel(B.name)
    S.oo_.posterior_mean.parameters.(B.name{j})=diagnostics.posterior_mean(j);
    S.oo_.posterior_median.parameters.(B.name{j})=diagnostics.posterior_median(j);
    S.oo_.posterior_hpdinf.parameters.(B.name{j})=diagnostics.hpd90(j,1);
    S.oo_.posterior_hpdsup.parameters.(B.name{j})=diagnostics.hpd90(j,2);
end
S.options_=options; S.options_.mh_replic=diagnostics.raw_draws_per_chain(1);
S.options_.shock_decomp.parameter_set='posterior_mean';
S.options_.shock_decomp.nograph=true; S.options_.nograph=true;
[S.oo_,S.M_]=shock_decomposition(S.M_,S.oo_,S.options_, ...
    {'GDPobs','Cobs','Iobs','N_obs','Piobs','Robs'},B,E);
S.options_.irf=40; S.options_.order=1; S.options_.noprint=true;
[info,S.oo_,S.options_,S.M_]=stoch_simul(S.M_,S.options_,S.oo_, ...
    {'GDPobs','Cobs','Iobs','N_obs','Piobs','Robs','wobs','Kobs','uobs','qobs'});
assert(info(1)==0,'Posterior mean does not solve the model.');
S.bayestopt_=B; S.estim_params_=E;
S.sampler='Frozen Student t mixture independence MH (85%) and random walk MH (15%) in transformed parameters';
S.chain_work=mixtureWork;
S.finished_at=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
save(fullfile(root,'results','estimated_model.mat'),'-struct','S','-v7.3');
DSGE_analyze(root);
end

function mixture_chain(b,draws,initial,proposal,unitMask,args,folder,work)
maxNumCompThreads(1);
DSGE_setup; addpath(work);
cd(work);
file=fullfile(folder,sprintf('medium_dsge_mh1_blck%d.mat',b));
if isfile(file)
    old=load(file); rng(old.rngState); initial=old.x2(end,:);
    history=old.x2; oldLog=old.logpo2;
else
    rng(20261120+b,'twister'); history=[]; oldLog=[];
end
z=to_unconstrained(initial,unitMask)';
[target,theta,logpost]=evaluate_target(z,unitMask,args);
assert(isfinite(target),'Invalid initial point.');
qCurrent=mixture_logpdf(z,proposal);
xnew=zeros(draws,numel(z)); lognew=zeros(draws,1);
accepted=0; independentAccepted=0; independentTried=0;
fid=fopen(fullfile(folder,sprintf('chain%d_progress.log',b)),'a');
cleanup=onCleanup(@() fclose(fid));
for t=1:draws
    independent=rand<.85;
    if independent
        independentTried=independentTried+1;
        candidate=mixture_draw(proposal);
    else
        candidate=z+proposal.rwScale*proposal.globalChol*randn(numel(z),1);
    end
    [proposed,thetaNew,postNew]=evaluate_target(candidate,unitMask,args);
    qNew=mixture_logpdf(candidate,proposal);
    ratio=proposed-target;
    if independent, ratio=ratio+qCurrent-qNew; end
    if isfinite(proposed) && isfinite(qNew) && log(rand)<min(0,ratio)
        z=candidate; target=proposed; theta=thetaNew; logpost=postNew;
        qCurrent=qNew; accepted=accepted+1;
        independentAccepted=independentAccepted+double(independent);
    end
    xnew(t,:)=theta'; lognew(t)=logpost;
    if mod(t,5000)==0
        fprintf(fid,'%s: %d/%d, acceptance %.4f, independent acceptance %.4f\n', ...
            char(datetime('now','Format','HH:mm:ss')),t,draws,accepted/t,independentAccepted/max(1,independentTried));
        % Closing/reopening flushes progress so external monitoring sees it.
        fclose(fid); fid=fopen(fullfile(folder,sprintf('chain%d_progress.log',b)),'a');
        x2=[history;xnew(1:t,:)]; logpo2=[oldLog;lognew(1:t)]; rngState=rng;
        save(fullfile(folder,sprintf('checkpoint_blck%d.mat',b)),'x2','logpo2','rngState');
    end
end
x2=[history;xnew]; logpo2=[oldLog;lognew]; rngState=rng;
save(file,'x2','logpo2','rngState','accepted','independentAccepted','independentTried');
end

function z=to_unconstrained(theta,unitMask)
z=log(theta);
z(:,unitMask)=log(theta(:,unitMask))-log1p(-theta(:,unitMask));
end

function [target,theta,logpost]=evaluate_target(z,unitMask,args)
theta=exp(z);
theta(unitMask)=1./(1+exp(-z(unitMask)));
target=-Inf; logpost=-Inf;
if any(~isfinite(theta)) || any(theta<=0) || any(theta(unitMask)>=1), return; end
[value,~,valid]=dsge_likelihood(theta,args{:});
if valid~=1 || ~isfinite(value), return; end
jacobian=sum(log(theta))+sum(log1p(-theta(unitMask)));
logpost=-value; target=logpost+jacobian;
end

function proposal=fit_proposal(x,K)
% K-means initialisation followed by Gaussian EM, then heavy-tailed components.
[n,d]=size(x); globalCov=cov(x)+1e-6*eye(d); L=chol(globalCov,'lower');
w=(x-mean(x,1))/L'; centers=w(randperm(n,K),:);
for iter=1:30
    distance=zeros(n,K);
    for k=1:K, distance(:,k)=sum((w-centers(k,:)).^2,2); end
    [~,label]=min(distance,[],2);
    for k=1:K
        if sum(label==k)>d+2, centers(k,:)=mean(w(label==k,:),1); end
    end
end
mu=zeros(d,K); covariance=zeros(d,d,K); weight=zeros(1,K);
for k=1:K
    select=label==k; weight(k)=max(sum(select),d+2);
    if sum(select)>d+2
        mu(:,k)=mean(x(select,:),1)'; covariance(:,:,k)=cov(x(select,:))+.03*globalCov;
    else
        mu(:,k)=mean(x,1)'; covariance(:,:,k)=globalCov;
    end
end
weight=weight/sum(weight);
for iter=1:35
    logResponsibility=zeros(n,K);
    for k=1:K
        C=chol(covariance(:,:,k),'lower'); v=(x-mu(:,k)')/C';
        logResponsibility(:,k)=log(weight(k))-sum(log(diag(C)))-.5*sum(v.^2,2);
    end
    logResponsibility=logResponsibility-max(logResponsibility,[],2);
    responsibility=exp(logResponsibility); responsibility=responsibility./sum(responsibility,2);
    mass=sum(responsibility,1); weight=mass/sum(mass);
    for k=1:K
        mu(:,k)=(responsibility(:,k)'*x/mass(k))';
        v=x-mu(:,k)'; covariance(:,:,k)=v'*(v.*responsibility(:,k))/mass(k)+.02*globalCov;
    end
end
proposal=struct('weight',weight,'mu',mu,'df',8,'dimension',d,'globalChol',L, ...
    'rwScale',.5*2.38/sqrt(d),'chol',zeros(d,d,K),'logConstant',zeros(1,K));
for k=1:K
    proposal.chol(:,:,k)=chol(covariance(:,:,k),'lower');
    proposal.logConstant(k)=log(weight(k))+gammaln((8+d)/2)-gammaln(4)-d/2*log(8*pi) ...
        -sum(log(diag(proposal.chol(:,:,k))));
end
end

function z=mixture_draw(q)
k=find(rand<=cumsum(q.weight),1);
z=q.mu(:,k)+q.chol(:,:,k)*randn(q.dimension,1)/sqrt(2*randg(q.df/2)/q.df);
end

function value=mixture_logpdf(z,q)
terms=zeros(size(q.weight));
for k=1:numel(terms)
    v=q.chol(:,:,k)\(z-q.mu(:,k));
    terms(k)=q.logConstant(k)-(q.df+q.dimension)/2*log1p(sum(v.^2)/q.df);
end
maximum=max(terms); value=maximum+log(sum(exp(terms-maximum)));
end
