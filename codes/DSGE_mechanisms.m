function result=DSGE_mechanisms(root,number)
% Paired contractionary monetary policy IRFs, one mechanism off, no reestimation.
% Keep baseline posterior draws and recalibrate chi to the same steady state.
if nargin<1 || isempty(root), root=fileparts(mfilename('fullpath')); end
if nargin<2, number=2000; end
DSGE_setup;
previous=pwd; cleanup=onCleanup(@() cd(previous));
work=fullfile(root,'results','mechanism_work');
if ~isfolder(work), mkdir(work); end
addpath(work);
S=load(fullfile(root,'results','estimated_model.mat'));
pack=load(fullfile(root,'results','mcmc_diagnostics.mat'),'diagnostics');
bandsPack=load(fullfile(root,'results','irf_posterior_bands.mat'),'bands');
d=pack.diagnostics; pool=reshape(d.chains,[],size(d.chains,3));
drawIndices=bandsPack.bands.draw_indices;
number=min(number,numel(drawIndices)); drawIndices=drawIndices(1:number);

% Generate an isolated counterfactual model from the canonical baseline.
% util_switch=0 replaces the utilization FOC exactly by u=1.
source=fileread(fullfile(root,'medium_dsge.mod'));
tail=strfind(source,'@#if RUN_ESTIMATION'); assert(numel(tail)==1);
source=source(1:tail-1);
source=[source,sprintf('\nstoch_simul(order=1,irf=40,nograph,noprint) GDPobs Cobs Iobs N_obs wobs Piobs Robs Kobs uobs;\n')];
source=strrep(source,'parameters beta eta chi','parameters util_switch beta eta chi');
source=strrep(source,'load calibration_input.mat;',sprintf('util_switch=1;\nload calibration_input.mat;'));
old=' alpha*mc*Yw=u*(a1+a2*(u-1))*K(-1);';
assert(numel(strfind(source,old))==1);
source=strrep(source,old,' util_switch*(alpha*mc*Yw-u*(a1+a2*(u-1))*K(-1))+(1-util_switch)*(u-1)=0;');
fid=fopen(fullfile(work,'mechanism_model.mod'),'w','n','UTF-8');
assert(fid>=0); fprintf(fid,'%s',source); fclose(fid);
copyfile(fullfile(root,'results','dynare_work','calibration_input.mat'),work);
cd(work); evalin('base','dynare mechanism_model.mod noclearall');
M0=evalin('base','M_'); O0=evalin('base','oo_'); opts=evalin('base','options_');
opts.order=1; opts.nograph=true;
assert(isequal(M0.endo_names,S.M_.endo_names));
for j=1:numel(S.M_.param_names)
    M0.params(strcmp(M0.param_names,S.M_.param_names{j}))=S.M_.params(j);
end
vars={'GDPobs','Cobs','Iobs','N_obs','wobs','Piobs','Robs','Kobs','uobs'};
labels={'GDP','消费','投资','总劳动工时','真实工资','通胀率','名义利率','资本','资本利用率'};
cases={'baseline','no_habit','no_investment_cost','fixed_utilization','flexible_wages'};
caseLabels={'基准设定','关闭消费习惯','关闭投资调整成本','固定资本利用率','关闭工资粘性'};
indices=cellfun(@(name)find(strcmp(M0.endo_names,name)),vars);
paramNames=d.parameter_names; if ischar(paramNames), paramNames=cellstr(paramNames); end
parameterIndices=cellfun(@(name)find(strcmp(M0.param_names,name)),paramNames);
responses=zeros(40,9,5,number);
meanParameterResponses=zeros(40,9,5);
checks=zeros(number+1,5,3); % static residual, steady state difference, identities
eigenvalues=cell(5,1);
linearEquationErrors=zeros(5,1);

for b=0:number
    M=M0;
    if b==0, values=d.posterior_mean; else, values=pool(drawIndices(b),:)'; end
    M.params(parameterIndices)=values;
    baselineSS=[];
    for c=1:5
        Mc=M;
        if c==2, Mc.params(strcmp(Mc.param_names,'h'))=0; end
        if c==3, Mc.params(strcmp(Mc.param_names,'phiI'))=0; end
        if c==4, Mc.params(strcmp(Mc.param_names,'util_switch'))=0; end
        if c==5, Mc.params(strcmp(Mc.param_names,'tw'))=0; end
        [dr,info,Mc.params]=resol(0,Mc,opts,O0.dr,O0.steady_state, ...
            O0.exo_steady_state,O0.exo_det_steady_state);
        assert(info(1)==0,'Case %s, draw %d failed to solve: %d',cases{c},b,info(1));
        if b==0, eigenvalues{c}=dr.eigval; end
        impulse=zeros(Mc.exo_nbr,1); impulse(strcmp(Mc.exo_names,'eR'))=1;
        path=irf(Mc,opts,dr,impulse,40,0,1,1);
        if b==0
            linearEquationErrors(c)=verify_equilibrium(Mc,dr,path,impulse);
            assert(linearEquationErrors(c)<1e-6,'Linearized equilibrium residual failed.');
        end
        current=path(indices,:)'; assert(all(isfinite(current),'all'));
        if b==0, meanParameterResponses(:,:,c)=current; else, responses(:,:,c,b)=current; end
        residual=feval('mechanism_model.static_resid',dr.ys,zeros(1,Mc.exo_nbr),Mc.params);
        invariant={'GDP','Y','C','I','K','N','L','w','wh','Pi','R','u','A','G'};
        invariantIdx=cellfun(@(name)find(strcmp(Mc.endo_names,name)),invariant);
        if c==1, baselineSS=dr.ys(invariantIdx); end
        checks(b+1,c,1)=max(abs(residual));
        checks(b+1,c,2)=max(abs(dr.ys(invariantIdx)-baselineSS));
        checks(b+1,c,3)=verify_path(Mc,dr,path,c);
        assert(max(checks(b+1,c,:))<1e-7,'Counterfactual identity/steady state check failed.');
    end
    if mod(b,200)==0, fprintf('Mechanism comparisons: %d/%d draws\n',b,number); end
end
posteriorMean=mean(responses,4);
baselineError=NaN;
if number==bandsPack.bands.number_of_draws
    policyIndex=find(strcmp(S.M_.exo_names,'eR'));
    baselineError=max(abs(posteriorMean(:,:,1)-bandsPack.bands.posterior_mean(:,:,policyIndex)),[],'all');
    assert(baselineError<1e-8,'Baseline IRFs differ from existing posterior mean IRFs.');
end
result=struct('variables',{vars},'case_names',{cases},'case_labels',{caseLabels}, ...
    'posterior_mean',posteriorMean,'paired_responses',responses, ...
    'irfs_at_posterior_mean_parameters',meanParameterResponses, ...
    'draw_indices',drawIndices,'number_of_draws',number,'checks',checks, ...
    'mean_parameter_eigenvalues',{eigenvalues},'baseline_crosscheck_error',baselineError, ...
    'parameter_names',{paramNames},'parameter_draws',pool(drawIndices,:), ...
    'parameter_mean',d.posterior_mean,'generated_at',char(datetime('now')));
result.mean_parameter_linear_equation_errors=linearEquationErrors;
result.shock_name='eR';
result.shock_label='一标准差正向（紧缩性）货币政策冲击';
% Pilot calls check the full workflow without replacing published outputs.
if number<bandsPack.bands.number_of_draws
    save(fullfile(work,'pilot_results.mat'),'result','-v7.3');
    disp(squeeze(posteriorMean(1,:,:))'); return
end
save(fullfile(root,'results','mechanism_comparison.mat'),'result','-v7.3');
% Long format preserves all horizons and variables for reproduction.
[tt,jj,cc]=ndgrid(0:39,1:9,1:5);
tableLong=table(tt(:),reshape(string(vars(jj(:))),[],1),reshape(string(cases(cc(:))),[],1),posteriorMean(:), ...
    'VariableNames',{'Quarter','Variable','Specification','PosteriorMeanIRF'});
writetable(tableLong,fullfile(root,'results','mechanism_comparison.csv'));
disp(squeeze(posteriorMean(1,:,:))');
fprintf('PASS: %d paired draws, all five models solved; baseline max difference %.3g.\n',number,baselineError);
end

function error=verify_equilibrium(M,dr,path,impulse)
% Independently verify every dynamic equation along the deterministic IRF.
% Dynare 7 uses a full lag/current/lead stack for dynamic_resid.
step=1e-4; error=0;
for t=1:size(path,2)-1
    previous=zeros(M.endo_nbr,1); if t>1, previous=path(:,t-1); end
    direction=[previous;path(:,t);path(:,t+1)];
    exo=zeros(M.exo_nbr,1); if t==1, exo=impulse; end
    yp=repmat(dr.ys,3,1)+step*direction;
    ym=repmat(dr.ys,3,1)-step*direction;
    plus=feval('mechanism_model.dynamic_resid',yp,step*exo,M.params,dr.ys);
    minus=feval('mechanism_model.dynamic_resid',ym,-step*exo,M.params,dr.ys);
    error=max(error,max(abs((plus-minus)/(2*step))));
end
end

function error=verify_path(M,dr,path,c)
idx=@(name)find(strcmp(M.endo_names,name));
p=@(name)M.params(strcmp(M.param_names,name));
r=@(name)path(idx(name),:); ss=@(name)dr.ys(idx(name));
errors=[r('GDP')-r('C')-r('I')-r('G'); ...
    r('Y')-r('GDP')-p('a1')*ss('K')*r('u'); ...
    r('GDPobs')-100*r('GDP')/ss('GDP')];
% Align the lag explicitly, including the zero initial deviation.
policy=r('R')/ss('R')-p('rhoR')*[0,path(idx('R'),1:end-1)]/ss('R') ...
    -(1-p('rhoR'))*(p('phipi')*r('Pi')/ss('Pi')+p('phiy')*r('GDP')/ss('GDP'));
policy(1)=policy(1)-p('sigR'); % Unit positive eR innovation on impact.
errors=[errors;policy];
if c==3, errors=[errors;r('q')]; end % With monetary policy only, zeta=1 and q=1.
if c==4, errors=[errors;r('u')]; end
if c==5, errors=[errors;r('w')/ss('w')-r('wh')/ss('wh')]; end
error=max(abs(errors),[],'all');
end

