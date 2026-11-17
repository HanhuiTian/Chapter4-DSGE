function result=DSGE_mcmc_diagnostics(work, parameterNames, burnFraction)
% Rank normalized split Rhat and bulk/tail ESS, without extra toolboxes.
% Reference: Vehtari et al. (2021), Bayesian Analysis, 16(2), 667-718.
% The input chains include repeated states after rejected proposals.
folder=fullfile(work,'medium_dsge','metropolis');
files=dir(fullfile(folder,'medium_dsge_mh*_blck*.mat'));
records=[];
for i=1:numel(files)
    tok=regexp(files(i).name,'_mh(\d+)_blck(\d+)\.mat$','tokens','once');
    if ~isempty(tok), records(end+1,:)=[str2double(tok{2}),str2double(tok{1}),i]; end %#ok<AGROW>
end
assert(~isempty(records),'No completed MCMC chains found.');
records=sortrows(records,[1,2]); blocks=unique(records(:,1));
chains=cell(numel(blocks),1); acceptance=zeros(numel(blocks),1);
rawLength=zeros(numel(blocks),1);
for b=1:numel(blocks)
    subset=records(records(:,1)==blocks(b),:); x=[];
    for k=1:size(subset,1)
        part=load(fullfile(folder,files(subset(k,3)).name),'x2');
        x=[x;part.x2]; %#ok<AGROW>
    end
    rawLength(b)=size(x,1);
    acceptance(b)=mean(any(diff(x,1,1)~=0,2));
    chains{b}=x(floor(size(x,1)*burnFraction)+1:end,:);
end
assert(numel(unique(rawLength))==1,'Chains have unequal completed lengths.');
n=size(chains{1},1); p=size(chains{1},2); m=numel(chains);
assert(p==numel(parameterNames));
allChains=zeros(n,m,p);
for b=1:m, allChains(:,b,:)=reshape(chains{b},n,1,p); end
pooled=vertcat(chains{:});
rh=zeros(p,1); bulk=zeros(p,1); tail=zeros(p,1); meanESS=zeros(p,1); interval=zeros(p,2);
for j=1:p
    x=allChains(:,:,j); half=floor(n/2);
    split=[x(1:half,:),x(end-half+1:end,:)];
    ranked=rank_normalize(split);
    folded=rank_normalize(abs(split-median(split(:))));
    rh(j)=max(classic_rhat(ranked),classic_rhat(folded));
    bulk(j)=effective_size(ranked);
    meanESS(j)=effective_size(split);
    sorted=sort(pooled(:,j)); count=numel(sorted);
    lo=sorted(max(1,ceil(.05*count))); hi=sorted(min(count,ceil(.95*count)));
    tail(j)=min(effective_size(double(split<=lo)),effective_size(double(split>=hi)));
    width=floor(.9*count);
    [~,ind]=min(sorted(width+1:end)-sorted(1:end-width));
    interval(j,:)=[sorted(ind),sorted(ind+width)];
end
result.parameter_names=parameterNames(:);
result.posterior_mean=mean(pooled,1)';
result.posterior_median=median(pooled,1)';
result.hpd90=interval;
result.rhat=rh; result.bulk_ess=bulk; result.tail_ess=tail;
result.mean_ess=meanESS;
result.mcse_mean=std(pooled,0,1)'./sqrt(meanESS);
result.acceptance=acceptance; result.raw_draws_per_chain=rawLength;
result.retained_draws=n*m; result.burn_fraction=burnFraction;
result.passed=all(rh<1.01 & bulk>=400 & tail>=400);
result.chains=allChains;
result.table=table(parameterNames(:),result.posterior_mean,interval(:,1),interval(:,2), ...
    rh,bulk,tail,result.mcse_mean,'VariableNames', ...
    {'Parameter','PosteriorMean','HPD90Low','HPD90High','Rhat','BulkESS','TailESS','MCSE'});
end

function z=rank_normalize(x)
[values,order]=sort(x(:)); n=numel(values);
starts=[1;find(diff(values)~=0)+1]; ends=[starts(2:end)-1;n];
ranks=zeros(n,1);
for k=1:numel(starts), ranks(starts(k):ends(k))=(starts(k)+ends(k))/2; end
original=zeros(n,1); original(order)=ranks;
z=reshape(sqrt(2)*erfinv(2*(original-3/8)/(n+1/4)-1),size(x));
end

function value=classic_rhat(x)
n=size(x,1); W=mean(var(x,0,1)); B=n*var(mean(x,1),0,2);
value=sqrt(((n-1)/n*W+B/n)/W);
end

function value=effective_size(x)
[n,m]=size(x); centered=x-mean(x,1);
fftLength=2^nextpow2(2*n);
transform=fft(centered,fftLength,1);
acov=real(ifft(abs(transform).^2,[],1)); acov=acov(1:n,:)/n;
W=mean(var(x,0,1)); B=n*var(mean(x,1),0,2);
variancePlus=(n-1)/n*W+B/n;
if variancePlus<=0, value=0; return; end
rho=1-(W-mean(acov,2))/variancePlus; rho(1)=1;
pairs=rho(1:2:end-1)+rho(2:2:end);
stop=find(pairs<=0,1); if ~isempty(stop), pairs=pairs(1:stop-1); end
if isempty(pairs), value=n*m; return; end
for k=2:numel(pairs), pairs(k)=min(pairs(k),pairs(k-1)); end
tau=max(1/log10(n*m),-1+2*sum(pairs));
value=min(n*m,n*m/tau);
end
