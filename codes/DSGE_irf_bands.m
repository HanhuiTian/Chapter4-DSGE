function bands=DSGE_irf_bands(S,d,variables,number)
% Parameter uncertainty: recompute the full model at retained posterior draws.
[n,m,p]=size(d.chains);
pool=reshape(d.chains,n*m,p);
state=rng; cleanup=onCleanup(@() rng(state)); rng(20261111,'twister');
drawIndices=randperm(n*m,min(number,n*m)); number=numel(drawIndices);
response=zeros(40,numel(variables),S.M_.exo_nbr,number);
indices=zeros(numel(variables),1);
for j=1:numel(variables), indices(j)=find(strcmp(S.M_.endo_names,variables{j})); end
options=S.options_; options.order=1;
for b=1:number
    M=set_all_parameters(pool(drawIndices(b),:)',S.estim_params_,S.M_);
    [dr,info,M.params]=resol(0,M,options,S.oo_.dr,S.oo_.steady_state,S.oo_.exo_steady_state,S.oo_.exo_det_steady_state);
    assert(info(1)==0,'Retained posterior draw failed to solve: %d',info(1));
    for k=1:M.exo_nbr
        impulse=zeros(M.exo_nbr,1); impulse(k)=1;
        path=irf(M,options,dr,impulse,40,0,1,1);
        response(:,:,k,b)=path(indices,:)';
    end
end
ordered=sort(response,4); span=floor(.9*number);
width=ordered(:,:,:,span+1:end)-ordered(:,:,:,1:end-span);
[~,location]=min(width,[],4);
low=zeros(size(location)); high=low;
for t=1:40
    for j=1:numel(variables)
        for k=1:S.M_.exo_nbr
            index=location(t,j,k);
            low(t,j,k)=ordered(t,j,k,index);
            high(t,j,k)=ordered(t,j,k,index+span);
        end
    end
end
bands=struct('low',low,'high',high,'posterior_mean',mean(response,4), ...
    'posterior_median',median(response,4), ...
    'draw_indices',drawIndices,'number_of_draws',number,'probability',.90);
end
