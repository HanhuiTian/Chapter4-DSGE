function DSGE_slides_figures(selection,root)
% Redraw Chapter 4 evidence for 4:3 Beamer slides, using saved baseline results.
% Pass 'prior_families' to redraw only the prior distribution examples.
% Pass 'prior_posterior' to redraw only the three prior/posterior comparison slides.
% Pass 'prior_posterior_4' to redraw only the fourth comparison figure.
% Pass 'mechanisms' to redraw only the four mechanism comparison slides.
if nargin<1, selection='all'; end
prior_only=strcmp(selection,'prior_families');
posterior_only=ismember(selection,{'prior_posterior','prior_posterior_4'});
posterior_fourth=strcmp(selection,'prior_posterior_4');
mechanisms_only=strcmp(selection,'mechanisms');
assert(nargin>=2,'Use DSGE_main(''figures'') to select the output workspace.');
DSGE_setup;
out=fullfile(root,'figures'); if ~isfolder(out), mkdir(out); end
S=load(fullfile(root,'results','estimated_model.mat'));
blue=[.12 .32 .62];
if ~prior_only
if ~posterior_only
T=load(fullfile(root,'results','analysis_results.mat')); a=T.analysis;
T=load(fullfile(root,'results','mechanism_comparison.mat')); mechanisms=T.result;
D=load(fullfile(root,'data','dsge_observables.mat'));
blue=[.12 .32 .62]; orange=[.82 .32 .12];
colors=[.18 .43 .67;.85 .38 .15;.94 .70 .16;.45 .65 .32;.53 .40 .69;.25 .67 .72;.70 .34 .52;.72 .72 .72];
labels={'GDP','消费','投资','劳动工时','真实工资','通胀率','名义利率','资本','资本利用率'};
shockLabels={'TFP','外生支出','货币政策','风险溢价','投资效率','价格加成','工资加成'};
primary=[1 2 3 4 6 7]; factors=[5 8 9]; time=D.year_quarter(:,1)+(D.year_quarter(:,2)-1)/4;

% All nine reported IRF variables: six macro variables plus three factor variables.
if ~mechanisms_only
for k=1:7
    [f,tl]=canvas(2,3);
    for jj=1:6
        j=primary(jj); ax=nexttile(tl); draw_irf(ax,a,j,k,blue);
        title(ax,labels{j});
    end
    output(f,out,['irf_',a.shock_names{k}]);
end
for group=1:4
    shocks=(2*group-1):min(2*group,7); [f,tl]=canvas(numel(shocks),3);
    for k=shocks
        for j=factors
            ax=nexttile(tl); draw_irf(ax,a,j,k,blue);
            title(ax,[shockLabels{k},'：',labels{j}]);
        end
    end
    output(f,out,sprintf('irf_factors_%d',group));
end
end
selection={[1 2 3 4 6 7],[1 2 3 4 8 9],[1 2 3 4 8 9],[1 2 3 4 5 6]};
for c=2:5
    [f,tl]=canvas(2,3);
    for j=selection{c-1}
        ax=nexttile(tl); hold(ax,'on');
        h1=plot(ax,0:39,mechanisms.posterior_mean(:,j,1),'Color',blue,'LineWidth',2.8);
        h2=plot(ax,0:39,mechanisms.posterior_mean(:,j,c),'--','Color',orange,'LineWidth',2.8);
        yline(ax,0,':','Color',[.5 .5 .5]); title(ax,labels{j}); axes_irf(ax);
        if j==selection{c-1}(1)
            lg=legend(ax,[h1,h2],{'基准设定',mechanisms.case_labels{c}},'Orientation','horizontal','Box','off');
            lg.Layout.Tile='south';
        end
    end
    output(f,out,['mechanism_',mechanisms.case_names{c}]);
end
if mechanisms_only, return; end
% Observations use two 2x2 figures, keeping source units explicit.
dataFields={'GDPobs','Cobs','Iobs','Nobs','wobs','Piobs','Robs'};
dataLabels={'GDP','消费','投资','劳动工时','真实每小时劳动报酬','通胀率','名义利率'};
for group=1:2
    [f,tl]=canvas(2,2);
    for j=(group-1)*4+1:min(group*4,7)
        ax=nexttile(tl); plot(ax,time,D.(dataFields{j}),'Color',blue,'LineWidth',2.1);
        yline(ax,0,':','Color',[.5 .5 .5]); title(ax,dataLabels{j});
        set(ax,'XLim',[1984 2007.3],'XTick',[1984 1990 1996 2002 2007]);
        if j<=5, ylabel(ax,'对数周期项（%）'); else, ylabel(ax,'季度对数百分点'); end
        axes_style(ax);
    end
    output(f,out,sprintf('observations_%d',group));
end

end
T=load(fullfile(root,'results','mcmc_diagnostics.mat')); d=T.diagnostics;
% Prior curves use Dynare's actual hyperparameters; posterior KDE uses 20000
% evenly spaced retained draws, with support-aware boundary correction.
pool=reshape(d.chains,[],size(d.chains,3));
keep=round(linspace(1,size(pool,1),min(20000,size(pool,1))));
symbols={'h','\phi_I','a_2','\theta_p','\theta_w','\rho_R','\phi_\pi','\phi_y', ...
    '\rho_A','\rho_G','\rho_b','\rho_\zeta','\rho_p','\rho_w', ...
    '\sigma_A','\sigma_G','\sigma_R','\sigma_b','\sigma_\zeta','\sigma_p','\sigma_w'};
posterior_groups=1:4;
if posterior_only, posterior_groups=1:3; end
if posterior_fourth, posterior_groups=4; end
for group=posterior_groups
    inds=(group-1)*6+1:min(group*6,21); rows=ceil(numel(inds)/3);
    [f,tl]=canvas(rows,3);
    for j=inds
        ax=nexttile(tl); sample=pool(keep,j); b=S.bayestopt_;
        [~,~,priorX]=draw_prior_density(j,b);
        if b.pshape(j)==1
            lower=0.0001; upper=.9999; support=[0 1];
        else
            lower=max(1e-8,min(priorX(1),quantile(sample,.001))*.5);
            upper=max(priorX(end),quantile(sample,.999)*1.08); support='positive';
        end
        grid=linspace(lower,upper,450); prior=density(j,grid,b);
        post=ksdensity(sample,grid,'Support',support,'BoundaryCorrection','reflection');
        hold(ax,'on'); h1=plot(ax,grid,prior,'--','Color',[.5 .5 .5],'LineWidth',2.5);
        h2=plot(ax,grid,post,'Color',blue,'LineWidth',2.6);
        title(ax,['$',symbols{j},'$'],'Interpreter','latex','FontSize',22);
        xlabel(ax,'参数取值'); ylabel(ax,'概率密度'); axes_style(ax); xlim(ax,[lower upper]);
        if j==inds(1)
            legend_labels={'先验分布','后验分布'};
            lg=legend(ax,[h1,h2],legend_labels,'Orientation','horizontal','Box','off'); lg.Layout.Tile='south';
        end
    end
    output(f,out,sprintf('prior_posterior_%d',group));
end
if posterior_only, return; end
end
[f,tl]=canvas(2,2);
demoNames={'Beta(0.6, 0.15)','Gamma(4, 2)','IG1(0.01, 0.005)','N(0, 1)'};
for j=1:4
    ax=nexttile(tl);
    if j==1, grid=linspace(.001,.999,400); curve=density(1,grid,S.bayestopt_);
    elseif j==2, grid=linspace(.001,12,400); curve=density(2,grid,S.bayestopt_);
    elseif j==3, grid=linspace(.001,.04,400); curve=density(15,grid,S.bayestopt_);
    else, grid=linspace(-4,4,400); curve=exp(-grid.^2/2)/sqrt(2*pi); end
    plot(ax,grid,curve,'Color',blue,'LineWidth',2.8); title(ax,demoNames{j});
    xlabel(ax,'参数取值'); ylabel(ax,'概率密度'); axes_style(ax);
end
output(f,out,'prior_families');
if prior_only, return; end

% Variance decomposition and history use the same colors for the seven shocks.
[f,tl]=canvas(1,1); ax=nexttile(tl);
bars=barh(ax,a.variance_decomposition.shares_percent,'stacked','EdgeColor','none');
for k=1:7, bars(k).FaceColor=colors(k,:); end
set(ax,'YTick',1:6,'YTickLabel',labels(primary),'YDir','reverse','XLim',[0 100]);
xlabel(ax,'无条件方差贡献（%）'); axes_style(ax);
lg=legend(ax,bars,shockLabels,'NumColumns',4,'Location','southoutside','Box','off');
lg.FontSize=15; output(f,out,'variance_decomposition');
for group=1:3
    [f,tl]=canvas(2,1);
    for j=(2*group-1):2*group
        ax=nexttile(tl); hold(ax,'on'); pieces=a.historical_contributions(:,:,j);
        bp=bar(ax,time,max(pieces,0),1,'stacked','EdgeColor','none');
        bn=bar(ax,time,min(pieces,0),1,'stacked','EdgeColor','none');
        for k=1:8, bp(k).FaceColor=colors(k,:); bn(k).FaceColor=colors(k,:); end
        obs=plot(ax,time,a.observed(:,j),'k','LineWidth',2.5);
        title(ax,labels{primary(j)}); axes_style(ax);
        set(ax,'XLim',[1983.8 2007.4],'XTick',[1984 1990 1996 2002 2007]);
        if j<=4, ylabel(ax,'对数周期项（%）'); else, ylabel(ax,'季度对数百分点'); end
        if j==2*group-1
            lg=legend(ax,[bp,obs],[shockLabels,{'初始状态','观测序列'}], ...
                'NumColumns',5,'Box','off','Orientation','horizontal'); lg.Layout.Tile='south';
        end
    end
    output(f,out,sprintf('history_%d',group));
end
fprintf('Chapter 4 figures regenerated in %s\n',out);
end

function curve=density(j,x,b)
switch b.pshape(j)
    case 1, curve=betapdf((x-b.p3(j))/(b.p4(j)-b.p3(j)),b.p6(j),b.p7(j))/(b.p4(j)-b.p3(j));
    case 2, curve=gampdf(x-b.p3(j),b.p6(j),b.p7(j));
    case 3, curve=normpdf(x,b.p6(j),b.p7(j));
    case 4, curve=exp(lpdfig1(x-b.p3(j),b.p6(j),b.p7(j)));
    otherwise, error('Unsupported prior family');
end
end
function [f,tl]=canvas(rows,cols)
if rows==1 && cols>1, h=470; else, h=900; end
f=figure('Visible','off','Color','w','Position',[30 30 1600 h]);
tl=tiledlayout(f,rows,cols,'TileSpacing','compact','Padding','compact');
set(f,'DefaultAxesFontName','Microsoft YaHei','DefaultAxesFontSize',15, ...
    'DefaultTextFontName','Microsoft YaHei','DefaultTextFontSize',17);
end
function axes_style(ax)
set(ax,'FontName','Microsoft YaHei','FontSize',15,'Box','off','YGrid','on','GridAlpha',.13);
ax.Title.FontSize=19; ax.XLabel.FontSize=16; ax.YLabel.FontSize=16;
end
function axes_irf(ax)
xlabel(ax,'季度'); ylabel(ax,'偏离稳态（%）'); axes_style(ax);
set(ax,'XLim',[0 39],'XTick',[0 10 20 30]);
end
function draw_irf(ax,a,j,k,blue)
hold(ax,'on'); low=a.irf_bands.low(:,j,k); high=a.irf_bands.high(:,j,k);
fill(ax,[0:39,39:-1:0],[low',flipud(high)'],[.79 .86 .95],'EdgeColor','none');
plot(ax,0:39,a.irfs(:,j,k),'Color',blue,'LineWidth',2.7);
yline(ax,0,':','Color',[.5 .5 .5]); axes_irf(ax);
end
function output(f,out,name)
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',200); close(f);
end
