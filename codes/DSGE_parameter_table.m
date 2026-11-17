function summary=DSGE_parameter_table(root,S,d)
% The parameter table in the slides reports the joint posterior mode.
package=fileparts(fileparts(mfilename('fullpath')));
[~,source]=fileparts(root);
if strcmp(source,'published')
    file=fullfile(package,'estimates','posterior_mode.mat');
else
    file=fullfile(root,'results','dynare_work','medium_dsge','Output','medium_dsge_mode.mat');
end
mode=load(file,'xparam1');
summary=table(string(d.parameter_names(:)),S.bayestopt_.pshape(:), ...
    S.bayestopt_.p1(:),S.bayestopt_.p2(:),mode.xparam1(:),d.posterior_mean(:), ...
    d.hpd90(:,1),d.hpd90(:,2),'VariableNames', ...
    {'Parameter','PriorFamilyCode','PriorMean','PriorSD','PosteriorMode','PosteriorMean','HPD90Low','HPD90High'});
writetable(summary,fullfile(root,'results','parameter_table.csv'));
end
