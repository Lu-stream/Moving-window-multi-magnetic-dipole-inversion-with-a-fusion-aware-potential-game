function [Fit,Info,Raw]=mw_primin_baseline(Data,Pool,O)
% Explicit baseline: raw Primin equations, full-map min residual, DBSCAN.
% No rejection, player interactions, continuous inversion or true source count.
t=tic;
if Pool.branch_count>0
 Raw=struct('positions',vertcat(Pool.position{:}),'moments',vertcat(Pool.moment{:}));
else,Raw=struct('positions',zeros(0,3),'moments',zeros(0,3));end
Scores=pg_scores(Pool,Data.Bz,Data.X,Data.Y,O,'G');
actions=pg_min_actions(Pool,Scores,[]);
F=cluster_primin_actions(Pool,actions,O);
if isempty(F.clusters),p=zeros(0,3);m=p;support=zeros(0,1);
else,p=vertcat(F.clusters.position);m=vertcat(F.clusters.moment);support=[F.clusters.count].';end
Fit=struct('position',p,'moment',m,'support',support,'all_positions',F.points,...
 'all_moments',F.point_moments,'labels',F.labels,'clusters',F.clusters,'actions',actions);
Info=struct('elapsed_seconds',toc(t),'raw_branches',Pool.branch_count,...
 'retained_points',numel(F.node_ids),'score_samples',Scores.sample_count,...
 'definition','Primin analytic candidates + full-map minimum residual + unchanged DBSCAN');
end
