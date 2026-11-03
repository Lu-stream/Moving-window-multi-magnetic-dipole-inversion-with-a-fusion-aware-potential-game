function [actions,D]=pg_min_actions(pool,scores,depth_prior)
% One min-residual branch per observation node; no residual rejection.
% 1=reject/no admissible branch, j+1=branch j. Shared application prior allowed.
actions=ones(pool.node_count,1);best_error=nan(pool.node_count,1);
for i=1:pool.node_count
 e=scores.error{i};
 if ~isempty(depth_prior)
  allow=pool.position{i}(:,3)>=depth_prior(1) & pool.position{i}(:,3)<=depth_prior(2);
  e(~allow)=Inf;
 end
 [value,j]=min(e);best_error(i)=value;
 if isfinite(value),actions(i)=j+1;end
end
D=struct('best_error',best_error,'kept',sum(actions>1));
end
