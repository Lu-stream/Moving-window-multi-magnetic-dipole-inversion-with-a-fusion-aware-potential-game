function fit=cluster_primin_actions(pool,actions,opt)
%CLUSTER_PRIMIN_ACTIONS Identical output stage for every comparison method.
% One retained candidate per observation. No source truth is used here.
    use=find(actions>1); p=zeros(numel(use),3); m=p;
    for k=1:numel(use)
        i=use(k); a=actions(i)-1; p(k,:)=pool.position{i}(a,:); m(k,:)=pool.moment{i}(a,:);
    end
    assert(any(strcmp(opt.cluster_space,{'xy','xyz'})),'cluster_space must be xy or xyz.');
    if strcmp(opt.cluster_space,'xy'), coords=p(:,1:2); else, coords=p; end
    labels=local_dbscan(coords,opt.cluster_epsilon_m,opt.cluster_min_points);
    ids=unique(labels(labels>0));
    template=struct('id',0,'count',0,'position',nan(1,3),'moment',nan(1,3),'depth_iqr_m',NaN);
    clusters=repmat(template,0,1);
    for k=1:numel(ids)
        selected=labels==ids(k); c=template; c.id=ids(k); c.count=sum(selected);
        c.position=median(p(selected,:),1); c.moment=median(m(selected,:),1);
        z=sort(p(selected,3)); c.depth_iqr_m=percentile(z,.75)-percentile(z,.25);
        clusters(end+1,1)=c;
    end
    if isempty(clusters), primary=template; else
        [~,ix]=sort([clusters.count],'descend'); clusters=clusters(ix); primary=clusters(1);
    end
    fit=struct('has_estimate',~isempty(clusters),'position',primary.position, ...
        'moment',primary.moment,'clusters',clusters,'points',p,'point_moments',m, ...
        'node_ids',use,'labels',labels,'primary_cluster_id',primary.id,'actions',actions);
end

function labels=local_dbscan(P,epsilon,min_points)
% Exact DBSCAN, same input/neighbor order as the old dense implementation.
% Spatial bins bound memory; no dense n-by-n distance matrix is allocated.
assert(epsilon>0 && min_points>=1 && min_points==round(min_points));
n=size(P,1);labels=zeros(n,1);visited=false(n,1);cluster=0;
if n==0,return;end
dim=size(P,2);bins=floor(P/epsilon);[occupied,~,group]=unique(bins,'rows');
nb=size(occupied,1);members=accumarray(group,(1:n).',[],@(x){x});
lookup=containers.Map('KeyType','char','ValueType','double');
for k=1:nb,lookup(bin_key(occupied(k,:)))=k;end
if dim==2,[a,b]=ndgrid(-1:1);offset=[a(:),b(:)];
else,[a,b,c]=ndgrid(-1:1);offset=[a(:),b(:),c(:)];end
near_cells=cell(nb,1);
for k=1:nb
 ids=zeros(0,1);
 for t=1:size(offset,1)
  key=bin_key(occupied(k,:)+offset(t,:));
  if isKey(lookup,key),ids=[ids;members{lookup(key)}];end %#ok<AGROW>
 end
 near_cells{k}=sort(ids); % preserves dense DBSCAN's ascending neighbor order
end
for i=1:n
 if visited(i),continue;end
 visited(i)=true;seed=neighbors(i,P,group,near_cells,epsilon);
 if numel(seed)<min_points,labels(i)=-1;continue;end
 cluster=cluster+1;labels(i)=cluster;
 queue=seed(:);in_queue=false(n,1);in_queue(queue)=true;head=1;
 while head<=numel(queue)
  j=queue(head);head=head+1;
  if ~visited(j)
   visited(j)=true;near=neighbors(j,P,group,near_cells,epsilon);
   if numel(near)>=min_points
    new=near(~in_queue(near));queue=[queue;new(:)];in_queue(new)=true; %#ok<AGROW>
   end
  end
  if labels(j)<=0,labels(j)=cluster;end
 end
end
end
function key=bin_key(v)
v(v==0)=0; % canonicalize signed zero before string keys
key=sprintf('%.0f,',v);
end
function ids=neighbors(i,P,group,near_cells,epsilon)
ids=near_cells{group(i)};d=bsxfun(@minus,P(ids,:),P(i,:));
ids=ids(sum(d.^2,2)<=epsilon^2);
end

function q=percentile(v,p)
    z=1+(numel(v)-1)*p; a=floor(z); b=ceil(z); q=v(a)+(z-a)*(v(b)-v(a));
end
