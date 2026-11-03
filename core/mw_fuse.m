function [Fit,Info]=mw_fuse(Windows,actions,Data,O) %#ok<INUSD>
%MW_FUSE Merge repeated window estimates without a prescribed source count.
% One report per window per group; complete-link diameter and vector-moment
% gates prevent long chains. At least one report must normally place the
% source inside the trusted window interior, to avoid boundary-only fits.
% A joint position/moment medoid, chosen from interior reports when present,
% represents each accepted group. Support counts distinct windows.
% Data is reserved for an optional later joint physical refinement. This
% version deliberately keeps identical, unrefined fusion for all methods.
O=defaults(O);t=tic;
P=zeros(0,3);M=P;wid=zeros(0,1);unary=wid;weights=wid;core=false(0,1);
assert(numel(actions)==numel(Windows));
for i=1:numel(Windows)
    assert(actions(i)>=1&&actions(i)<=numel(Windows(i).candidates));
    c=Windows(i).candidates(actions(i));
    if c.K==0||(isfield(c,'valid')&&~c.valid),continue;end
    good=all(isfinite(c.position),2)&all(isfinite(c.moment),2);
    p=c.position(good,:);m=c.moment(good,:);k=size(p,1);
    P=[P;p];M=[M;m];wid=[wid;repmat(i,k,1)]; %#ok<AGROW>
    unary=[unary;repmat(c.cost,k,1)]; %#ok<AGROW>
    weights=[weights;repmat(exp(-min(max(c.cost,0),30)),k,1)]; %#ok<AGROW>
    box=Windows(i).bounds_xy;inset=O.core_margin_fraction*(box(2,:)-box(1,:));
    lo=box(1,:)+inset;hi=box(2,:)-inset;
    inside=p(:,1)>=lo(1)&p(:,2)>=lo(2)&p(:,1)<=hi(1)&p(:,2)<=hi(2);
    core=[core;inside]; %#ok<AGROW>
end
n=size(P,1);labels=zeros(n,1);groups=cell(0,1);representatives=zeros(0,1);
% High-quality reports seed groups first; original report index breaks ties.
[~,visit]=sortrows([unary,(1:n).'],[1 2]);
for q=reshape(visit,1,[])
    best=0;best_distance=Inf;
    for g=1:numel(groups)
        members=groups{g};r=representatives(g);
        if any(wid(members)==wid(q)),continue;end
        d=sqrt(sum(bsxfun(@minus,P(members,:),P(q,:)).^2,2));
        if any(d>O.max_cluster_diameter_m),continue;end
        if norm(P(q,:)-P(r,:))>O.position_radius_m,continue;end
        dm=relative_moments(M(members,:),M(q,:));
        if any(dm>O.moment_relative_gate),continue;end
        score=(norm(P(q,:)-P(r,:))/O.position_radius_m)^2+ ...
            (relative_moments(M(r,:),M(q,:))/O.moment_relative_gate)^2;
        if score<best_distance,best=g;best_distance=score;end
    end
    if best==0
        groups{end+1,1}=q;representatives(end+1,1)=q; %#ok<AGROW>
    else
        groups{best}(end+1,1)=q;
        representatives(best)=medoid(groups{best},P,M,weights,core,O);
    end
end
template=struct('id',0,'count',0,'position',zeros(1,3), ...
    'moment',zeros(1,3),'window_ids',zeros(0,1),'report_ids',zeros(0,1), ...
    'core_window_ids',zeros(0,1),'core_count',0, ...
    'representative_report_id',0,'diameter_m',0,'mean_unary_cost',NaN);
clusters=repmat(template,0,1);rejected_groups=zeros(0,1);
rejected_no_core=0;rejected_insufficient_support=0;
for g=1:numel(groups)
    members=groups{g};support=numel(unique(wid(members)));
    core_windows=unique(wid(members(core(members))));
    if support<O.minimum_window_support
        rejected_groups(end+1,1)=g; %#ok<AGROW>
        rejected_insufficient_support=rejected_insufficient_support+1;
        continue;
    end
    if numel(core_windows)<O.minimum_core_window_support
        rejected_groups(end+1,1)=g; %#ok<AGROW>
        rejected_no_core=rejected_no_core+1;
        continue;
    end
    r=representatives(g);c=template;c.id=numel(clusters)+1;
    c.count=support;c.position=P(r,:);c.moment=M(r,:);
    c.window_ids=unique(wid(members));c.report_ids=members(:);
    c.core_window_ids=core_windows;c.core_count=numel(core_windows);
    c.representative_report_id=r;c.mean_unary_cost=mean(unary(members));
    for a=1:numel(members)
        delta=bsxfun(@minus,P(members,:),P(members(a),:));
        c.diameter_m=max(c.diameter_m,max(sqrt(sum(delta.^2,2))));
    end
    clusters(end+1,1)=c;labels(members)=c.id; %#ok<AGROW>
end
if isempty(clusters)
    positions=zeros(0,3);moments=zeros(0,3);support=zeros(0,1);support_core=zeros(0,1);
else
    positions=vertcat(clusters.position);moments=vertcat(clusters.moment);
    support=reshape([clusters.count],[],1);
    support_core=reshape([clusters.core_count],[],1);
end
Fit=struct('has_estimate',~isempty(clusters),'position',positions, ...
    'moment',moments,'support',support,'support_core',support_core, ...
    'all_positions',P,'all_moments',M, ...
    'window_ids',wid,'labels',labels,'clusters',clusters,'actions',actions);
Info=struct('fusion_seconds',toc(t),'selected_reports',n, ...
    'groups_before_support',numel(groups),'accepted_groups',numel(clusters), ...
    'rejected_groups',rejected_groups,'rejected_reports',sum(labels==0), ...
    'rejected_no_core',rejected_no_core, ...
    'rejected_insufficient_support',rejected_insufficient_support, ...
    'core_reports',core, ...
    'weights',weights,'unary_costs',unary,'refinement','none');
end

function O=defaults(O)
d=struct('position_radius_m',2,'moment_relative_gate',.75, ...
    'minimum_window_support',2,'refinement','none','max_cluster_diameter_m',4, ...
    'core_margin_fraction',.10,'minimum_core_window_support',1);
f=fieldnames(d);
for k=1:numel(f),if ~isfield(O,f{k}),O.(f{k})=d.(f{k});end,end
assert(O.position_radius_m>0&&O.moment_relative_gate>0);
assert(O.max_cluster_diameter_m>0);
assert(O.minimum_window_support>=1&& ...
    O.minimum_window_support==round(O.minimum_window_support));
assert(O.core_margin_fraction>=0&&O.core_margin_fraction<.5);
assert(O.minimum_core_window_support>=0&& ...
    O.minimum_core_window_support==round(O.minimum_core_window_support));
assert(strcmp(O.refinement,'none'),'This version supports refinement=''none'' only.');
end

function d=relative_moments(A,b)
scale=max(.5*(sqrt(sum(A.^2,2))+norm(b)),eps);
d=sqrt(sum(bsxfun(@minus,A,b).^2,2))./scale;
end

function report=medoid(ids,P,M,w,core,O)
eligible=ids;
if O.minimum_core_window_support>0&&any(core(ids)),eligible=ids(core(ids));end
cost=zeros(numel(eligible),1);v=w(ids);v=v/sum(v);
for k=1:numel(eligible)
    dp=sum(bsxfun(@minus,P(ids,:),P(eligible(k),:)).^2,2)/O.position_radius_m^2;
    dm=relative_moments(M(ids,:),M(eligible(k),:))/O.moment_relative_gate;
    cost(k)=sum(v.*(dp+dm.^2));
end
[~,q]=min(cost);report=eligible(q);
end
