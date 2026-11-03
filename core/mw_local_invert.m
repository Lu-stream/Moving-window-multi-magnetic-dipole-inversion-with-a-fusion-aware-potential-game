function [Windows,Info]=mw_local_invert(Data,Pool,O)
%MW_LOCAL_INVERT Overlapping-window one/two-dipole variable-projection fits.
% Only measured Data and analytic Pool candidates enter this function.
% Positive depth is below the sensor plane; moments use z positive upward.
% Action 1 means reject the window and has a fixed cost of 1. All other
% actions jointly fit the SUM of K dipole fields plus a local background.
% fminsearch is used, so Optimization/Statistics Toolboxes are not required.
started=tic;
O=defaults(O);
assert(isequal(size(Data.Bz),size(Data.X),size(Data.Y)),...
    'Bz, X and Y must have identical shapes.');
assert(O.max_sources>=1 && O.max_sources<=2,'max_sources must be 1 or 2.');
r=O.half_width_pixels;
assert(r>=2 && r==floor(r),'half_width_pixels must be an integer >= 2.');
[ny,nx]=size(Data.Bz);
assert(ny>=2*r+1 && nx>=2*r+1,'The map is smaller than a full window.');
rr=(1+r):O.stride_pixels:(ny-r);cc=(1+r):O.stride_pixels:(nx-r);
% Optional data-only boundary coverage for rectangular QDM maps.
if isfield(O,'include_last_window') && O.include_last_window
    rr=unique([rr ny-r]);cc=unique([cc nx-r]);
end
global_rms=sqrt(mean(Data.Bz(:).^2));
field_floor=max(O.relative_field_floor*global_rms,eps(max(global_rms,1)));
template=struct('id',0,'center_xy',[0 0],'bounds_xy',zeros(2,2),...
    'linear_indices',zeros(0,1),'sigma_scale',0,'candidates',struct([]));
Windows=repmat(template,numel(rr)*numel(cc),1);
Info=struct('elapsed_seconds',0,'windows_total',numel(Windows),...
    'windows_fitted',0,'fit_attempts',0,'objective_evaluations',0,...
    'accepted_hypotheses',0,'accepted_nonconverged_hypotheses',0,...
    'field_floor_nT',field_floor,...
    'residual_reference',O.residual_reference,...
    'minimum_component_fraction',O.minimum_component_fraction);
w=0;
for ir=1:numel(rr)
    for ic=1:numel(cc)
        w=w+1;
        [jj,ii]=meshgrid((cc(ic)-r):(cc(ic)+r),(rr(ir)-r):(rr(ir)+r));
        ids=sub2ind([ny,nx],ii(:),jj(:));
        x=Data.X(ids);y=Data.Y(ids);b=Data.Bz(ids);
        bounds=[min(x),min(y);max(x),max(y)];
        center=[Data.X(rr(ir),cc(ic)),Data.Y(rr(ir),cc(ic))];
        scale=max(sqrt(mean(b.^2)),field_floor);
        H=background_matrix(x,y,center,bounds,O.background_order);
        reject=blank_action();reject.K=0;reject.background=H\b;
        reject.prediction=H*reject.background;
        reject.relative_residual=sqrt(mean((reject.prediction-b).^2))/scale;
        reject.cost=1;reject.condition=1;reject.valid=true;
        Windows(w)=template;Windows(w).id=w;
        Windows(w).center_xy=center;Windows(w).bounds_xy=bounds;
        Windows(w).linear_indices=ids;Windows(w).sigma_scale=scale;
        Windows(w).candidates=reject;
        % Eligibility uses the anomaly after the chosen background removal;
        % normalization deliberately remains max(total window RMS, floor).
        anomaly=b-reject.prediction;
        if sqrt(mean(anomaly.^2))<=O.signal_fraction_floor*global_rms
            continue;
        end
        lo=max(O.position_bounds(1,:),[bounds(1,:),O.position_bounds(1,3)]);
        hi=min(O.position_bounds(2,:),[bounds(2,:),O.position_bounds(2,3)]);
        if any(hi<=lo),continue;end
        Info.windows_fitted=Info.windows_fitted+1;
        seeds=make_seeds(Data,Pool,ids,x,y,b,H,scale,lo,hi,O);
        singles=repmat(blank_action(),0,1);
        nstarts=min(O.starts_per_order,size(seeds,1));
        for s=1:nstarts
            [a,ncall]=optimize(seeds(s,:),x,y,b,H,scale,lo,hi,O);
            Info.fit_attempts=Info.fit_attempts+1;
            Info.objective_evaluations=Info.objective_evaluations+ncall;
            singles=append_distinct(singles,a,O);
        end
        singles=sort_actions(singles,O.max_hypotheses_per_order);
        for s=1:numel(singles)
            Windows(w).candidates(end+1,1)=singles(s); %#ok<AGROW>
        end
        if O.max_sources>=2
            % Good optimized singles are useful starts, but retain the raw
            % alternatives: a single-source compromise need not equal a source.
            for s=1:numel(singles)
                seeds=[singles(s).position;seeds]; %#ok<AGROW>
            end
            seeds=unique_seeds(seeds,O.seed_merge_distance_m,...
                O.max_seed_candidates);
            pairs=pair_starts(seeds,x,y,b,H,scale,O);
            doubles=repmat(blank_action(),0,1);
            for s=1:min(O.starts_per_order,numel(pairs))
                [a,ncall]=optimize(pairs{s},x,y,b,H,scale,lo,hi,O);
                Info.fit_attempts=Info.fit_attempts+1;
                Info.objective_evaluations=Info.objective_evaluations+ncall;
                doubles=append_distinct(doubles,a,O);
            end
            doubles=sort_actions(doubles,O.max_hypotheses_per_order);
            for s=1:numel(doubles)
                Windows(w).candidates(end+1,1)=doubles(s); %#ok<AGROW>
            end
        end
        Info.accepted_hypotheses=Info.accepted_hypotheses+...
            numel(Windows(w).candidates)-1;
        Info.accepted_nonconverged_hypotheses=...
            Info.accepted_nonconverged_hypotheses+...
            sum(~[Windows(w).candidates(2:end).optimizer_converged]);
    end
    if O.verbose
        fprintf('  Local inversion: %d/%d windows; %.1f s\n',...
            w,numel(Windows),toc(started));
    end
end
Info.elapsed_seconds=toc(started);
end

function O=defaults(O)
names={'half_width_pixels','stride_pixels','max_sources',...
    'max_hypotheses_per_order','starts_per_order','max_iterations',...
    'max_function_evaluations','position_bounds',...
    'minimum_source_separation_m','model_penalty','relative_field_floor',...
    'signal_fraction_floor','background_order','moment_ridge',...
    'max_condition','seed_merge_distance_m','max_seed_candidates',...
    'minimum_component_fraction','verbose','residual_reference'};
values={10,10,2,2,3,180,1000,[-24 -24 .5;24 24 16],2,.015,.02,...
    .03,1,1e-8,1e10,2,8,0,false,.15};
for k=1:numel(names)
    if ~isfield(O,names{k}),O.(names{k})=values{k};end
end
assert(O.stride_pixels>=1 && O.stride_pixels==floor(O.stride_pixels),...
    'stride_pixels must be a positive integer.');
assert(O.background_order==0 || O.background_order==1,...
    'background_order must be 0 (constant) or 1 (plane).');
assert(isfinite(O.residual_reference) && O.residual_reference>0,...
    'residual_reference must be positive.');
end

function a=blank_action()
a=struct('position',zeros(0,3),'moment',zeros(0,3),...
    'background',zeros(0,1),'prediction',zeros(0,1),...
    'relative_residual',Inf,'cost',Inf,'K',0,'condition',Inf,'valid',false,...
    'optimizer_exitflag',NaN,'optimizer_converged',false);
end

function H=background_matrix(x,y,center,bounds,order)
H=ones(numel(x),1);
if order==1
    width=max(bounds(2,:)-bounds(1,:),eps);
    H=[H,(x-center(1))/width(1),(y-center(2))/width(2)];
end
end

function seeds=make_seeds(Data,Pool,ids,x,y,b,H,scale,lo,hi,O)
raw=zeros(0,3);
if isfield(Pool,'node_count') && Pool.node_count>0
    local=ismember(Pool.linear_indices,ids);
    raw=vertcat(Pool.position{local});
    if isempty(raw),raw=zeros(0,3);end
end
raw=raw(all(isfinite(raw),2),:);
raw=raw(all(bsxfun(@gt,raw,lo),2)&all(bsxfun(@lt,raw,hi),2),:);
% Deterministically cap screening work. No random draws or truth are used.
if size(raw,1)>100
    raw=raw(unique(round(linspace(1,size(raw,1),100))),:);
end
res=b-H*(H\b);[~,order]=sort(abs(res),'descend');
peak_xy=zeros(0,2);
step=sqrt(abs((Data.X(1,2)-Data.X(1,1))*...
    (Data.Y(2,1)-Data.Y(1,1))));
for j=1:numel(order)
    p=[x(order(j)),y(order(j))];
    if isempty(peak_xy) || all(sqrt(sum(bsxfun(@minus,peak_xy,p).^2,2))>...
            max(2*step,O.minimum_source_separation_m))
        peak_xy(end+1,:)=p; %#ok<AGROW>
    end
    if size(peak_xy,1)>=4,break;end
end
depths=lo(3)+(hi(3)-lo(3))*[.2 .45 .7];
for j=1:size(peak_xy,1)
    for d=1:numel(depths)
        raw(end+1,:)=[peak_xy(j,:),depths(d)]; %#ok<AGROW>
    end
end
raw=[raw;mean([lo;hi],1)];
raw=min(max(raw,lo+1e-4*(hi-lo)),hi-1e-4*(hi-lo));
score=Inf(size(raw,1),1);
for j=1:size(raw,1)
    a=linear_fit(raw(j,:),x,y,b,H,scale,O);
    score(j)=a.cost;
end
[~,ord]=sort(score);seeds=unique_seeds(raw(ord,:),...
    O.seed_merge_distance_m,O.max_seed_candidates);
end

function seeds=unique_seeds(raw,distance,max_count)
seeds=zeros(0,3);
for j=1:size(raw,1)
    if isempty(seeds) || all(sqrt(sum(bsxfun(@minus,seeds,raw(j,:)).^2,2))>distance)
        seeds(end+1,:)=raw(j,:); %#ok<AGROW>
    end
    if size(seeds,1)>=max_count,break;end
end
end

function pairs=pair_starts(seeds,x,y,b,H,scale,O)
pairs={};scores=zeros(0,1);
for a=1:size(seeds,1)-1
    for z=a+1:size(seeds,1)
        p=seeds([a,z],:);
        if norm(p(1,:)-p(2,:))<O.minimum_source_separation_m,continue;end
        fit=linear_fit(p,x,y,b,H,scale,O);
        if fit.valid
            pairs{end+1,1}=p;scores(end+1,1)=fit.cost; %#ok<AGROW>
        end
    end
end
[~,ord]=sort(scores);pairs=pairs(ord);
end

function [a,ncalls]=optimize(p0,x,y,b,H,scale,lo,hi,O)
K=size(p0,1);ncalls=0;
q=bsxfun(@rdivide,bsxfun(@minus,p0,lo),hi-lo);
q=min(max(q,1e-4),1-1e-4);u0=log(q./(1-q));
opts=optimset('Display','off','MaxIter',O.max_iterations,...
    'MaxFunEvals',O.max_function_evaluations,'TolX',1e-5,'TolFun',1e-8);
[u,~,exitflag]=fminsearch(@objective,u0(:),opts);
p=decode(u);a=linear_fit(p,x,y,b,H,scale,O);
a.optimizer_exitflag=exitflag;a.optimizer_converged=exitflag>0;
% Preserve a good initial fit if an iteration-limited run ended worse.
initial=linear_fit(p0,x,y,b,H,scale,O);
if initial.valid && (~a.valid || initial.cost<a.cost),a=initial;end
% Optional visibility restriction is OFF by default: a weak real source must
% not be silently excluded. If requested, apply after continuous fitting.
if a.valid && O.minimum_component_fraction>0
    for j=1:a.K
        component=dipole_bz_kernel(x,y,a.position(j,:))*a.moment(j,:)';
        if sqrt(mean(component.^2))/scale<O.minimum_component_fraction
            a.valid=false;break;
        end
    end
end
    function p=decode(v)
        qv=1./(1+exp(-min(max(reshape(v,K,3),-35),35)));
        p=bsxfun(@plus,lo,bsxfun(@times,qv,hi-lo));
    end
    function value=objective(v)
        ncalls=ncalls+1;
        pp=decode(v);
        if K==2 && norm(pp(1,:)-pp(2,:))<O.minimum_source_separation_m
            value=1e6+1e4*(1-norm(pp(1,:)-pp(2,:))/...
                O.minimum_source_separation_m)^2;
            return;
        end
        fit=linear_fit(pp,x,y,b,H,scale,O);
        if fit.valid,value=fit.cost;else,value=1e6;end
    end
end

function a=linear_fit(p,x,y,b,H,scale,O)
a=blank_action();K=size(p,1);a.K=K;a.position=p;
if K==2 && norm(p(1,:)-p(2,:))<O.minimum_source_separation_m,return;end
A=zeros(numel(b),3*K+size(H,2));
for j=1:K,A(:,3*j-2:3*j)=dipole_bz_kernel(x,y,p(j,:));end
A(:,3*K+1:end)=H;
norms=sqrt(sum(A.^2,1));
if any(~isfinite(norms)) || any(norms<=eps),return;end
An=bsxfun(@rdivide,A,norms);Gram=An'*An;
rc=rcond(Gram);
if ~isfinite(rc) || rc<=0,return;end
a.condition=1/sqrt(rc);
if a.condition>O.max_condition,return;end
ridge=zeros(size(Gram));ridge(1:3*K,1:3*K)=O.moment_ridge*eye(3*K);
coefficient=(Gram+ridge)\(An'*b);
coefficient=coefficient./norms';
prediction=A*coefficient;
if any(~isfinite(coefficient)) || any(~isfinite(prediction)),return;end
a.moment=reshape(coefficient(1:3*K),3,K)';
a.background=coefficient(3*K+1:end);a.prediction=prediction;
a.relative_residual=sqrt(mean((prediction-b).^2))/scale;
a.cost=(a.relative_residual/O.residual_reference)^2+O.model_penalty*K;
a.valid=isfinite(a.cost);
end

function list=append_distinct(list,a,O)
if ~a.valid,return;end
for j=1:numel(list)
    old=list(j);
    if a.K~=old.K,continue;end
    if a.K==1
        d=norm(a.position-old.position);
    else
        d=min(max(sqrt(sum((a.position-old.position).^2,2))),...
            max(sqrt(sum((a.position-old.position([2 1],:)).^2,2))));
    end
    if d<O.seed_merge_distance_m
        if a.cost<old.cost,list(j)=a;end
        return;
    end
end
list(end+1,1)=a;
end

function list=sort_actions(list,max_count)
if isempty(list),return;end
[~,ord]=sort([list.cost]);list=list(ord(1:min(numel(ord),max_count)));
end
