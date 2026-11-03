function S=pg_scores(pool,B,X,Y,P,domain)
% Same batched forward kernel for N and G. NO parameter refit or truth access.
% N: 9x9 minus generating pixel. G: whole 61x61 observed survey.
n=pool.node_count;S.error=cell(n,1);S.sample_count=0;S.domain=domain;
if n==0,S.pixels_per_branch=0;return;end
sz=cellfun(@(p)size(p,1),pool.position);owner=repelem((1:n).',sz);
p=vertcat(pool.position{:});m=vertcat(pool.moment{:});score=zeros(size(p,1),1);
r=P.patch_radius_pixels;[dr,dc]=ndgrid(-r:r,-r:r);take=dr~=0|dc~=0;
offset=dr(take)+size(B,1)*dc(take);centers=pool.linear_indices(:);
batch=P.batch_size;
if isfield(P,'max_score_elements')
 if strcmp(domain,'N'),np=(2*r+1)^2-1;else,np=numel(B);end
 batch=max(1,min(batch,floor(P.max_score_elements/max(np,1))));
end
for first=1:batch:size(p,1)
 j=first:min(first+batch-1,size(p,1));
 if strcmp(domain,'N')
  ix=bsxfun(@plus,offset(:),centers(owner(j)).');
  x=X(ix);y=Y(ix);b=B(ix);
 else
  assert(strcmp(domain,'G'));x=X(:);y=Y(:);b=B(:);
 end
 dx=bsxfun(@minus,x,p(j,1).');dy=bsxfun(@minus,y,p(j,2).');d=p(j,3).';
 r2=bsxfun(@plus,dx.^2+dy.^2,d.^2);
 md=bsxfun(@times,dx,m(j,1).')+bsxfun(@times,dy,m(j,2).');
 md=bsxfun(@plus,md,d.*m(j,3).');
 pred=100*(3*bsxfun(@times,md,d)./r2.^2.5-bsxfun(@rdivide,m(j,3).',r2.^1.5));
 den=max(sqrt(sum(b.^2,1)),pool.field_floor_nT*sqrt(size(b,1)));den=max(den,realmin);
 score(j)=sqrt(sum(bsxfun(@minus,pred,b).^2,1)).'./den(:);
 S.sample_count=S.sample_count+numel(pred);
end
stop=cumsum(sz);start=1;
for i=1:n,S.error{i}=score(start:stop(i));start=stop(i)+1;end
S.pixels_per_branch=S.sample_count/max(pool.branch_count,1);
end
