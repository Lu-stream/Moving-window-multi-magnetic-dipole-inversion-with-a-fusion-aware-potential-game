function State=mw_physical_setup(Data,Windows,Fusion,Physical)
%MW_PHYSICAL_SETUP Fixed measured-data-only global scoring ingredients.
% The original candidates and the original mw_fuse rules are not changed.
% A shared affine background is used once over the entire measured map.
% tau=multiplier*3*log(N) is a predeclared BIC-inspired source-count penalty;
% it is not claimed to be a calibrated BIC for data-selected positions.
if nargin<4,Physical=struct();end
if ~isfield(Physical,'tau_multiplier'),Physical.tau_multiplier=1;end
if ~isfield(Physical,'svd_relative_tolerance'),Physical.svd_relative_tolerance=1e-12;end
assert(Physical.tau_multiplier>=0&&isfinite(Physical.tau_multiplier));
assert(Physical.svd_relative_tolerance>0&&Physical.svd_relative_tolerance<1);
assert(isequal(size(Data.X),size(Data.Y),size(Data.Bz)));
assert(all(isfinite(Data.Bz(:)))&&all(isfinite(Data.X(:)))&&all(isfinite(Data.Y(:))));
assert(min(size(Data.Bz))>=2);
y=Data.Bz(:);N=numel(y);z=Data.Bz;
mixed=z(1:end-1,1:end-1)-z(2:end,1:end-1)-...
 z(1:end-1,2:end)+z(2:end,2:end);
kappa_MAD=.6744897501960817; % Standard-normal 75th percentile; same value as before.
sigma_hat=median(abs(mixed(:)-median(mixed(:))))/(2*kappa_MAD);
if sigma_hat==0,sigma_hat=eps(max(1,max(abs(y))));end
assert(sigma_hat>0&&isfinite(sigma_hat));
xc=mean(Data.X(:));yc=mean(Data.Y(:));
xs=max(abs(Data.X(:)-xc));ys=max(abs(Data.Y(:)-yc));
if xs==0,xs=1;end
if ys==0,ys=1;end
H=[ones(N,1),(Data.X(:)-xc)/xs,(Data.Y(:)-yc)/ys];
% Precompute each possible source field.  Reports remain distinct even when
% positions happen to coincide.  Their identities map exactly to mw_fuse's
% selected-report indices after its finite-position/moment filtering.
report_ids=cell(numel(Windows),1);positions=zeros(0,3);next=0;
valid_actions=cell(numel(Windows),1);initial=ones(numel(Windows),1);
for i=1:numel(Windows)
 C=Windows(i).candidates;report_ids{i}=cell(numel(C),1);
 costs=reshape([C.cost],[],1);
 for a=1:numel(C)
  valid=isfinite(costs(a));
  if isfield(C,'valid'),valid=valid&&C(a).valid;end
  if ~valid,costs(a)=Inf;continue;end
  valid_actions{i}(end+1)=a; %#ok<AGROW>
  if C(a).K==0,continue;end
  good=all(isfinite(C(a).position),2)&all(isfinite(C(a).moment),2);
  p=C(a).position(good,:);k=size(p,1);
  assert(all(p(:,3)>0));
  report_ids{i}{a}=(next+1:next+k).';
  positions=[positions;p];next=next+k; %#ok<AGROW>
 end
 assert(~isempty(valid_actions{i}),'A window has no valid finite candidate.');
 [~,initial(i)]=min(costs);
end
kernel=zeros(N,3*next);
for j=1:next
 kernel(:,3*j-2:3*j)=dipole_bz_kernel(Data.X,Data.Y,positions(j,:));
end
State=struct('Data',Data,'Windows',Windows,'Fusion',Fusion,...
 'Physical',Physical,'y',y,'N',N,'background_design',H,...
 'background_center_xy',[xc yc],'background_scale_xy',[xs ys],...
 'sigma_hat_nT',sigma_hat,'tau',Physical.tau_multiplier*3*log(N),...
 'report_ids',{report_ids},'report_positions',positions,'kernel',kernel,...
 'valid_actions',{valid_actions},'initial_actions',initial);
end
