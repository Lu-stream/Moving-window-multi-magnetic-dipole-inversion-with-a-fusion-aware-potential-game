function S=mw_scene(Grid,Source,noise_fraction,seed)
% Forward simulation only. Inversion functions never receive this truth.
assert(Grid.n>=9 && Grid.n==round(Grid.n) && Grid.scan_half_width_m>0);
assert(noise_fraction>=0 && isfinite(noise_fraction));
p=Source.position;u=Source.direction;a=Source.magnitude(:);
assert(size(p,2)==3 && isequal(size(p),size(u)) && numel(a)==size(p,1));
assert(all(p(:,3)>0) && all(a>0));
un=sqrt(sum(u.^2,2));assert(all(un>0));
m=bsxfun(@times,bsxfun(@rdivide,u,un),a);
v=linspace(-Grid.scan_half_width_m,Grid.scan_half_width_m,Grid.n);[X,Y]=meshgrid(v,v);
B=zeros(size(X));parts=zeros([size(X),size(p,1)]);
for k=1:size(p,1)
 parts(:,:,k)=reshape(dipole_bz_kernel(X,Y,p(k,:))*m(k,:).',size(X));B=B+parts(:,:,k);
end
sigma=noise_fraction*sqrt(mean(B(:).^2));
old=rng;cleanup=onCleanup(@()rng(old));rng(seed,'twister'); %#ok<NASGU>
Noise=sigma*randn(size(B));
Truth=struct('position',p,'moment',m);
S=struct('X',X,'Y',Y,'Bclean',B,'Bobserved',B+Noise,'Noise',Noise,...
 'Truth',Truth,'SourceParts',parts,'Sigma_nT',sigma,'NoiseFraction',noise_fraction,...
 'NoiseSeed',seed,'GridStep_m',v(2)-v(1));
end
