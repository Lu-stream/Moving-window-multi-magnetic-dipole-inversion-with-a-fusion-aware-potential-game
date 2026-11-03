function pool=pg_pool(B,X,Y,opt)
% Same candidate equations and bounds as earlier a3_build_pool; NO scores here.
% Separating this stage prevents double-counting local scoring in timing.
dx=X(1,2)-X(1,1);dy=Y(2,1)-Y(1,1);
T=primin_fft_tensor(B,dx,dy,opt.pad_factor);
distance=max(abs(X-opt.center_m(1)),abs(Y-opt.center_m(2)));
outer=distance<=max(opt.window_widths_m)/2;
eligible=outer & abs(B)>opt.min_b_fraction*max(abs(B(outer)));
stride=false(size(B));stride(1:opt.candidate_stride:end,1:opt.candidate_stride:end)=true;
eligible=eligible & stride;r=opt.patch_radius_pixels;
eligible(1:r,:)=false;eligible(end-r+1:end,:)=false;
eligible(:,1:r)=false;eligible(:,end-r+1:end)=false;
ids=find(eligible);positions={};moments={};indices=[];
for k=1:numel(ids)
 index=ids(k);G=[T.xx(index),T.xy(index),T.xz(index);...
 T.xy(index),T.yy(index),T.yz(index);T.xz(index),T.yz(index),T.zz(index)];
 c=primin_point_candidates(B(index),G,[X(index);Y(index);0]);
 good=all(bsxfun(@gt,c.position,opt.position_bounds(1,:)),2) & ...
  all(bsxfun(@lt,c.position,opt.position_bounds(2,:)),2);
 p=c.position(good,:);m=c.moment(good,:);if isempty(p),continue;end
 [p,ord]=sortrows(p,[1 2 3]);positions{end+1,1}=p;moments{end+1,1}=m(ord,:); %#ok<AGROW>
 indices(end+1,1)=index; %#ok<AGROW>
end
pool=struct('position',{positions},'moment',{moments},'linear_indices',indices,...
 'node_count',numel(indices),'branch_count',sum(cellfun(@(p)size(p,1),positions)),...
 'field_floor_nT',opt.relative_field_floor*sqrt(mean(B(outer).^2)));
end
