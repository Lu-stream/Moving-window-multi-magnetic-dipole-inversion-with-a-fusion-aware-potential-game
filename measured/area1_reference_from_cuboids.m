function Reference=area1_reference_from_cuboids(input_dir,output_dir,xy_offset_um,sensor_height_um)
% Evaluation-only CT reference. The half-sizes a,b,c give volume 8*a*b*c.
cubes=dlmread(fullfile(input_dir,'FWInput-FineCuboids-A1.txt'));
assert(size(cubes,2)==7&&all(all(cubes(:,4:6)>0)),'Invalid CT cuboid table.');
ids=unique(cubes(:,7));centers=zeros(numel(ids),3);volume=zeros(numel(ids),1);
counts=zeros(numel(ids),1);
for j=1:numel(ids)
 q=cubes(cubes(:,7)==ids(j),:);v=8*prod(q(:,4:6),2);
 volume(j)=sum(v);counts(j)=size(q,1);
 centers(j,:)=sum(bsxfun(@times,q(:,1:3),v),1)/volume(j);
end
Reference=struct('position',bsxfun(@plus,centers,[-xy_offset_um sensor_height_um]),...
 'ids',ids,'volume_um3',volume,...
 'meaning','CT volume-weighted centers; geometric references, not magnetic ground truth');
fid=fopen(fullfile(output_dir,'CT_reference_centers.csv'),'w');
assert(fid>=0,'Cannot create CT reference table.');
fprintf(fid,'ID,cuboids,volume_um3,registered_x_um,registered_y_um,surface_depth_um\n');
for j=1:numel(ids)
 fprintf(fid,'%d,%d,%.15g,%.15g,%.15g,%.15g\n',ids(j),counts(j),volume(j),centers(j,:));
end
fclose(fid);
end
