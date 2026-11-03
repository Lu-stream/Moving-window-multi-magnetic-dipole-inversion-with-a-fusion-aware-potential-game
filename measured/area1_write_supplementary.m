function area1_write_supplementary(Results,output_dir)
% Transparent supplementary tables from exactly the same estimated sets.
c=Results.Metadata.evaluation_tolerance_um;
fid=fopen(fullfile(output_dir,'supplementary_metrics.csv'),'w');
assert(fid>=0,'Cannot create supplementary metrics file.');
fprintf(fid,'method,Recall,Precision,F1,F2,GOSPA_um,asymmetric_GOSPA_r2_um,asymmetric_GOSPA_r4_um\n');
for j=1:numel(Results.Methods)
 M=Results.Methods{j}.Metrics;
 F2=5*M.TP/(5*M.TP+M.FP+4*M.FN);
 delta=M.output_count-size(Results.Reference.position,1);
 w2=M.GOSPA+c*(1/3-1/2)*delta;
 w4=M.GOSPA+c*(1/5-1/2)*delta;
 fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',...
  Results.Methods{j}.name,M.Recall,M.Precision,M.F1,F2,M.GOSPA,w2,w4);
end
fclose(fid);
fid=fopen(fullfile(output_dir,'tolerance_comparison.csv'),'w');
assert(fid>=0,'Cannot create tolerance table.');
fprintf(fid,'tolerance_um,method,outputs,CT_matches,unmatched_outputs,unmatched_CT,Recall,Precision,F1,GOSPA_um\n');
for t=1:numel(Results.Metadata.tolerance_grid_um)
 for j=1:numel(Results.Methods)
  M=Results.ToleranceResults{j,t}.Metrics;
  fprintf(fid,'%.12g,%s,%d,%d,%d,%d,%.12g,%.12g,%.12g,%.12g\n',...
   Results.Metadata.tolerance_grid_um(t),Results.Methods{j}.name,M.output_count,...
   M.TP,M.FP,M.FN,M.Recall,M.Precision,M.F1,M.GOSPA);
 end
end
fclose(fid);
fid=fopen(fullfile(output_dir,'source_estimates.csv'),'w');
assert(fid>=0,'Cannot create source estimate table.');
fprintf(fid,'method,estimate_id,registered_x_um,registered_y_um,surface_depth_um,matched_CT_ID,matched_distance_um,nearest_CT_distance_um,mx_Am2,my_Am2,mz_Am2\n');
for j=1:numel(Results.Methods)
 M=Results.Methods{j};n=size(M.Fit.position,1);
 matched_id=nan(n,1);matched_distance=nan(n,1);
 e=M.Details.matched_estimate_indices(:);t=M.Details.matched_truth_indices(:);
 matched_id(e)=Results.Reference.ids(t);
 matched_distance(e)=M.Details.matched_errors(:);
 reg=bsxfun(@plus,M.Fit.position,[Results.Metadata.xy_offset_um -Results.Metadata.sensor_height_um]);
 nearest=min(M.Details.distance_matrix,[],1);
 for k=1:n
  fprintf(fid,'%s,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',...
   M.name,k,reg(k,:),matched_id(k),matched_distance(k),nearest(k),M.Fit.moment_Am2(k,:));
 end
end
fclose(fid);
end
