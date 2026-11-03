function R=wg_compact_result(R)
% Keep final numerical outputs while dropping reproducible working arrays.
% This does not modify any estimate, evaluation, action or scalar diagnostic.
if isfield(R,'Base'),R=rmfield(R,'Base');end
if isfield(R.Scene,'SourceParts'),R.Scene=rmfield(R.Scene,'SourceParts');end
large_fields={'fitted_background','prediction','source_field','residual',...
 'predicted_field_nT','source_field_nT','background_field_nT'};
for j=1:numel(R.Methods)
 if iscell(R.Methods),M=R.Methods{j};else,M=R.Methods(j);end
 if isfield(M,'PhysicalInfo')&&isstruct(M.PhysicalInfo)
  names=large_fields(isfield(M.PhysicalInfo,large_fields));
  if ~isempty(names),M.PhysicalInfo=rmfield(M.PhysicalInfo,names);end
 end
 if iscell(R.Methods),R.Methods{j}=M;else,R.Methods(j)=M;end
end
R.PackageCacheVersion='demo-stability-compact-v1';
end
