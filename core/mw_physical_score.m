function [Fit,Info,FusionInfo]=mw_physical_score(actions,State)
%MW_PHYSICAL_SCORE Evaluate the actual fused output against every datum once.
% Positions and inclusion are the unchanged mw_fuse output.  Magnetic moments
% and ONE global affine background are fitted jointly and deterministically.
% Column-normalized SVD supplies the minimum-norm solution in scaled variables
% with a fixed relative rank threshold, without an additional ridge penalty.
started=tic;actions=actions(:);
assert(numel(actions)==numel(State.Windows));
for i=1:numel(actions)
 assert(any(State.valid_actions{i}==actions(i)),'Invalid source-set action.');
end
[Fit,FusionInfo]=mw_fuse(State.Windows,actions,State.Data,State.Fusion);
selected_ids=zeros(0,1);
for i=1:numel(actions)
 selected_ids=[selected_ids;State.report_ids{i}{actions(i)}(:)]; %#ok<AGROW>
end
K=size(Fit.position,1);
if K>0
 representatives=reshape([Fit.clusters.representative_report_id],[],1);
 global_ids=selected_ids(representatives);
 assert(isequal(Fit.position,State.report_positions(global_ids,:)),...
  'Global kernel cache does not correspond to fused source representatives.');
 cols=reshape(bsxfun(@plus,3*(global_ids(:).'-1),(1:3).'),1,[]);
 G=State.kernel(:,cols);
else
 global_ids=zeros(0,1);G=zeros(State.N,0);
end
A=[G,State.background_design];
scales=sqrt(sum(A.^2,1));
assert(all(scales>0)&&all(isfinite(scales)));
An=bsxfun(@rdivide,A,scales);
[U,D,V]=svd(An,'econ');s=diag(D);
threshold=State.Physical.svd_relative_tolerance*s(1);
keep=s>threshold;nrank=sum(keep);
scaled_coefficients=V(:,keep)*((U(:,keep)'*State.y)./s(keep));
coefficients=scaled_coefficients./scales.';
prediction=A*coefficients;residual=State.y-prediction;
SSE=sum(residual.^2);data_loss=.5*SSE/State.sigma_hat_nT^2;
if K>0,moments=reshape(coefficients(1:3*K),3,K).';else,moments=zeros(0,3);end
original_moments=Fit.moment;Fit.moment=moments;
for k=1:K,Fit.clusters(k).moment=moments(k,:);end
beta=coefficients(3*K+1:end);source_field=G*coefficients(1:3*K);
background=State.background_design*beta;
if nrank<size(A,2),condition=Inf;else,condition=s(1)/s(end);end
Info=struct('J',data_loss+State.tau*K,'data_loss',data_loss,...
 'source_penalty',State.tau*K,'tau',State.tau,'K',K,...
 'sigma_hat_nT',State.sigma_hat_nT,'SSE',SSE,'RMSE',sqrt(SSE/State.N),...
 'relative_residual',sqrt(SSE)/max(norm(State.y),eps),...
 'rank',nrank,'number_of_columns',size(A,2),'condition',condition,...
 'singular_values',s,'svd_threshold',threshold,...
 'fitted_background_coefficients',beta,...
 'fitted_background',reshape(background,size(State.Data.Bz)),...
 'prediction',reshape(prediction,size(State.Data.Bz)),...
 'source_field',reshape(source_field,size(State.Data.Bz)),...
 'residual',reshape(residual,size(State.Data.Bz)),...
 'original_moments',original_moments,'global_report_ids',global_ids,...
 'background_center_xy',State.background_center_xy,...
 'background_scale_xy',State.background_scale_xy,...
 'score_seconds',toc(started));
% Descriptive aliases used by the comparison/analysis driver.
Info.objective=Info.J;Info.condition_number=Info.condition;
Info.predicted_field_nT=Info.prediction;
Info.source_field_nT=Info.source_field;
Info.background_field_nT=Info.fitted_background;
FusionInfo.refinement='global fixed-position moments and affine background';
end
