function [M,D,E]=multi_evaluate_output(Fit,Scene,PI,is_refit,c)
% Evaluation only: source truth is never passed to inference.
[M,D]=pilot_evaluate(Fit,Scene.Truth,c);
source=zeros(numel(Scene.X),1);
for k=1:size(Fit.position,1)
 source=source+dipole_bz_kernel(Scene.X,Scene.Y,Fit.position(k,:))*Fit.moment(k,:)';
end
if is_refit,prediction=PI.prediction(:);else,prediction=source;end
clean=Scene.Bclean(:);obs=Scene.Bobserved(:);
merr=zeros(numel(D.matched_truth_indices),1);
magerr=merr;pospct=merr;depthpct=merr;
for k=1:numel(merr)
 ti=D.matched_truth_indices(k);ei=D.matched_estimate_indices(k);
 den=norm(Scene.Truth.moment(ti,:));
 merr(k)=norm(Fit.moment(ei,:)-Scene.Truth.moment(ti,:))/den;
 magerr(k)=abs(norm(Fit.moment(ei,:))-den)/den;
 pospct(k)=100*D.matched_errors(k)/norm(Scene.Truth.position(ti,:));
 depthpct(k)=100*D.matched_errors(k)/Scene.Truth.position(ti,3);
end
if isempty(merr),mm=NaN;mg=NaN;mp=NaN;md=NaN;
else,mm=mean(merr);mg=mean(magerr);mp=mean(pospct);md=mean(depthpct);end
E=struct('matched_moment_relative_errors',merr,'mean_matched_moment_relative_error',mm,...
 'matched_magnitude_relative_errors',magerr,'mean_matched_magnitude_relative_error',mg,...
 'matched_position_errors_pct',pospct,'mean_matched_position_error_pct',mp,...
 'matched_depth_normalized_position_errors_pct',depthpct,...
 'mean_matched_depth_normalized_position_error_pct',md,...
 'source_field_relative_error',norm(source-clean)/norm(clean),...
 'total_field_relative_error',norm(prediction-clean)/norm(clean),...
 'observed_residual_RMSE_nT',sqrt(mean((prediction-obs).^2)),...
 'observed_relative_residual',norm(prediction-obs)/norm(obs));
end
