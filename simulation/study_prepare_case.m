function [Metrics,Fields]=study_prepare_case(Result)
%STUDY_PREPARE_CASE Reconstruct the three saved predictions without refitting.
% Uses saved moments and background coefficients, and checks saved residuals.
assert(isstruct(Result)&&isfield(Result,'Meta')&&isfield(Result,'Scene')&&isfield(Result,'Methods'));
names={'primin_direct','independent_refit','physical_game'};
labels={'Primin-based','Independent + fusion','Fusion-aware game'};
S=Result.Scene;assert(isequal(size(S.X),size(S.Y),size(S.Bobserved),size(S.Bclean)));
Fields=struct('X',S.X,'Y',S.Y,'measured_nT',S.Bobserved,'clean_nT',S.Bclean,...
    'method_names',{names},'source_nT',{{}},'background_nT',{{}},'prediction_nT',{{}},'residual_nT',{{}});
proto=struct('method_name','','method_label','','TP',NaN,'FP',NaN,'FN',NaN,...
    'Recall',NaN,'Precision',NaN,'F1',NaN,'GOSPA_m',NaN,'matched_position_error_m',NaN,...
    'matched_magnitude_error_percent',NaN,'residual_RMSE_nT',NaN,'residual_relative',NaN,...
    'saved_residual_check_relative_error',NaN,'position_tolerance_m',NaN);
Metrics=repmat(proto,1,numel(names));
for j=1:numel(names)
    M=get_method(Result,names{j});assert(~isempty(M),'Missing method %s in %s.',names{j},Result.Meta.name);
    q=proto;q.method_name=names{j};q.method_label=labels{j};
    copied={'TP','FP','FN','Recall','Precision','F1'};
    for n=1:numel(copied),q.(copied{n})=M.Metrics.(copied{n});end
    q.GOSPA_m=M.Metrics.GOSPA;
    q.matched_position_error_m=M.Metrics.mean_matched_position_error;
    q.matched_magnitude_error_percent=100*M.Evaluation.mean_matched_magnitude_relative_error;
    q.position_tolerance_m=M.Details.cutoff;
    source=zeros(numel(S.X),1);
    for k=1:size(M.Fit.position,1)
        source=source+dipole_bz_kernel(S.X,S.Y,M.Fit.position(k,:))*M.Fit.moment(k,:).';
    end
    background=zeros(size(source));
    if isfield(M,'PhysicalInfo')&&isfield(M.PhysicalInfo,'fitted_background_coefficients')
        PI=M.PhysicalInfo;beta=PI.fitted_background_coefficients(:);
        assert(numel(beta)==3,'Expected the saved affine background coefficients.');
        center=PI.background_center_xy(:).';scale=PI.background_scale_xy(:).';
        assert(numel(center)==2&&numel(scale)==2&&all(scale>0),'Invalid saved background normalization.');
        H=[ones(numel(S.X),1),(S.X(:)-center(1))/scale(1),(S.Y(:)-center(2))/scale(2)];
        background=H*beta;
    elseif ~strcmp(names{j},'primin_direct')
        error('Saved affine background is missing for %s; cannot reconstruct the original prediction.',names{j});
    end
    prediction=source+background;residual=S.Bobserved(:)-prediction;
    q.residual_RMSE_nT=sqrt(mean(residual.^2));
    q.residual_relative=norm(residual)/max(norm(S.Bobserved(:)),eps);
    expected=M.Evaluation.observed_residual_RMSE_nT;
    q.saved_residual_check_relative_error=abs(q.residual_RMSE_nT-expected)/max(1,abs(expected));
    assert(q.saved_residual_check_relative_error<1e-9,...
        'Reconstructed residual disagrees with saved Evaluation for %s (relative discrepancy %.3g).',...
        names{j},q.saved_residual_check_relative_error);
    expected_relative=M.Evaluation.observed_relative_residual;
    assert(abs(q.residual_relative-expected_relative)<1e-9*max(1,abs(expected_relative)),...
        'Relative residual does not reproduce the saved evaluation for %s.',names{j});
    Fields.source_nT{j}=reshape(source,size(S.X));
    Fields.background_nT{j}=reshape(background,size(S.X));
    Fields.prediction_nT{j}=reshape(prediction,size(S.X));
    Fields.residual_nT{j}=reshape(residual,size(S.X));
    Metrics(j)=q;
end
end

function M=get_method(R,name)
M=[];A=R.Methods;if ~iscell(A),A=num2cell(A);end
for k=1:numel(A),if strcmp(A{k}.name,name),M=A{k};return;end,end
end
