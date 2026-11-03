function [Data,Summary,Rows]=study_summarize(Results,layout_ids,noise_fractions,output_dir)
%STUDY_SUMMARIZE Paired layout-by-noise results for the three named methods.
% No reinversion or selection is performed. Each noise-level summary uses the
% same layouts separately; no pooling across noise levels is performed.
% Data metric arrays have dimensions layout x noise x method. Fractions stay
% in [0,1] in the data and CSV files; the plotting helper displays percentages.
if nargin<4,output_dir='';end
layout_ids=layout_ids(:).';noise_fractions=noise_fractions(:).';
assert(iscell(Results),'Results must be a cell array.');
assert(~isempty(layout_ids)&&numel(unique(layout_ids))==numel(layout_ids),'Layout IDs must be unique.');
assert(~isempty(noise_fractions)&&all(isfinite(noise_fractions))&&all(noise_fractions>=0),...
    'Noise fractions must be finite and nonnegative.');
assert(numel(unique(noise_fractions))==numel(noise_fractions),'Noise levels must be unique.');
nl=numel(layout_ids);nn=numel(noise_fractions);
assert(numel(Results)==nl*nn,'Expected one completed Result for every layout/noise combination.');
names={'primin_direct','independent_refit','physical_game'};
labels={'Primin-based','Independent + fusion','Fusion-aware game'};
metrics={'Recall','Precision','F1','TP','FP','FN','GOSPA_m','output_count',...
    'matched_position_error_m','residual_RMSE_nT','residual_relative'};
Data=struct('layout_ids',layout_ids,'noise_fractions',noise_fractions,...
    'method_names',{names},'method_labels',{labels},'metric_fields',{metrics},...
    'case_names',{cell(nl,nn)},'original_layout_ids',NaN(1,nl),...
    'original_case_names',{cell(nl,nn)},'source_count',NaN,...
    'selection_note','These fixed layouts were selected using results at 5% noise.');
for k=1:numel(metrics),Data.(metrics{k})=NaN(nl,nn,3);end
seen=false(nl,nn);reference_truth=cell(nl,1);
proto=struct('layout_index',0,'realization',0,'original_realization',0, ...
    'noise_fraction',0,'noise_percent',0,...
    'case_name','','original_case_name','','method','','method_label','','source_count',0);
for k=1:numel(metrics),proto.(metrics{k})=NaN;end
Rows=repmat(proto,nl*nn*3,1);nr=0;
for index=1:numel(Results)
    R=Results{index};assert(isstruct(R)&&isfield(R,'Meta')&&isfield(R,'Methods'),...
        'Every cell must contain a completed Result structure.');
    il=find(layout_ids==R.Meta.realization,1);
    in=find(abs(noise_fractions-R.Meta.noise_fraction)<=1e-12,1);
    assert(~isempty(il)&&~isempty(in),'Result metadata does not match the requested layout/noise grid.');
    assert(~seen(il,in),'Duplicate layout/noise result.');seen(il,in)=true;
    nt=size(R.Scene.Truth.position,1);
    if isnan(Data.source_count),Data.source_count=nt;end
    assert(nt==Data.source_count,'Every layout must contain the same number of true sources.');
    truth=R.Scene.Truth;
    if isempty(reference_truth{il}),reference_truth{il}=truth;
    else
        previous=reference_truth{il};a=truth.moment;b=previous.moment;
        same_moment=isequal(size(a),size(b))&&all(isfinite(a(:)))&&all(isfinite(b(:)))&&...
            max(abs(a(:)-b(:)))<=1e-12*max(1,max(abs(b(:))));
        assert(isequal(truth.position,previous.position)&&same_moment,...
            'The source layout or magnetic moments changed between noise levels.');
    end
    Data.case_names{il,in}=R.Meta.name;
    if ~isnan(Data.original_layout_ids(il))
        assert(Data.original_layout_ids(il)==R.Meta.original_realization,...
            'Original layout ID changed between noise levels.');
    end
    Data.original_layout_ids(il)=R.Meta.original_realization;
    Data.original_case_names{il,in}=R.Meta.original_name;
    for jm=1:3
        M=lookup_method(R,names{jm});q=proto;
        q.layout_index=il;q.realization=layout_ids(il);
        q.original_realization=R.Meta.original_realization;
        q.noise_fraction=noise_fractions(in);
        q.noise_percent=100*noise_fractions(in);q.case_name=R.Meta.name;
        q.original_case_name=R.Meta.original_name;
        q.method=names{jm};q.method_label=labels{jm};q.source_count=nt;
        copy={'Recall','Precision','F1','TP','FP','FN','output_count'};
        for k=1:numel(copy),q.(copy{k})=M.Metrics.(copy{k});end
        q.GOSPA_m=M.Metrics.GOSPA;
        q.matched_position_error_m=M.Metrics.mean_matched_position_error;
        q.residual_RMSE_nT=M.Evaluation.observed_residual_RMSE_nT;
        q.residual_relative=M.Evaluation.observed_relative_residual;
        assert(q.TP+q.FN==nt&&q.TP+q.FP==q.output_count,'Inconsistent detection counts.');
        assert(abs(q.Recall-q.TP/nt)<1e-10,'Recall does not agree with TP/true-source count.');
        assert(abs(q.F1-2*q.TP/(nt+q.output_count))<1e-10,'F1 does not agree with counts.');
        for k=1:numel(metrics),Data.(metrics{k})(il,in,jm)=q.(metrics{k});end
        nr=nr+1;Rows(nr)=q;
    end
end
assert(all(seen(:)),'A layout/noise result is missing.');
% Reorder the long table to the user-facing selected-layout order, then noise,
% then method. Original IDs are retained in original_realization.
[~,order]=sortrows([[Rows.layout_index].',[Rows.noise_fraction].'],[1 2]);
Rows=Rows(order);
sp=struct('noise_fraction',0,'noise_percent',0,'method','','method_label','',...
    'n_layouts',nl,'source_count',Data.source_count);
for k=1:numel(metrics)
    field=metrics{k};sp.([field '_mean'])=NaN;sp.([field '_median'])=NaN;
    sp.([field '_SD'])=NaN;sp.([field '_n_finite'])=0;
end
Summary=repmat(sp,nn*3,1);ns=0;
for in=1:nn
    for jm=1:3
        q=sp;q.noise_fraction=noise_fractions(in);q.noise_percent=100*noise_fractions(in);
        q.method=names{jm};q.method_label=labels{jm};
        for k=1:numel(metrics)
            field=metrics{k};values=Data.(field)(:,in,jm);values=values(isfinite(values));
            q.([field '_n_finite'])=numel(values);
            if ~isempty(values)
                q.([field '_mean'])=mean(values);q.([field '_median'])=median(values);
                q.([field '_SD'])=std(values,0); % Sample SD, denominator n-1.
            end
        end
        assert(q.Recall_n_finite==nl&&q.F1_n_finite==nl&&q.FP_n_finite==nl&&q.GOSPA_m_n_finite==nl,...
            'Primary comparison metrics must be available for every selected layout.');
        ns=ns+1;Summary(ns)=q;
    end
end
if ~isempty(output_dir)
    if ~exist(output_dir,'dir'),mkdir(output_dir);end
    write_struct_csv(fullfile(output_dir,'all_case_metrics.csv'),Rows);
    write_struct_csv(fullfile(output_dir,'noise_summary.csv'),Summary);
    save(fullfile(output_dir,'ensemble_metrics.mat'),'Data','Summary','Rows','-v7');
end
end

function M=lookup_method(R,name)
A=R.Methods;if ~iscell(A),A=num2cell(A);end
which=[];for k=1:numel(A),if strcmp(A{k}.name,name),which(end+1)=k;end,end %#ok<AGROW>
assert(numel(which)==1,'Expected exactly one method named %s.',name);M=A{which};
end

function write_struct_csv(path,rows)
fields=fieldnames(rows);fid=fopen(path,'w');assert(fid>=0,'Cannot write %s.',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
for k=1:numel(fields)
    if k>1,fprintf(fid,',');end;fprintf(fid,'%s',fields{k});
end
fprintf(fid,'\n');
for i=1:numel(rows)
    for k=1:numel(fields)
        if k>1,fprintf(fid,',');end;v=rows(i).(fields{k});
        if ischar(v),fprintf(fid,'"%s"',strrep(v,'"','""'));
        elseif isnumeric(v)&&isscalar(v),fprintf(fid,'%.17g',v);
        else,error('CSV fields must be character vectors or numeric scalars.');end
    end
    fprintf(fid,'\n');
end
end
