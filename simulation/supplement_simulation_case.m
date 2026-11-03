function [Result,RunInfo] = supplement_simulation_case(root_dir,id,noise_fraction,S,mode)
% Explicit recomputation or archived-result reading. No cache fallback.
assert(isscalar(noise_fraction)&&isfinite(noise_fraction)&&noise_fraction>=0);
assert(any(strcmp(mode,{'recompute','archived'})));
p=fullfile(root_dir,'data','synthetic',sprintf('realization_%02d.mat',id));
assert(exist(p,'file')==2,'Missing frozen input: %s',p);
f=load(p);
id_map=dlmread(fullfile(root_dir,'data','synthetic','layout_id_map.csv'),',',1,0);
assert(isscalar(id)&&id==floor(id)&&id>=1&&id<=size(id_map,1));
assert(id_map(id,1)==id&&f.Metadata.realization==id_map(id,2), ...
    'Frozen input does not match its public/original ID mapping.');
original_id=f.Metadata.realization;
assert(isequaln(f.Grid,S.Grid)&&isequaln(f.SourceOptions,S.SourceOptions), ...
    'Restore the fixed grid/source settings before using these inputs.');
assert(size(f.Source.position,1)==S.SourceCount);
Scene=mw_scene(S.Grid,f.Source,noise_fraction,f.Metadata.noise_seed);
assert(isequal(size(f.StandardNoise),size(Scene.Bclean)));
Scene.Noise=Scene.Sigma_nT*f.StandardNoise;
Scene.Bobserved=Scene.Bclean+Scene.Noise;
name=sprintf('realization_%02d_noise_%03d',id,round(1000*noise_fraction));
assert(abs(noise_fraction-round(1000*noise_fraction)/1000)<1e-12, ...
    'Use noise fractions in steps of 0.001 for unique case filenames.');
original_name=sprintf('realization_%02d_noise_%03d',original_id,round(1000*noise_fraction));
Meta=struct('name',name,'realization',id, ...
    'original_realization',original_id,'original_name',original_name, ...
    'split','selected_illustration', ...
    'geometry_seed',f.Metadata.geometry_seed,'noise_seed',f.Metadata.noise_seed, ...
    'noise_fraction',noise_fraction);
RunInfo=struct('mode',mode,'case_name',name,'realization',id, ...
    'original_realization',original_id,'original_case_name',original_name, ...
    'computed',false, ...
    'loaded_archived_result',false,'software',version,'archive_comparison',struct());
archive=fullfile(root_dir,'archive','synthetic',[name '.mat']);
if strcmp(mode,'archived')
    assert(exist(archive,'file')==2,'No archived case: %s. No recompute fallback is used.',archive);
    old=load(archive,'Result','Settings');
    assert(old.Result.Meta.realization==original_id, 'Archive original ID differs.');
    assert(wg_same_scene(old.Result.Scene,Scene),'Archived and frozen scene differ.');
    keys={'Grid','SourceOptions','Candidate','Primin','Local','Game','Fusion','Physical','Evaluation'};
    for k=1:numel(keys)
        assert(isequaln(old.Settings.(keys{k}),S.(keys{k})), ...
            'Archive settings differ in %s. Archived mode cannot use changed settings.',keys{k});
    end
    Result=select_three_methods(old.Result);
    RunInfo.loaded_archived_result=true;
    fprintf('%s: reading archived estimates.\n',name);
else
    fprintf('%s: computing all three methods from scratch.\n',name);
    Result=supplement_run_three_methods(Scene,f.SourceInfo,Meta,S);
    RunInfo.computed=true;
    % Differences are reported, never replaced with favorable saved values.
    if exist(archive,'file')==2
        old=load(archive,'Result');
        RunInfo.archive_comparison=compare_results(Result,select_three_methods(old.Result));
    end
end
% Public run labels are consecutive; original MAT payloads remain unchanged.
Result.Meta.original_realization=original_id;
Result.Meta.original_name=original_name;
Result.Meta.realization=id;Result.Meta.name=name;
if isfield(Result,'Summary')
    Result.Summary.original_realization=original_id;
    Result.Summary.original_name=original_name;
    Result.Summary.realization=id;Result.Summary.name=name;
end
end

function R=select_three_methods(R)
names={'primin_direct','independent_refit','physical_game'};
A=R.Methods;if ~iscell(A),A=num2cell(A);end
B=cell(3,1);
for k=1:3
    ix=find(cellfun(@(m)strcmp(m.name,names{k}),A));
    assert(numel(ix)==1,'Expected one method named %s.',names{k});B{k}=A{ix};
end
R.Methods=B;
if isfield(R,'Summary')&&isfield(R.Summary,'Methods')
    A=R.Summary.Methods;if ~iscell(A),A=num2cell(A);end
    B=cell(3,1);
    for k=1:3
        ix=find(cellfun(@(m)strcmp(m.name,names{k}),A));
        assert(numel(ix)==1);B{k}=A{ix};
    end
    R.Summary.Methods=B;
end
R=wg_compact_result(R);
end

function Comparison=compare_results(R,A)
proto=struct('method','','TP_new',NaN,'TP_archived',NaN,'FP_new',NaN, ...
    'FP_archived',NaN,'GOSPA_new_m',NaN,'GOSPA_archived_m',NaN, ...
    'RMS_new_nT',NaN,'RMS_archived_nT',NaN,'same_positions',false);
Comparison=repmat(proto,3,1);
fprintf('Fresh / archived comparison (different results are kept):\n');
for j=1:3
    n=R.Methods{j};a=A.Methods{j};assert(strcmp(n.name,a.name));
    q=proto;q.method=n.name;q.TP_new=n.Metrics.TP;q.TP_archived=a.Metrics.TP;
    q.FP_new=n.Metrics.FP;q.FP_archived=a.Metrics.FP;
    q.GOSPA_new_m=n.Metrics.GOSPA;q.GOSPA_archived_m=a.Metrics.GOSPA;
    q.RMS_new_nT=n.Evaluation.observed_residual_RMSE_nT;
    q.RMS_archived_nT=a.Evaluation.observed_residual_RMSE_nT;
    q.same_positions=isequaln(n.Fit.position,a.Fit.position);Comparison(j)=q;
    fprintf('  %s: TP %d/%d; FP %d/%d; RMS %.6g/%.6g nT\n', ...
        q.method,q.TP_new,q.TP_archived,q.FP_new,q.FP_archived,q.RMS_new_nT,q.RMS_archived_nT);
end
end
