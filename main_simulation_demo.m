% MAIN_SIMULATION_DEMO  Fixed-input magnetic-source inversion demonstration.
% Run this SCRIPT. Default: a fresh inversion of ONE ten-source configuration.
% The original ten layouts were selected using their results at 5% noise;
% this is a controlled illustration, not a held-out random validation sample.
%
% recompute: run candidates, local inversion, fusion/refit and game from scratch.
% archived:  load the original saved estimates and redraw; never reinvert.
% Both modes use the same frozen source configurations and Gaussian noise maps.

%% User settings
run_mode = 'recompute';             % 'recompute' or 'archived'
run_batch = false;                  % true: all 10 layouts x 5 noise levels
case_index = 4;                     % public realization_04; original ID 7
noise_fraction = 0.05;              % single-demo noise / clean-map RMS
make_figures = true;
save_figures = true;                % PNG and PDF
save_editable_fig = false;          % optional MATLAB .fig files
figure_visible = 'on';              % 'off' for unattended exports
figure_dpi = 300;
residual_color_limit_nT = 50;        % display only; [] uses full common range
model_glyph_radius_m = 2;           % display only; sources are point dipoles
model_depth_exaggeration = 5;       % display only; depth labels stay physical

% Optional unattended validation overrides; ordinary users edit values above.
v=getenv('WG_SUPPLEMENT_MODE');if ~isempty(v),run_mode=v;end
v=getenv('WG_SUPPLEMENT_VISIBLE');if ~isempty(v),figure_visible=v;end
if strcmp(getenv('WG_SUPPLEMENT_NO_PLOTS'),'1'),make_figures=false;end
if strcmp(getenv('WG_SUPPLEMENT_BATCH'),'1'),run_batch=true;end

%% Fixed experiment and output locations
package_dir = fileparts(mfilename('fullpath'));
if isempty(package_dir), package_dir = pwd; end
addpath(fullfile(package_dir,'core'),'-begin');
addpath(fullfile(package_dir,'simulation'),'-begin');
Settings = study_settings();
id_map = dlmread(fullfile(package_dir,'data','synthetic','layout_id_map.csv'),',',1,0);
assert(isequal(id_map(:,1),(1:10).'),'Invalid public layout ID map.');
layout_ids = id_map(:,1).';
noise_fraction_list = [0 0.01 0.05 0.10 0.15];
box_noise_fraction = 0.05;
assert(any(strcmp(run_mode,{'recompute','archived'})), 'Unknown run_mode.');
assert(isscalar(case_index) && case_index==floor(case_index) && ...
    case_index>=1 && case_index<=numel(layout_ids));
if ~run_batch
    layout_ids = layout_ids(case_index);
    noise_fraction_list = noise_fraction;
end
% tempname creates a unique directory name; existing outputs are never reused.
output_parent = fullfile(package_dir,'output','simulation',run_mode);
if ~exist(output_parent,'dir'), mkdir(output_parent); end
output_dir = tempname(output_parent); mkdir(output_dir);
RunManifest = struct('mode',run_mode,'software',version, ...
    'computer',computer,'timestamp',datestr(now,30), ...
    'layout_ids',layout_ids,'original_layout_ids',id_map(layout_ids,2).', ...
    'layout_id_map','data/synthetic/layout_id_map.csv', ...
    'noise_fractions',noise_fraction_list, ...
    'selection_note','Layouts were selected using the original 5% results.', ...
    'source_archive','WindowGame_TwoMains_REPACK.zip');
fprintf('\nMode: %s | %d layout(s) x %d noise level(s).\n', ...
    run_mode,numel(layout_ids),numel(noise_fraction_list));
fprintf('Results are written to %s\n',output_dir);
if strcmp(run_mode,'recompute')
    fprintf('Fresh inversion: no saved estimates or candidate caches are reused.\n');
else
    fprintf('Archived mode: saved estimates are loaded; no inversion is performed.\n');
end

%% Run every requested case in one explicitly chosen mode
Results = cell(numel(layout_ids),numel(noise_fraction_list));
RunInfo = cell(size(Results));
for il = 1:numel(layout_ids)
    for in = 1:numel(noise_fraction_list)
        [Results{il,in},RunInfo{il,in}] = supplement_simulation_case( ...
            package_dir,layout_ids(il),noise_fraction_list(in),Settings,run_mode);
        Result = Results{il,in}; %#ok<NASGU>
        save(fullfile(output_dir,[Result.Meta.name '.mat']), ...
            'Result','Settings','RunManifest','-v7');
    end
end
[Data,Summary,CaseMetrics] = study_summarize( ...
    Results,layout_ids,noise_fraction_list,output_dir);
save(fullfile(output_dir,'Summary.mat'), ...
    'Data','Summary','CaseMetrics','RunInfo','Settings','RunManifest','-v7');

%% Figures use only the selected run's estimates; they never refit
Figures = struct();
PlotOptions = struct('visible',figure_visible,'save_figures',save_figures, ...
    'save_editable_fig',save_editable_fig,'font_name','Times New Roman', ...
    'dpi',figure_dpi,'colors',[.36 .44 .57;.84 .55 .19;.03 .49 .55], ...
    'model_glyph_radius_m',model_glyph_radius_m, ...
    'model_depth_exaggeration',model_depth_exaggeration, ...
    'residual_color_limit_nT',residual_color_limit_nT, ...
    'box_noise_fraction',box_noise_fraction,'show_sd',true,'pair_lines',false, ...
    'show_titles',false);
if run_batch
    if make_figures
        Figures = study_plot_ensemble(Data,Summary,PlotOptions,output_dir);
    end
else
    Result = Results{1};
    Source = Result.Scene.Truth;
    [Metrics,FieldMaps] = study_prepare_case(Result);
    if make_figures
        Figures = study_plot_single(Result,Metrics,FieldMaps,PlotOptions,output_dir);
    end
    fprintf('\n%-23s TP FP FN Recall Precision F1 GOSPA(m) RMS(nT)\n','Method');
    for j=1:numel(Metrics)
        q=Metrics(j);
        fprintf('%-23s %2d %2d %2d %.3f %.3f %.3f %.6g %.6g\n', ...
            q.method_name,q.TP,q.FP,q.FN,q.Recall,q.Precision,q.F1, ...
            q.GOSPA_m,q.residual_RMSE_nT);
    end
end
fprintf('\nFinished. Numerical outputs and settings remain in the workspace.\n');
fprintf('Output folder: %s\n',output_dir);
