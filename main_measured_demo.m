%% MAIN_MEASURED_DEMO: three methods on the measured Area 1 SQUID map.
% Run this SCRIPT. Results and Figures remain in the workspace.
% Modes are explicit; no cache fallback occurs:
%   'recompute'        : raw Bz -> analytical/local candidates -> all methods
%   'fixed_candidates' : supplied candidate library -> fusion/refit/game
%   'archived'         : load the supplied historical output and plot only
% Recompute may take a long time. For the published historical figures use
% archived mode; a fresh nonlinear candidate search can differ by runtime.

%% User parameters
run_mode = 'recompute';
make_figures = true;
show_figures = true;
save_figures = true;          % PNG and vector PDF, in a new output folder
save_editable_fig = false;    % MATLAB .fig; unavailable in Octave
show_titles = false;
export_dpi = 300;
axis_font_size = 14;
label_font_size = 16;
border_width = 2.5;
residual_limit_uT = [];       % []: full shared range; 50: display [-50,50]

%% Inference settings; archived mode retains its own stored settings.
Options=struct('tau_multiplier',1,'svd_relative_tolerance',1e-12,...
 'max_sweeps',50,'evaluation_tolerance_um',10,'tolerance_grid_um',[1 2 5 10]);
package_dir=fileparts(mfilename('fullpath'));
if isempty(package_dir),package_dir=pwd;end
addpath(fullfile(package_dir,'core'),'-begin');
addpath(fullfile(package_dir,'measured'),'-begin');
input_dir=fullfile(package_dir,'data','measured');
% Optional noninteractive overrides, also used by the verification checks.
env_mode=getenv('AREA1_RUN_MODE');
if ~isempty(env_mode),run_mode=env_mode;end
if strcmp(getenv('AREA1_NO_PLOTS'),'1'),make_figures=false;end
assert(any(strcmp(run_mode,{'recompute','fixed_candidates','archived'})),...
 'run_mode must be recompute, fixed_candidates, or archived.');
output_parent=fullfile(package_dir,'output','measured');
if ~exist(output_parent,'dir'),mkdir(output_parent);end
output_dir=tempname(output_parent);mkdir(output_dir);

if strcmp(run_mode,'archived')
 result_file=fullfile(package_dir,'archive','measured','area1_latest_result.mat');
 assert(exist(result_file,'file')==2,'The archived result is missing: %s',result_file);
 saved=load(result_file,'Results');Results=saved.Results;
 fprintf('ARCHIVED output loaded; no inversion or candidate search is performed.\n');
 fprintf('Archive: %s\n',result_file);
else
 candidate_mode='rebuild_candidates';
 if strcmp(run_mode,'fixed_candidates'),candidate_mode='fixed_candidates';end
 Results=area1_run_comparison(input_dir,output_dir,candidate_mode,Options);
end

% Identify methods explicitly. Never infer their identities from a legacy
% four-method array or silently replace a newly computed result with a cache.
expected={'primin_direct','independent_refit','physical_game'};
names=cellfun(@(m)m.name,Results.Methods,'UniformOutput',false);
assert(isequal(names(:),expected(:)), 'Unexpected method set or order.');
RunInfo=struct('mode',run_mode,'recomputed',strcmp(run_mode,'recompute'),...
 'computed_with_fixed_candidates',strcmp(run_mode,'fixed_candidates'),...
 'archived_plot_only',strcmp(run_mode,'archived'),'runtime_version',version,...
 'timestamp_local',datestr(now,30),'output_dir',output_dir);
save(fullfile(output_dir,'run_result.mat'),'Results','RunInfo','-v7');
fid=fopen(fullfile(output_dir,'run_mode.txt'),'w');
assert(fid>=0,'Cannot write run provenance.');
fprintf(fid,'Mode: %s\nMATLAB/Octave version: %s\n',run_mode,version);
fprintf(fid,'Archived numerical output is never overwritten.\n');fclose(fid);
area1_write_supplementary(Results,output_dir);
fid=fopen(fullfile(output_dir,'summary.csv'),'w');
assert(fid>=0,'Cannot write summary.csv.');
fprintf(fid,'method,outputs,CT_matches,unmatched_outputs,unmatched_CT,recall,precision,F1,GOSPA_um,mean_matched_error_um,RMSE_nT,relative_residual\n');
for j=1:3
 fprintf(fid,'%s',Results.Methods{j}.name);
 fprintf(fid,',%.12g',Results.Summary(j,:));fprintf(fid,'\n');
end
fclose(fid);

Figures=struct();
if make_figures
 visible='off';if show_figures,visible='on';end
 Plot=struct('visible',visible,'show_titles',show_titles,...
  'save_figures',save_figures,'save_editable_fig',save_editable_fig,...
  'export_dpi',export_dpi,'font_name','Times New Roman',...
  'axis_font_size',axis_font_size,'label_font_size',label_font_size,...
  'border_width',border_width,'residual_limit_uT',residual_limit_uT,...
  'input_dir',input_dir,'make_counts',true);
 Figures=plot_area1_paper(Results,output_dir,Plot);
end
disp('Columns: outputs, CT matches, unmatched outputs, missed CT, recall, precision, F1, GOSPA (um), matched error (um), RMS (nT), relative residual');
disp(Results.Summary);
fprintf('Results, RunInfo and Figures remain in the workspace.\nOutput: %s\n',output_dir);
