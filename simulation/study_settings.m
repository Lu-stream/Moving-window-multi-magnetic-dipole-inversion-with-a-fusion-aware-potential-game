function Settings = study_settings()
% Original controlled numerical settings. Shared by both main scripts.
% Keep unchanged when reproducing the supplied fixed-layout experiment.
source_count = 10; % forward simulation/evaluation ONLY
original_geometry_seed = 2026102001;
original_noise_seed = 2026102101;
geometry_seed_base = 2026102500;
noise_seed_base = 2026102600;
use_frozen_inputs = true; % identical source/noise arrays in MATLAB and Octave
Grid=struct(); SourceOptions=struct(); Candidate=struct(); Primin=struct();
Local=struct(); Game=struct(); Fusion=struct(); Physical=struct(); Evaluation=struct();
Grid.n = 145;
Grid.scan_half_width_m = 72;
SourceOptions.layout = 'volume';
SourceOptions.source_bounds_xy = [-50 -50;50 50];
SourceOptions.depth_range_m = [3 15];
SourceOptions.moment_range_Am2 = [800 1200];
SourceOptions.minimum_separation_m = 4; % XYZ, not horizontal separation
SourceOptions.direction_tilt_deg = 20; % uniform solid angle in cone
SourceOptions.max_placement_attempts = 100000;

Candidate.center_m = [0 0];
Candidate.window_widths_m = 128;
Candidate.position_bounds = [-64 -64 .5;64 64 20];
Candidate.candidate_stride = 1;
Candidate.min_b_fraction = .01;
Candidate.relative_field_floor = .02;
Candidate.pad_factor = 4;
Candidate.patch_radius_pixels = 4;
Primin.patch_radius_pixels = 4;
Primin.batch_size = 128;
Primin.max_score_elements = 2000000;
Primin.cluster_space = 'xyz';
Primin.cluster_epsilon_m = .75;
Primin.cluster_min_points = 10;

Local.half_width_pixels = 10;
Local.stride_pixels = 10;
Local.max_sources = 2;
Local.max_hypotheses_per_order = 2;
Local.starts_per_order = 2;
Local.max_iterations = 1000; % current controlled experiment budget
Local.max_function_evaluations = 4000; % historical main used 200/800
Local.position_bounds = Candidate.position_bounds;
Local.minimum_source_separation_m = 2;
Local.model_penalty = .015;
Local.residual_reference = .15;
Local.relative_field_floor = .02;
Local.signal_fraction_floor = .03;
Local.background_order = 1;
Local.moment_ridge = 1e-8;
Local.max_condition = 1e10;
Local.seed_merge_distance_m = 2;
Local.max_seed_candidates = 8;
Local.minimum_component_fraction = 0;
Local.verbose = true;
Game.lambda = .15;
Game.position_scale_m = 2;
Game.moment_scale = .50;
Game.max_sweeps = 50;
Game.improvement_tolerance = 1e-10;
Game.min_overlap_fraction = .15;
Game.missing_penalty = .25;
Game.order = 'forward';
Game.verbose = true;
Fusion.position_radius_m = 2;
Fusion.moment_relative_gate = .75;
Fusion.minimum_window_support = 2;
Fusion.max_cluster_diameter_m = 4;
Fusion.core_margin_fraction = .10;
Fusion.minimum_core_window_support = 1;
Fusion.refinement = 'none';
Physical.tau_multiplier = 1; % tau = 3*log(number of observations)
Physical.svd_relative_tolerance = 1e-12;
Physical.version = 'fusion-aware-physical-v1';
Evaluation.position_tolerance_m = 2;


Settings=struct('Version','ten-source-physical-frozen-v2-matlab-source-rng',...
 'SourceCount',source_count,'OriginalGeometrySeed',original_geometry_seed,...
 'OriginalNoiseSeed',original_noise_seed,'GeometrySeedBase',geometry_seed_base,...
 'NoiseSeedBase',noise_seed_base,'Grid',Grid,'SourceOptions',SourceOptions,...
 'Candidate',Candidate,'Primin',Primin,'Local',Local,'Game',Game,...
 'Fusion',Fusion,'Physical',Physical,'Evaluation',Evaluation,...
 'UseFrozenInputs',use_frozen_inputs);
end
