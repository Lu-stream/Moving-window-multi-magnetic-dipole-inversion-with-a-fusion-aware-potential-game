% MAIN_PLOT_ARCHIVED_DEMO  Redraw saved synthetic estimates; no inversion.
% Run this SCRIPT from the repository containing main_simulation_demo.m.
% Uses the case/noise/plot settings at the top of main_simulation_demo.m.
% Distributed defaults: run_batch=false, case_index=4, noise_fraction=0.05.
% Public layout 04 is original layout 7. Keep run_batch=false for one case.
% Requires archive/synthetic/realization_04_noise_050.mat for the default.
% The saved independent-method counts are TP=8, FP=0, FN=2; these differ
% from Figure 4(b) of the submitted manuscript. See README.md for details.
% This entry preserves the real recomputation entry and never substitutes
% archived estimates into a run labelled as recompute.

wg_ad_names = {'WG_SUPPLEMENT_MODE', ...
               'WG_SUPPLEMENT_BATCH', ...
               'WG_SUPPLEMENT_NO_PLOTS'};
wg_ad_values = {'archived', '0', '0'};
wg_ad_previous = cellfun(@getenv, wg_ad_names, 'UniformOutput', false);

try
    for wg_ad_k = 1:numel(wg_ad_names)
        setenv(wg_ad_names{wg_ad_k}, wg_ad_values{wg_ad_k});
    end
    run(fullfile(fileparts(mfilename('fullpath')), 'main_simulation_demo.m'));
catch wg_ad_error
    for wg_ad_k = 1:numel(wg_ad_names)
        setenv(wg_ad_names{wg_ad_k}, wg_ad_previous{wg_ad_k});
    end
    rethrow(wg_ad_error);
end

for wg_ad_k = 1:numel(wg_ad_names)
    setenv(wg_ad_names{wg_ad_k}, wg_ad_previous{wg_ad_k});
end
clear wg_ad_names wg_ad_values wg_ad_previous wg_ad_k
