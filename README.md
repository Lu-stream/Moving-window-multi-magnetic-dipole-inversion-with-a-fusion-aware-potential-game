# Fusion-aware window game: supplementary MATLAB package

This package contains two computation entry scripts, an archived-plotting
shortcut, the numerical functions they use, fixed synthetic inputs, the Area 1
measured inputs, and two demonstration archives.
It assembles the existing implementation; it does not introduce a new inversion
method or retune parameters to improve the reported comparison.

## Quick start

Download and extract the **entire** repository, then set MATLAB's Current Folder
to the folder containing `main_simulation_demo.m`. Open an entry script and
click Run. The computation settings are at the top of the two computation scripts.
Keep the subfolders with the scripts. Do not add an older WindowGame package to
the MATLAB path in the same session.

| Entry script | Default computation | Alternative |
|---|---|---|
| `main_simulation_demo.m` | Recompute one 10-source scene: realization 04 (original ID 7), 5% noise | Set `run_batch = true` for all 10 layouts at five noise levels |
| `main_measured_demo.m` | Rebuild candidates from the Area 1 measured field and run all three methods | Use `fixed_candidates` to rerun coordination from supplied candidates |
| `main_plot_archived_demo.m` | Redraw saved synthetic results using the simulation script's settings; default: realization 04, 5% noise | Does not perform inversion; keep `run_batch = false` in `main_simulation_demo.m` for a single case |

The two computation scripts default to `run_mode = 'recompute'`. This performs an actual
inversion and can take considerably longer than plotting. For a quick view of
the archived synthetic results, run `main_plot_archived_demo.m`. It temporarily
selects archived mode, restores its environment overrides even on error, and
uses the case/noise/plot settings in `main_simulation_demo.m`. It does not edit
that script. Alternatively, in either computation script, set:

```matlab
run_mode = 'archived';
```

Archived mode reads saved estimates; it does **not** run inversion. Missing
archives produce an error rather than triggering a different calculation.
To reduce download size, only two saved demonstrations are distributed:

| Saved result | Demonstration |
|---|---|
| `archive/synthetic/realization_04_noise_050.mat` | One 10-source synthetic layout: public ID 04, 5% noise |
| `archive/measured/area1_latest_result.mat` | Area 1 measured SQUID field, recovered sources, and CT-reference comparison |

For the synthetic archive, run `main_plot_archived_demo.m` with the existing
defaults `run_batch = false`, `case_index = 4`, and `noise_fraction = 0.05` in
`main_simulation_demo.m`. For the measured archive, set `run_mode = 'archived'`
at the top of `main_measured_demo.m` and run that script.

Other synthetic layouts and noise levels have no saved result in this reduced
distribution. The ten frozen inputs and all inversion code remain available:
use `run_mode = 'recompute'` and `run_batch = true` to calculate the full batch.
Full-batch `archived` mode is not available with these two archives.

All new outputs go into a unique subfolder under `output/`; supplied archives
are never overwritten. The output path is printed in the Command Window.

For calculation without graphics, set `make_figures = false`. To view figures
without exporting them, set `save_figures = false`. Editable `.fig` export is
off by default. No Python installation is needed to run the MATLAB scripts.

## What is compared

The three implemented methods are selected by their stored names, not array index:

1. `primin_direct`: this project's implementation of the Primin-based baseline.
2. `independent_refit`: independent local-cost selection, source fusion, and
   full-map refitting of moments and affine background at the fused positions.
3. `physical_game`: the same candidate library, fusion and refit, with sequential
   window updates evaluated on the fused full-map objective.

The original synthetic MAT files also retain an earlier graphical-game result.
It is preserved for archive integrity and is excluded from the three-method
comparison. Its array position must not be mistaken for the fusion-aware game.

The common numerical kernels are in `core/`. Synthetic parameters are in
`simulation/study_settings.m`; measured candidate parameters are in
`measured/area1_candidate_settings.m`, with global settings at the top of the
measured main. The game uses a finite candidate set and strict objective
decreases. A sweep-limit exit is distinguished from termination with no accepted
update; neither condition certifies globally optimal physical source recovery.

## Synthetic inputs and outputs

`data/synthetic/` contains ten frozen source configurations and standard-noise
arrays named `realization_01.mat` through `realization_10.mat`.
`case_index` selects the same consecutive public ID. The archived-result
filenames use these public IDs too. The public/original ID mapping is recorded
in `data/synthetic/layout_id_map.csv`.
All supplied MAT payloads are unchanged and retain their original metadata.
The main script labels newly loaded/computed results by public ID and retains
`original_realization` in the run metadata and per-case CSV. Consecutive
numbering is for file organization and does not change the selection history. There are 10 point dipoles per layout. Lengths are
in metres, fields in nT, and moments in A m^2.

The recomputed batch uses noise fractions `[0 0.01 0.05 0.10 0.15]`, with noise standard
deviation equal to the selected fraction of the clean full-map RMS. Each layout
uses the same frozen noise array scaled across the five levels. The arrays,
rather than a newly generated MATLAB random sequence, define the inputs.

**These layouts were selected using the original 5% noise results.** They form a
controlled illustration and noise-sensitivity study, not a held-out random
validation set. The five noise levels are paired within each layout. The source
count and true coordinates are used for forward simulation and evaluation,
not supplied to the window-selection algorithm.

The single example produces source-model, field, residual, recall and source
recovery figures. Recomputed batch mode produces separate metric boxplots and noise curves,
plus CSV and MAT summaries. Standard GOSPA is reported alongside recall,
precision, F1 and unmatched-output counts. Matching is one-to-one in 3D with a
2 m tolerance. GOSPA uses p=1, alpha=2 and c=2 m, with its own cost-minimizing
assignment. Residual color clipping is for display only; RMS uses the full array.


Archived mode redraws the saved estimates; recompute mode retains the newly
calculated estimates. These are distinct execution modes. A fresh nonlinear
inversion may produce different candidates or estimates, so outputs from
different runs should be compared explicitly before being combined.

## Measured Area 1 inputs and modes

The measured data are the Area 1 scanning SQUID/SSM map and CT-derived grain
geometry from de Groot and colleagues, not a QDM measurement:

<https://doi.org/10.1594/PANGAEA.886724>

`data/measured/Area1-90-fig2MMT.txt` contains a 201-by-351 field map in tesla.
The code converts it to nT. Horizontal spacing is 1 micrometre; the sensor is
2 micrometres above the sample surface. Internal lengths are numerical
micrometres, even in inherited field names ending in `_m`. Internally scaled
moments are multiplied by 1e-18 to obtain A m^2. Figures display fields in
microtesla. CT volume-weighted centers are computed from the supplied cuboids
and used only for reference matching after inference.

Keep `data/measured/FWInput-FineCuboids-A1.txt`: the measured plotting script
reads these cuboids to draw the actual CT grain geometry. The archived MAT
contains the fitted results and CT-reference centers; the input cuboids remain
necessary for displaying the full CT shapes. Archived CSV/JSON sidecar files
are not required; the main script regenerates its output tables from the MAT.

The measured main supports three explicitly different modes:

| Mode | Operations performed |
|---|---|
| `recompute` | Raw map -> analytical candidates -> local nonlinear fitting -> independent fusion/refit and fusion-aware game -> evaluation |
| `fixed_candidates` | Load the original candidate library and baseline; recompute independent fusion/refit and fusion-aware game -> evaluation |
| `archived` | Load the saved final estimates and fields, then write tables and plots |

Measured source matching uses a 10 micrometre three-dimensional tolerance.
The scripts write recall, precision, F1, standard GOSPA, residual RMS, and
source-estimate tables to the output folder. CT centers are evaluation
references; an unmatched estimated source is not, by itself, evidence of an
additional physical grain. Optional F2/asymmetric-loss columns are exploratory.

## Software and provenance

The scripts use MATLAB R2021b-compatible base-language functions such as
`fminsearch`, FFT and SVD. No Optimization, Statistics or Parallel Computing
Toolbox is required. Native Windows MATLAB execution has not been verified in this
packaging environment.

Primary source packages are `WindowGame_TwoMains_REPACK.zip` and
`MMT_MATLAB_and_Paper_Revision.zip`. Numerical kernels and settings were retained;
entry scripts, paths, explicit run modes and output handling were reorganized.
Synthetic archives and frozen inputs preserve their original bytes. Historical
metadata and candidate provenance are retained as records of their original
generation, not as a claim that all those computations were rerun for this package.

## Data attribution

de Groot, L. V.; Fabian, K.; Beguin, A.; Reith, P.; Barnhoorn, A.; Hilgenkamp, H.
(2018). List of grains and results of the Scanning SQUID Magnetometer (SSM) scan.
PANGAEA. <https://doi.org/10.1594/PANGAEA.886724>. The source dataset is distributed
under CC BY 3.0. The extracted raw numerical files are unchanged. Original
provider notes are retained in `data/measured/README_original_V2.txt`.

Related article: Determining Individual Particle Magnetizations in Assemblages
of Micrograins. Geophysical Research Letters 45, 2995-3000 (2018).
<https://doi.org/10.1002/2017GL076634>.
