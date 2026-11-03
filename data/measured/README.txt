AREA 1 MEASURED INPUTS
=====================

Area1-90-fig2MMT.txt is the registered normal-field map in the unknown/random
magnetic state from the de Groot et al. sample. This is scanning SQUID
microscopy (SSM) data, not quantum diamond microscope (QDM) data.

Acquisition: 201 rows by 351 columns; 1 micrometre pixel spacing. The text
matrix is in tesla. The MATLAB loader converts tesla to nanotesla.
Inference grid: x = -175:175, y = -100:100 in numerical micrometres.
Reference sensor-to-surface distance: 2 micrometres.

FWInput-FineCuboids-A1.txt describes registered CT cuboids in micrometres:
xc yc zc a b c grain_ID, where a, b, c are HALF sizes. Volume is 8*a*b*c.
Volume-weighted CT centres are recomputed from this file for evaluation.
CT coordinates are not passed to the analytical or local source inversion.

Underlying study:
de Groot, L. V., Fabian, K., Beguin, A., Reith, P., Barnhoorn, A., and
Hilgenkamp, H. (2018). Determining Individual Particle Magnetizations in
Assemblages of Micrograins. Geophysical Research Letters, 45, 2995-3000.
https://doi.org/10.1002/2017GL076634

Underlying raw dataset:
de Groot et al. (2018). List of grains and results of the Scanning SQUID
Magnetometer (SSM) scan. PANGAEA.
https://doi.org/10.1594/PANGAEA.886724

The supplied files are the registered Area 1 derivative already used in the
project. README_original_V2.txt documents the registration and original
attribution; keep that file with the data. This package does not claim that
the registered matrix is byte-identical to the untransformed PANGAEA scan.

area1_fixed_candidates.mat is a historical algorithm intermediate, not a
second measurement. It supplies the retained local candidates and the
Primin baseline for the explicit fixed_candidates mode. Settings are checked
against area1_candidate_settings.m before this cache can be used. Recompute
mode does not read the fixed candidate file.

historical_candidate_provenance.json is the retained provenance record from
the earlier candidate/results package. Its helper hashes describe that
historical package, not the newly organized supplement. Likewise,
archive/measured/historical_reproducibility.json describes the earlier
GNU Octave 8.4.0 fixed-candidate run, not a new execution of these entry scripts.

Important unit convention: the shared kernel receives numerical lengths in
micrometres in this experiment, while the field is in nanotesla. Its internal
moment coordinates are converted to physical A m^2 by multiplication by
1e-18. Legacy setting names ending in _m represent numerical micrometres
for this measured experiment. Depth is positive below the sensor; subtract
2 micrometres to report depth below the sample surface.
