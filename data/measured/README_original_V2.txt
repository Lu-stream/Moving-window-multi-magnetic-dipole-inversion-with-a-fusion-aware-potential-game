MMT Areas data
==============

This data set contains three areas of the synthetic sample studied by de Groot
et al. in [1], where the magnetic signal is obtained by means of scanning
superconducting quantum interference device microscopy (SSM) and the grain
positions via X-ray computed tomography (microCT). The three regions are
obtained by corresponding the magnetic scan signal with the reference frame of
the magnetic grains. This is achieved by scaling, rotating and translating both
the grain positions and the scan signal. The raw data without transformation
has been published in [2].

The data set is organized in three directories for the corresponding areas.
Each of these directories contain:

    AreaX/AreaX-45-fig2MMT.txt  :: The magnetic scan signal in text format for
    AreaX/AreaX-45-fig2MMT.tif     the ARM state. A snapshot of the data is
                                   provided as a TIF image.

    AreaX/AreaX-90-fig2MMT.txt  :: The magnetic scan signal in text format for
    AreaX/AreaX-90-fig2MMT.tif     the unknown (random) magnetic state. A
                                   snapshot of the data is provided as a TIF
                                   image.

    AreaX/FWInput-FineCuboids-AX.txt :: The parameters of the aggregation of
                                        cuboids making the magnetic grains, in
                                        text format. Cuboids are specified with
                                        dimensions 2a, 2b and 2c in the x, y
                                        and z direction, respectively. Then the
                                        numerical data files contain 7 columns:

                                            xc  yc  zc  a  b  c  index

                                        where:
                                        - the first three columns refer to
                                        the center of the cuboids,
                                        - columns 4-6 to the cuboid half sizes in
                                        every dimension and
                                        - the last column refers to the grain
                                        index to which cuboids belong to.

Here, X = 1,2,3, which are the labels of the areas. The scan signal is provided
in Tesla units with a step size of 1 micrometer in both the x and y directions.
The cuboid parameters are provided in micrometer units.

References
==========

[1] de Groot, L. V., Fabian, K., Béguin, A., Reith, P., Barnhoorn, A. and
Hilgenkamp, H. (2018). Determining Individual Particle Magnetizations in
Assemblages of Micrograins. Geophysical Research Letters, 45(7), 2995-3000.

[2] de Groot, L. V., Fabian, K., Béguin, A., Reith, P., Barnhoorn, A. and
Hilgenkamp, H. (2018). List of grains and results of the Scanning SQUID
Magnetometer (SSM) scan. PANGAEA, https://doi.org/10.1594/PANGAEA.886724
