Radboud Polarized Integrator v1.0
Branch for emission models


Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)

Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi

This program integrates the equations of motion of General Relativity
to compute the trajectories of photon bundles (null geodesics); it then
performs radiative transfer along these geodesics to compute an image
or spectrum. The gravitational field is defined by the metric selected
by the user; plasma models can be supplied in the form of GRMHD
simulations or analytic models.

NOTE: RAPTOR is proprietary and private! This code is not yet publicly released and may not be used without permission of the authors.

HOW TO USE RAPTOR:

1) Download source code, compile by typing "make img".

2) Create output directory "RAPTORfolder/output", where RAPTORfolder is the folder in which you downloaded the source.

3) Run by typing "./RAPTOR model.in grmhd_file obs_inclination photon_t0", where model.in is the input file (which must be formatted exactly as the example provided here), grmhd_file is the GRMHD dump file to be used, obs_inclination is the inclination of the observer, and photon_t0 is the initial time of the photon.

Example of run command:

./RAPTOR model.in grmhd/dump019 90 0

This requires a HARM2D GRMHD file named "dump019" to be present in the subfolder "RAPTORfolder/grmhd" and produces two output files in the subdirectory "RAPTORfolder/output", viz. img_data.dat (a simple text file), and img_data.vtk (a VTK compatible file).

Note: HARM2D (Gammie et al. 2003) is available in the Illinois astrophysical code library, as is the dump file "dump019". See:

http://rainman.astro.illinois.edu/codelib/

(dump019 can be found under "grmonty".)