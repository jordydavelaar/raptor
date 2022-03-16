# RAPTOR++

This is the developer branch of RAPTOR.
Contributors; Thomas Bronzwaer, Jordy Davelaar, Ziri Younsi.

Code related papers;

Bronzwaer, Davelaar et al. 2018, A&A, 613, A2

Davelaar, Bronzwaer et al. 2018, CompAC, 5, 1, 1

Bronzwaer, Younsi, Davelaar et al. 2020, 641, A126

# Installation

RAPTOR needs hdf5, gsl, and openmp libraries.

RAPTOR is compilied via a makefile. This makefile needs the enviroment variable "RAPTOR" which can be set via,

```
export RAPTOR=/path/to/code
```

It is recommended to add this to your bash/zsh profiles.

To setup a RAPTOR run, first create a run directory. Copy setup.sh to the run directory. Run setup.sh


```
./setup.sh <raptor-model>
```

where ```raptor-model```, is the model corresponding with the folder name in *RAPTOR/model*. This setup should only be used for initialization of a setup folder, it should not run when a model file is present. For recompiling the code used ```make all```.   

# Running RAPTOR

RAPTOR run command is given by

```
./RAPTOR model.in <path/to/grmhd/file> output-index
```

```  model.in ```  contains model parameters. See below for an explanation

``` Inclination ``` set the inclination between  the line of sight of the observer and the angular momentum axis of the black hole and is in degrees.

``` output-index ``` is an integer that allows the user to set the output index of the writen hdf5 file.

# Model file

The model.in file allows to pass on code specific variables that are not needed during complication. This allows some flexibility in that the code does not have to be recompiled if one of these variables are changes.

The model file contains the following parameters

Model dependent parameters;

```MBH``` - Black hole mass

```M_UNIT``` - Mass scaling

```R_HIGH``` - Temperature ratio parameter based on Moscibrodzka et al. 2016

```R_LOW``` - Temperature ratio parameter based on Moscibrodzka et al. 2016


Camera dependent parameters

```INCLINATION``` - Observer inclination

```IMG_WIDTH``` - Amount of pixels in the x direction

```IMG_HEIGHT``` - Amount of pixels in the y direction

```CAM_SIZE_X``` - The FOV in the x direction

```CAM_SIZE_Y``` - The FOV in the y direction

```NUM_FREQ``` - Numbers of frequency in total

```FREQS_PER_DEC``` - The amount of frequencies per logaritmic decade

```FREQ_MIN``` - Starting frequencies of the frequency array

```STEPSIZE``` - Stepsize scaling

```MAX_LEVEL``` - Amount of adaptive levels allowed on the image domain

# Output

The output consist of an hdf5 file containing the images at all stokes parameters at all frequencies. Python scripts are provided in the *python* directory that explain how to read in and generate figures.

<p align="center">
  <img src="docs/output_example.png" width="500" title="hover text">
</p>
