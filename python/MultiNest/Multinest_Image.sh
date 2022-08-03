#! /bin/bash

#SBATCH -N 10 -n 120

#SBATCH -p proj_bhc
#sbatch --mem=64G

##MAIL SETTINGS



#SBATCH --time=999:00:00
mpiexec -n 10 -npernode 1 python3 -u Multinest_CentaurusA.py image
