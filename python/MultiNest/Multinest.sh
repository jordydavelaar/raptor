#! /bin/bash

#SBATCH -N 1 -n 24

#SBATCH -p proj_bhc
#sbatch --mem=64G

##MAIL SETTINGS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=renze.oosterhuis@ru.nl

#SBATCH --time=200:00:00
python3 -u Multinest.py antoon4
