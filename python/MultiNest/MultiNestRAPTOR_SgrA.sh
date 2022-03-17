#! /bin/bash
# Request 24 CPU cores
#SBATCH -N 1 -n 24

# Request maximum time
#SBATCH --time=720:00:00

#SBATCH -p proj_bhc
#sbatch --mem=64G

##MAIL SETTINGS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.dekleuver@science.ru.nl

python MultiNestRAPTOR_SgrA.py $1
