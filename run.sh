#!/bin/bash

# Request 24 CPU cores
#SBATCH -N 1 -n 24

# Request maximum time
#SBATCH --time=10:00:00

OMP_STACKSIZE=20m
export OMP_STACKSIZE

./RAPTOR model.in grmhd/dump019 75 0


