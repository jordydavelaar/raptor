#!/bin/bash

# Request 24 CPU cores
#SBATCH -N 1 -n 24

# Request maximum time
#SBATCH --time=7:59:00

OMP_STACKSIZE=20m
export OMP_STACKSIZE

./RAPTOR model.in harm3d.txt 0
