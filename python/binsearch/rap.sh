#!/bin/bash

#SBATCH --job-name="find_M"
#SBATCH --output="slurm/rap.%j.%N.out"
#SBATCH --partition=cca
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH -t 12:00:00
#SBATCH --mem=870G

module load hdf5

echo $1 $2 $3

folder='/mnt/home/jdavelaar/scratch/KN/KN-data/output_a09375q0'

for((d=$1; d<$2; d+=$3))
do
    if (($d <10)) ; then
	./RAPTOR model.in $folder/data000$d.dat $d
    elif(($d <100)); then
	./RAPTOR model.in $folder/data00$d.dat $d
    elif(($d<1000)); then
	./RAPTOR model.in $folder/data0$d.dat $d
    else
	./RAPTOR model.in $folder/data$d.dat $d

    fi
done


