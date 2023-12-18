#!/bin/bash
#SBATCH -p small
#SBATCH -J smth
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

export g16root=$HOME
export GAUSS_SCRDIR=/dev/shm
source $g16root/g16/bsd/g16.profile

g16 smth
