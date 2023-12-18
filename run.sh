#!/bin/bash
#SBATCH --job-name=g16.scan_pt0
#SBATCH --output=_job.out
#SBATCH --error=_job.err
#SBATCH --time=7-0:0:0
#SBATCH --partition=cpu
#SBATCH --ntasks=32

for file in *.gjf
do
  sbg16 q=cpu -8 "$file"
done


