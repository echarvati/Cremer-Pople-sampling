#!/bin/bash

declare -a ins
declare -a jobs

for file in *.gjf
do
  point=${file%.gjf}
  cp "job.sh" "${point}.sh"

  echo ${file}
  sed -i -e "s/smth/${file}/g" "${point}.sh"
  sbatch "${point}.sh"

done


#

##  sbatch "job_${file}.sh"

