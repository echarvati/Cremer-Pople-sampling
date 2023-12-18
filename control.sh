#!/bin/bash



for file in *.log
do
  if grep -q "Normal termination of Gaussian 16" "$file"
  then
    echo "$file" DONE

  elif grep -q "Error termination via Lnk1e" "$file"
  then
    echo "$file" ERROR

  else
    echo "$file" RUNNING
  fi
done



