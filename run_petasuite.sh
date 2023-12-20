#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1

threads=$1

module load petasuite
petasuite --md5match -d *.fasterq -t $threads
