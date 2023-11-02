#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 20


module load petasuite
petasuite --md5match -d *.fasterq -t 20
