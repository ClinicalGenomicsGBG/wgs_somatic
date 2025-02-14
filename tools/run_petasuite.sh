#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1

COMPRESSED=$1
THREADS=$2

module load petasuite
petasuite --md5match -d $COMPRESSED -t $THREADS
