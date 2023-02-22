#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 1

CREDENTIALS=$1
BUCKET=$2
FQ_FILE=$3
OUTDIR=$4

module load miniconda/4.8.3
source activate scruffy #iris_dev202210
iris -c $CREDENTIALS -b $BUCKET download -f $FQ_FILE -o $OUTDIR
