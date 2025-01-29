#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 1

SAMPLE="$1"
VCF="$2"
SIZE="$3"
REF="$4"

conda deactivate
module load anaconda2
source activate wopr_alissa

python workflows/scripts/Alissa_API_tools/alissa_API_tools.py -a $SAMPLE -v $VCF -s $SIZE -n $SAMPLE -i production -ref $REF

