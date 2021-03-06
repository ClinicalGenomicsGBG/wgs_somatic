#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 40

export LD_PRELOAD=/usr/lib/petalink.so
export PETASUITE_REFPATH=/seqstore/software/petagene/corpus:/opt/petagene/petasuite/species

#find args.outputdir -name "*.bam*” -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;

OUTPUTDIR=$(echo "$1")

find $OUTPUTDIR -name "*.bam" -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;

#echo $OUTPUTDIR
#find $OUTPUTDIR -name "*.bam" -exec echo {} \;
