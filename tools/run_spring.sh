#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1

COMPRESSED=$1
DECOMPRESSED=$2
THREADS=$3

module load spring
spring -d -i $COMPRESSED -o $DECOMPRESSED -t $THREADS -g
