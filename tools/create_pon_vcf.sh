#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 20
#$ -q wgs.q

# Script to create normal vcfs to be used in PoN for tnscope tumor only

module load anaconda2
source activate wgs_somatic

#Input realigned bam file
NORMAL_REALIGN_BAM=$1
#Input recalibrated data table
RECAL_TABLE=$2


#Get basename of input file
BNAME_FILE="$(basename -- $NORMAL_REALIGN_BAM)"

#Get sample name (DNAXXX_XXX_XXX)
SAMPLE=(${BNAME_FILE//_REALIGNED/ })

#Get sample name (DNAXXX)
SNAME=(${SAMPLE//_/ })


#Use singularity
#Create recalibrated bam file
singularity exec --cleanenv --bind /medstore --bind /seqstore --bind /apps --bind /home --bind /seqstore/software/petagene/corpus:/opt/petagene/corpus --bind /apps/bio/dependencies/wgs_somatic/hg38/sentieon/:/databases /apps/bio/singularities/wgs_somatic/sentieon_peta_201911_hg38_new220503.simg /sentieon-genomics-201911/libexec/driver -t 20 -r /databases/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -i $NORMAL_REALIGN_BAM -q $RECAL_TABLE --algo QualCal -k /databases/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz /medstore/Development/wgs_somatic_pon/recalibrated_bams/$SAMPLE"_RECAL_DATA.TABLE.POST" --algo ReadWriter /medstore/Development/wgs_somatic_pon/recalibrated_bams/$SAMPLE"_RECALIBRATED.bam"

#Create normal vcf to use for PoN
singularity exec --cleanenv --bind /medstore --bind /seqstore --bind /apps --bind /home --bind /seqstore/software/petagene/corpus:/opt/petagene/corpus --bind /apps/bio/dependencies/wgs_somatic/hg38/sentieon/:/databases /apps/bio/singularities/wgs_somatic/sentieon_peta_201911_hg38_new220503.simg /sentieon-genomics-201911/libexec/driver -t 20 -r /databases/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -i /medstore/Development/wgs_somatic_pon/recalibrated_bams/$SAMPLE"_RECALIBRATED.bam" --algo TNscope --tumor_sample $SNAME /medstore/Development/wgs_somatic_pon/normal_vcfs/$SAMPLE"_NORMAL.vcf"
