#!/bin/bash

# Input arguments
bcftools=$1
input_vcf=$2
genome_size=$3
min_mutant_allele_fraction=$4
min_mutant_allele_reads=$5
min_coverage=$6
output_file=$7
include_normal=$8

# Determine the FORMAT/AFDP value based on include_normalid
if [ "$include_normal" == "True" ]; then
    filt_norm="&& FORMAT/AFDP[1]>$min_coverage"  # Include both tumor and normal
else
    filt_norm=""  # Only include tumor
fi

# Run the bcftools and awk commands
$bcftools filter -i "FORMAT/AF[0]>$min_mutant_allele_fraction && FORMAT/AD[0:1]>$min_mutant_allele_reads && FORMAT/AFDP[0]>$min_coverage $filt_norm" $input_vcf | \
$bcftools stats | grep 'number of records:' | awk '{print $NF}' | \
awk -v gsize=$genome_size '{ printf "TMB\t%.2f\n", ($1 / gsize * 1e6) }' > $output_file

# Add parameters to the output file
echo -e "min_mutant_allele_fraction\t$min_mutant_allele_fraction" >> $output_file
echo -e "min_mutant_allele_reads\t$min_mutant_allele_reads"  >> $output_file
echo -e "min_coverage\t$min_coverage"  >> $output_file
echo -e "effective_genome_size\t$genome_size" >> $output_file
