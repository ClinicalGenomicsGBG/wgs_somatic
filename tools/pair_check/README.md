# Pair check

## Determine match

In the pipeline, a pair check is done to make sure that tumor and normal samples come from the same patient. Script [determine_match.py](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/workflows/scripts/determine_match.py) is used for this.

“SNVs Only” vcf:s for tumor and normal (outputs from DNAScope) are used to determine match for pair-check. (“{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf”)
 
Tumor vcf:
For each chromosome, starting from position 200k, dbsnp labeled variants are extracted. 20k variants are extracted from one chromosome until it switches to the next chromosome. In total, 400k dbsnp labeled variants are added to a list.
 
Normal vcf:
For each chromosome, starting from position 200k, dbsnp labeled variants are extracted. 40k variants are extracted from one chromosome until it switches to the next chromosome. In total, 800k dbsnp labeled variants are added to a list.
 
Each variant in the list of tumor variants that exists in the list of normal variants is labeled as a match and each variant that does not exist is labeled as a mismatch. If the match fraction is >=0.95, the match status between the tumor and normal samples will be “match”. If the match fraction is <0.95 but >=0.90, the match status will be “warning match lower than expected”. If the match fraction is <0.85, the match status will be “error, large mismatch”.

## Validation of pair check

To validate the pair check, and the match fraction intervals for match status, script [determine_match.py](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/workflows/scripts/determine_match.py) was run on 55 known pairs of tumor/normal (all samples available up until 2021-05-31), 55 known mismatch pairs (tumor vcf from one patient and normal vcf from another patient) and 55 pairs of family members (vcf from one family member as "tumor vcf" and vcf from another family member as "normal vcf"). Match fraction information was summarized in excel files for each group. These files were used as input to script Script [plot_matchfraction_groups.R](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/plot_matchfraction_groups.R). 

Script [plot_matchfraction_groups.R](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/plot_matchfraction_groups.R) was used to plot the match fractions of the three groups as normal distributions, and the following is the resulting plot [groups_match_fraction.pdf](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/groups_match_fraction.pdf). For each group, in the plot there are lines plus/minus one, two and three standard deviations from the mean value. There is also one plot per group: [match_plot.pdf](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/match_plot.pdf), [mismatch_plot.pdf](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/mismatch_plot.pdf) and [family_plot.pdf](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/pair_check/family_plot.pdf).

Previously, the interval for match status "warning match lower than expected" was match fraction <0.95 but >=0.85. As of 210621, the interval for match status "warning match lower than expected" is match fraction <0.95 but >=0.90.
