# The WGS Somatic Pipeline


## General description

 Pipeline started out as a giant shell-script developed specifically for analysing neuroblastoma WholeGenomeSequenced samples.

 Pipeline was converted to Snakemake and updated in various ways, including a transition from hg19 to hg38 as a reference genome in conjunction with the project from BarncancerFonden in association with Genomic Medicine Sweden. Clinical Genomics Göteborg has had a bioinformatician employed in this project at 80% (of FTE).

 The pipeline takes as input-data fastqfiles from a TUMOR and a NORMAL sample and generates a group of result and QC-files. The results contain Somatic AND Germline variantcalls of SNVs and InDels, as well as Structural Variants (SVs) and Copy Number Variants (CNVs) -- to provide the possibility of discovering both cancer pre-disposition variants but with a primary focus towards acquired mutations.

### How to install:

1. Clone the repository

`$ git clone https://github.com/ClinicalGenomicsGBG/wgs_somatic`

2. Install submodules

`$ git submodule update --init --recursive --remote`

Submodules used are annotate\_manta\_canvas, b\_allele\_igv\_plot, canvas\_to\_interpreter and Alissa\_API\_tools. They can all be found at [CGG](https://github.com/ClinicalGenomicsGBG).



You shouldn't have to build singularity images since paths to them are specified in the configs but if you would like to build them, you can use the definition files in the singularity subfolders.

Singularity images are located here: /apps/bio/singularities/wgs\_somatic/



### How to run manually:

1. Start a screen

2. Load anaconda2

`$ module load anaconda2`

3. Activate conda environment

`$ source activate wgs_somatic`

4. Run launch\_snakemake.py


```
$ ./launch_snakemake.py \
    --runnormal <ID of sequencing run> \
    --runtumor <ID of sequencing run> \
    --outputdir <Output directory> \
    --normalsample <Name of normal sample> \
    --normalfastqs <path to directory containing normal fastqs> \
    --tumorsample <Name of tumor sample> \
    --tumorfastqs <path to directory containing tumor fastqs> \
    --igvuser <name of igv user> \
    --hg38ref yes
```

For hg38ref, write 'yes' if you want this option. If you want to use hg19, simply don't use hg38ref argument.


If you want to run pipeline for normal only (run only germline steps of pipeline), simply don't use arguments runtumor, tumorsample and tumorfastqs.


If you want to run pipeline for tumor only, simply don't use arguments runnormal, normalsample and normalfastqs.


Runnormal and runtumor is only used to create a unique samplename based on the sequencing run the data comes from. Could probably be done in a better way.


Optional arguments to launch_snakemake.py:

```--noalissa``` Disables automatic upload of germline SNV_CNV vcf to Alissa. 

```--copyresults``` Use this argument if you want to automatically copy result files to result directory on seqstore webfolders.



Need to be on medair because python environment is hard-coded in launch_snakemake.py, also some tools in the snakemake workflow is hardcoded to install-locations on medair, and finally some dependencies such as references and databases are on medair.

 > \#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python


 Or you can just use an already set-up repository! Such as here:
 `/apps/bio/repos/wgs_somatic`


 ### Goal:


 * Have a high quality analysis of all types of genetic abberations (starting with SNVs and InDels, CNVs and SVs).
 * Which is fast (because some of these samples can help treatment of urgent pediatric cases).
 * Automated and connected to hospital systems so that it is not relied upon bioinformatician working-hours.
 * Packaged into Singularities so that it can in the future be run on a cloud-platform like AWS.
 * Connected to HCP to upload results (and download data for analysis?).
 * Robust and well-documented to fit into the clinical requirements.


 ### Yearly statistics

After running the pipeline for the first time, a yearly\_statistics text file is created in the repo. Every time the pipeline is run, sample name and date/time is added to this text file.


 ### Automatic start of pipeline

The pipeline is started automatically when new runs with GMS-BT/AL samples appear in novaseq_687_gc or novaseq_A01736 Demultiplexdirs.

Cron runs every 30 minutes (in crontab of cronuser)

Wrapper script ```wgs_somatic-run-wrapper.py``` looks for runs in Demultiplexdir. Every time there is a new run in Demultiplexdir, it is added to text file ```/apps/bio/repos/wgs_somatic/novaseq_runlist.txt``` to keep track of which runs that have already been analyzed. If a new run has GMS-BT/AL samples, the pipeline starts for these samples. Output is placed in working directory ```/medstore/results/wgs/wgs_somatic``` and the final result files are then copied to seqstore webfolders. You can find a more detailed description of the automation on GMS-BT confluence page.


 ### Status (2020-10-19):

 #### Reference genome

Can be run using either hg38 or hg19 reference genomes.

 #### Validation and optimisation:

 ##### Germline calling:
 * SNVs and InDel calling with DNAscope (was validated through WOPR in hg19, but not on hg38)
 * CNVs -- Canvas
 * SVs -- Manta


 ##### Somatic calling:
 * SNVs and InDel calling with TNscope and custom Bcftools filters. Optimised against artificial truthset by mixing coriell samples (see wiki)
 * CNVs -- Canvas
 * SVs -- Manta

 #### Code structure and organisation of dependencies

 Most of the pipeline is packaged in Singularities, Sentieon and Canvas (Sentieon is majority of pipeline).
 The remaining rules relies primarily on small python scripts, bcftools and Manta. Should not be too difficult to package these into singularities as well.

 #### Validation project 

 * Analyse all samples in the project with hg38 (nearly done)
 * Calculate sensitivity and precision with artificial truthset (nearly done)


