library(ASCAT)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option("--tumor-bam", type = "character", help = "Path to the tumor BAM file", metavar = "character"),
  make_option("--tumor-name", type = "character", help = "Name of the tumor sample", metavar = "character"),
  make_option("--normal-bam", type = "character", help = "Path to the normal BAM file (optional)", metavar = "character"),
  make_option("--normal-name", type = "character", help = "Name of the normal sample (optional)", metavar = "character"),
  make_option("--allelecounter-exe", type = "character", help = "Path to the allelecounter executable", metavar = "character"),
  make_option("--alleles-prefix", type = "character", help = "Prefix for alleles file", metavar = "character"),
  make_option("--loci-prefix", type = "character", help = "Prefix for loci file", metavar = "character"),
  make_option("--gender", type = "character", help = "Gender of the sample (e.g., XX or XY)", metavar = "character"),
  make_option("--genome-version", type = "character", help = "Genome version (e.g., hg38 or hg19)", metavar = "character"),
  make_option("--nthreads", type = "integer", help = "Number of threads to use", metavar = "integer"),
  make_option("--gc-content-file", type = "character", help = "Path to GC content file", metavar = "character"),
  make_option("--replic-timing-file", type = "character", help = "Path to replication timing file", metavar = "character"),
  make_option("--output-dir", type = "character", help = "Directory for output images and files", metavar = "character"),
  make_option("--tumoronly", type = "logical", default = TRUE, help = "Set to TRUE for tumor-only analysis, FALSE for tumor-and-normal analysis", metavar = "logical")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Convert input file paths to absolute paths
opt$`tumor-bam` <- normalizePath(opt$`tumor-bam`)
opt$`alleles-prefix` <- normalizePath(opt$`alleles-prefix`, mustWork = FALSE)
opt$`loci-prefix` <- normalizePath(opt$`loci-prefix`, mustWork = FALSE)
opt$`gc-content-file` <- normalizePath(opt$`gc-content-file`)
opt$`replic-timing-file` <- normalizePath(opt$`replic-timing-file`)
if (!opt$tumoronly) {
    opt$`normal-bam` <- normalizePath(opt$`normal-bam`)
}

# Map "male" and "female" to "XY" and "XX"
if (opt$gender == "male") {
  opt$gender <- "XY"
} else if (opt$gender == "female") {
  opt$gender <- "XX"
}

# Ensure the output directory exists
if (!dir.exists(opt$`output-dir`)) {
  dir.create(opt$`output-dir`, recursive = TRUE)
}

# Temporarily change the working directory to the output directory
original_wd <- getwd()
setwd(opt$`output-dir`)

# Prepare ASCAT input
if (opt$tumoronly) {
  # Tumor-only analysis
  ascat.prepareHTS(
    tumourseqfile = opt$`tumor-bam`,
    tumourname = opt$`tumor-name`,
    allelecounter_exe = opt$`allelecounter-exe`,
    alleles.prefix = opt$`alleles-prefix`,
    loci.prefix = opt$`loci-prefix`,
    gender = opt$gender,
    genomeVersion = opt$`genome-version`,
    nthreads = opt$nthreads,
    tumourLogR_file = "tumor_logr.txt",
    tumourBAF_file = "tumor_baf.txt"
  )

  ascat.bc <- ascat.loadData(
    Tumor_LogR_file = "tumor_logr.txt",
    Tumor_BAF_file = "tumor_baf.txt",
    gender = opt$gender,
    genomeVersion = opt$`genome-version`
  )

} else {
  # Tumor-and-normal analysis
  ascat.prepareHTS(
    tumourseqfile = opt$`tumor-bam`,
    normalseqfile = opt$`normal-bam`,
    tumourname = opt$`tumor-name`,  
    normalname = opt$`normal-name`,
    allelecounter_exe = opt$`allelecounter-exe`,
    alleles.prefix = opt$`alleles-prefix`,
    loci.prefix = opt$`loci-prefix`,
    gender = opt$gender,
    genomeVersion = opt$`genome-version`,
    nthreads = opt$nthreads,
    tumourLogR_file = "tumor_logr.txt",
    tumourBAF_file = "tumor_baf.txt",
    normalLogR_file = "normal_logr.txt",
    normalBAF_file = "normal_baf.txt"
  )

  ascat.bc <- ascat.loadData(
    Tumor_LogR_file = "tumor_logr.txt",
    Tumor_BAF_file = "tumor_baf.txt",
    Germline_LogR_file = "normal_logr.txt",
    Germline_BAF_file = "normal_baf.txt",
    gender = opt$gender,
    genomeVersion = opt$`genome-version`
  )

}

# Common steps for both tumor-only and tumor-and-normal
ascat.bc <- ascat.correctLogR(
  ascat.bc,
  GCcontentfile = opt$`gc-content-file`,
  replictimingfile = opt$`replic-timing-file`
)
if (opt$tumoronly) {
  gg <- ascat.predictGermlineGenotypes(ascat.bc, platform = "WGS_hg38_50X")
  ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg = gg)
} else {
  ascat.bc <- ascat.aspcf(ascat.bc)
}

# Restore the original working directory
setwd(original_wd)

ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, write_segments = TRUE, img.dir = opt$`output-dir`)

stats <- ascat.metrics(ascat.bc, ascat.output)
stats_output_file <- file.path(opt$`output-dir`, paste0(opt$`tumor-name`, "_ascat_stats.tsv"))
write.table(t(stats), stats_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

save(ascat.bc, ascat.output, stats, file = file.path(opt$`output-dir`, paste0(opt$`tumor-name`, "_ascat_bc.Rdata")))
