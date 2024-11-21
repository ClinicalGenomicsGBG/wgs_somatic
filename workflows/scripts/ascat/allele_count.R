library(ASCAT)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define variables
seqfile = args[1]
name = args[2]
CHR = args[3]
output = args[4]
loci.prefix = args[5]
allelecounter_exe = args[6]

min_base_qual = 20
min_map_qual = 35
additional_allelecounter_flags = NA

if (!dir.exists(dirname(output))) {
  dir.create(dirname(output))
}

cat("Running allele counting for ", name, " chromosome ", CHR, "\n")
cat("Output will be written to ", output, "\n")
cat("Using loci file ", paste0(loci.prefix, CHR, ".txt"), "\n")


ascat.getAlleleCounts(
  seq.file = seqfile, 
  output.file = output, 
  loci.file = paste0(loci.prefix, CHR, ".txt"), 
  min.base.qual = min_base_qual, 
  min.map.qual = min_map_qual, 
  allelecounter.exe = allelecounter_exe, 
  additional_allelecounter_flags = additional_allelecounter_flags)