library(ASCAT)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define variables
tumorid = args[1]
tumorLogR = args[2]
tumorBAF = args[3]
normalLogR = args[4]
normalBAF = args[5]
GCcontentfile = args[6]
replictimingfile = args[7]
output = dirname(args[8])
genomeVersion = args[9]

gender = "XX"

ascat.synchroniseFiles(samplename = tumorid, tumourLogR_file = tumorLogR, 
    tumourBAF_file = tumorBAF, normalLogR_file = normalLogR, 
    normalBAF_file = normalBAF)
ascat.bc = ascat.loadData(
  Tumor_LogR_file = tumorLogR, 
  Tumor_BAF_file = tumorBAF, 
  Germline_LogR_file = normalLogR, 
  Germline_BAF_file = normalBAF, 
  gender = gender, 
  genomeVersion = genomeVersion
  )
ascat.plotRawData(ascat.bc, img.dir = output, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(
  ascat.bc, 
  GCcontentfile = GCcontentfile, 
  replictimingfile = replictimingfile
  )
ascat.plotRawData(ascat.bc,  img.dir = output, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc, img.dir = output)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T, img.dir = output)
stats = ascat.metrics(ascat.bc,ascat.output)
stats_output_file <- file.path(output, paste0(tumorid, "_ascat_stats.tsv"))
write.table(t(stats), stats_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
