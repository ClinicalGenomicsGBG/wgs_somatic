library(ASCAT)
library(data.table)

# ASCAT functions
ascat.getBAFsAndLogRs_mod <- function (samplename, 
  tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, 
  out_BAF_file, out_LogR_file, stype, 
  alleles.prefix, gender, genomeVersion, chrom_names = c(1:22, 
    "X"), minCounts = 20, BED_file = NA, probloci_file = NA, 
  tumour_only_mode = FALSE, loci_binsize = 1, seed = as.integer(Sys.time())) 
{
  set.seed(seed)
  stopifnot(gender %in% c("XX", "XY"))
  stopifnot(genomeVersion %in% c("hg19", "hg38", "CHM13"))
  tumour_input_data = ASCAT:::readAlleleCountFiles(tumourAlleleCountsFile.prefix, 
    ".txt", chrom_names, 1)
  if (tumour_only_mode) {
    normal_input_data = tumour_input_data
  }
  else {
    normal_input_data = ASCAT:::readAlleleCountFiles(normalAlleleCountsFile.prefix, 
      ".txt", chrom_names, minCounts)
  }
  allele_data = ASCAT:::readAllelesFiles(alleles.prefix, ".txt", chrom_names)
  matched_data = Reduce(intersect, list(rownames(tumour_input_data), 
    rownames(normal_input_data), rownames(allele_data)))
  tumour_input_data = tumour_input_data[rownames(tumour_input_data) %in% 
    matched_data, ]
  normal_input_data = normal_input_data[rownames(normal_input_data) %in% 
    matched_data, ]
  allele_data = allele_data[rownames(allele_data) %in% matched_data, 
    ]
  rm(matched_data)
  if (!is.na(probloci_file)) {
    stopifnot(file.exists(probloci_file) && file.info(probloci_file)$size > 
      0)
    probloci = data.frame(fread(probloci_file, sep = "\t", 
      showProgress = FALSE, header = TRUE), stringsAsFactors = FALSE)
    probloci = paste0(gsub("^chr", "", probloci[, 1]), "_", 
      probloci[, 2])
    probloci = which(rownames(tumour_input_data) %in% probloci)
    if (length(probloci) > 0) {
      tumour_input_data = tumour_input_data[-probloci, 
        ]
      normal_input_data = normal_input_data[-probloci, 
        ]
      allele_data = allele_data[-probloci, ]
    }
    else {
      warning("The probloci did not remove any SNPs, it might be worth checking the data.")
    }
    rm(probloci)
  }
  stopifnot(isTRUE(all.equal(allele_data[, 1], tumour_input_data[, 
    1])) && isTRUE(all.equal(allele_data[, 1], normal_input_data[, 
    1])))
  stopifnot(isTRUE(all.equal(allele_data[, 2], tumour_input_data[, 
    2])) && isTRUE(all.equal(allele_data[, 2], normal_input_data[, 
    2])))
  tumour_input_data = tumour_input_data[, 3:6]
  normal_input_data = normal_input_data[, 3:6]
  if (!is.na(BED_file)) {
    stopifnot(file.exists(BED_file) && file.info(BED_file)$size > 
      0)
    BED = read.table(BED_file, sep = "\t", header = FALSE, 
      stringsAsFactors = FALSE)[, 1:3]
    colnames(BED) = c("chr", "start", "end")
    BED$chr = gsub("^chr", "", BED$chr)
    BED$start = BED$start + 1
    BED = BED[BED$chr %in% chrom_names, ]
    if (nrow(BED) == 0) 
      stop("Major issue with BED file, please double-check its content")
    requireNamespace("GenomicRanges")
    requireNamespace("IRanges")
    overlaps = findOverlaps(GRanges(seqnames = BED$chr, 
      ranges = IRanges(start = BED$start, end = BED$end)), 
      GRanges(seqnames = allele_data$chromosome, ranges = IRanges(start = allele_data$position, 
        end = allele_data$position)))
    if (length(overlaps) > 0) {
      tumour_input_data = tumour_input_data[unique(overlaps@to), 
        ]
      normal_input_data = normal_input_data[unique(overlaps@to), 
        ]
      allele_data = allele_data[unique(overlaps@to), ]
    }
    else {
      print(head(allele_data))
      print(head(BED))
      stop("The overlap between the BED file and loci is empty. Data must be checked!")
    }
    rm(BED, overlaps)
  }
  len = nrow(allele_data)
  tumour_input_data$REF = tumour_input_data[cbind(1:len, allele_data[, 
    3])]
  tumour_input_data$ALT = tumour_input_data[cbind(1:len, allele_data[, 
    4])]
  normal_input_data$REF = normal_input_data[cbind(1:len, allele_data[, 
    3])]
  normal_input_data$ALT = normal_input_data[cbind(1:len, allele_data[, 
    4])]
  if (tumour_only_mode) {
    TO_KEEP = which(tumour_input_data$REF + tumour_input_data$ALT >= 
      1)
  }
  else {
    TO_KEEP = which(tumour_input_data$REF + tumour_input_data$ALT >= 
      1 & normal_input_data$REF + normal_input_data$ALT >= 
      minCounts)
  }
  stopifnot(length(TO_KEEP) > 0)
  allele_data = allele_data[TO_KEEP, ]
  tumour_input_data = tumour_input_data[TO_KEEP, ]
  normal_input_data = normal_input_data[TO_KEEP, ]
  rm(TO_KEEP)
  len = nrow(allele_data)
  mutCount1 = tumour_input_data$REF
  mutCount2 = tumour_input_data$ALT
  totalTumour = mutCount1 + mutCount2
  normCount1 = normal_input_data$REF
  normCount2 = normal_input_data$ALT
  totalNormal = normCount1 + normCount2
  rm(tumour_input_data, normal_input_data)
  normalBAF = vector(length = len, mode = "numeric")
  tumourBAF = vector(length = len, mode = "numeric")
  normalLogR = vector(length = len, mode = "numeric")
  tumourLogR = vector(length = len, mode = "numeric")
  normalBAF_unmirrored = normCount2/totalNormal
  tumourBAF_unmirrored = mutCount2/totalTumour
  germline.BAF_unmirrored = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, baf = normalBAF_unmirrored, 
    ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
  tumor.BAF_unmirrored = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, baf = tumourBAF_unmirrored, 
    ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
  colnames(tumor.BAF_unmirrored)[3] = samplename
  colnames(germline.BAF_unmirrored)[3] = samplename
  if (stype == "tumor") {
    write.table(tumor.BAF_unmirrored, file = gsub("\\.txt$", 
      "_rawBAF.txt", out_BAF_file), row.names = TRUE, quote = FALSE, 
      sep = "\t", col.names = NA)
  }
  else {
    if (!tumour_only_mode) {
    write.table(germline.BAF_unmirrored, file = gsub("\\.txt$", 
      "_rawBAF.txt", out_BAF_file), row.names = TRUE, 
      quote = FALSE, sep = "\t", col.names = NA)
    }
  }
  rm(normalBAF_unmirrored, tumourBAF_unmirrored, germline.BAF_unmirrored, 
    tumor.BAF_unmirrored)
  selector = round(runif(len))
  normalBAF[which(selector == 0)] = normCount1[which(selector == 
    0)]/totalNormal[which(selector == 0)]
  normalBAF[which(selector == 1)] = normCount2[which(selector == 
    1)]/totalNormal[which(selector == 1)]
  tumourBAF[which(selector == 0)] = mutCount1[which(selector == 
    0)]/totalTumour[which(selector == 0)]
  tumourBAF[which(selector == 1)] = mutCount2[which(selector == 
    1)]/totalTumour[which(selector == 1)]
  rm(selector)
  if (tumour_only_mode) {
    tumourLogR = log2(totalTumour/median(totalTumour, na.rm = TRUE))
  }
  else {
    tumourLogR = totalTumour/totalNormal
    tumourLogR = log2(tumourLogR/mean(tumourLogR, na.rm = TRUE))
    if (gender == "XY") {
      if (genomeVersion == "hg19") {
        nonPAR = c(2699521, 154931043)
      }
      else if (genomeVersion == "hg38") {
        nonPAR = c(2781480, 155701382)
      }
      else if (genomeVersion == "CHM13") {
        nonPAR = c(2394411, 153925834)
      }
      nonPAR = which(allele_data$chromosome %in% c("X", 
        "chrX") & allele_data$position >= nonPAR[1] & 
        allele_data$position <= nonPAR[2])
      tumourLogR[nonPAR] = tumourLogR[nonPAR] - 1
    }
  }
  tumor.LogR = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, logr = tumourLogR, 
    ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
  tumor.BAF = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, baf = tumourBAF, ID = rownames(allele_data), 
    row.names = 4, stringsAsFactors = FALSE)
  germline.LogR = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, logr = normalLogR, 
    ID = rownames(allele_data), row.names = 4, stringsAsFactors = FALSE)
  germline.BAF = data.frame(Chromosome = allele_data$chromosome, 
    Position = allele_data$position, baf = normalBAF, ID = rownames(allele_data), 
    row.names = 4, stringsAsFactors = FALSE)
  colnames(tumor.LogR)[3] = samplename
  colnames(tumor.BAF)[3] = samplename
  colnames(germline.LogR)[3] = samplename
  colnames(germline.BAF)[3] = samplename
  if (loci_binsize > 1) {
    posbins <- paste(germline.BAF$Chromosome, round(germline.BAF$Position/loci_binsize), 
      sep = "_")
    hetidx <- which(germline.BAF$baf > 0.01 & germline.BAF$baf < 
      0.99)
    hetidx <- hetidx[which(!duplicated(posbins[hetidx]))]
    homidx <- which(!(duplicated(posbins) | posbins %in% 
      posbins[hetidx]))
    unifidx <- sort(x = union(homidx, hetidx), decreasing = FALSE)
  }
  else {
    unifidx <- 1:nrow(germline.BAF)
  }
  if (stype == "tumor") {
    write.table(tumor.LogR[unifidx, ], file = out_LogR_file, 
      row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
    write.table(tumor.BAF[unifidx, ], file = out_BAF_file, 
      row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
  }
  else {
    if (!tumour_only_mode) {
      write.table(germline.LogR[unifidx, ], file = out_LogR_file, 
        row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
      write.table(germline.BAF[unifidx, ], file = out_BAF_file, 
        row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA)
    }
  }
}


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define variables
tumorid = args[1]
in_tumor = args[2]
in_normal = args[3]
out_BAF = args[4]
out_LogR = args[5]
alleles.prefix = args[6]
stype = args[7]
tumoronly <- ifelse(length(args) >= 8, args[8], "FALSE")

# Convert tumoronly to logical
tumoronly <- ifelse(tumoronly == "TRUE", TRUE, FALSE)

tumor_prefix = sub("^(.*_alleleFrequencies_chr)[0-9XY]+\\.txt$", "\\1", in_tumor)
if (tumoronly) {
  normal_prefix = tumor_prefix
} else {
  normal_prefix = sub("^(.*_alleleFrequencies_chr)[0-9XY]+\\.txt$", "\\1", in_normal)
}

if (stype != "normal") {
  stype = "tumor"
}

gender = "XX"
genomeVersion = "hg38"
chrom_names = c(1:22, "X")
minCounts = 10
BED_file = NA
probloci_file = NA
loci_binsize = 1
seed = as.integer(Sys.time())

ascat.getBAFsAndLogRs_mod(samplename = tumorid, tumourAlleleCountsFile.prefix = tumor_prefix,
  normalAlleleCountsFile.prefix = normal_prefix, out_BAF_file = out_BAF, out_LogR_file = out_LogR,
  stype = stype, alleles.prefix = alleles.prefix, 
  gender = gender, genomeVersion = genomeVersion, chrom_names = chrom_names, 
  minCounts = minCounts, BED_file = BED_file, probloci_file = probloci_file, 
  tumour_only_mode = tumoronly, loci_binsize = loci_binsize, 
  seed = seed)
