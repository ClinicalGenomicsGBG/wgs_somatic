library(data.table)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option("--tumorname", type = "character", help = "Tumor sample name", metavar = "character"),
  make_option("--gender", type = "character", help = "Gender of the sample (e.g., XX or XY)", metavar = "character"),
  make_option("--genome-fai", type = "character", help = "Path to the genome FAI file", metavar = "character"),
  make_option("--segments", type = "character", help = "Path to the segments file", metavar = "character"),
  make_option("--Rdata-file", type = "character", help = "Path to the ascat run Rdata file", metavar = "character"),
  make_option("--output-plot", type = "character", help = "Path to output plot", metavar = "character"),
  make_option("--output-seg", type = "character", help = "Path to output segments file", metavar = "character"),
  make_option("--output-baf", type = "character", help = "Path to output BAF file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Map "male" and "female" to "XY" and "XX"
if (opt$gender == "male") {
  opt$gender <- "XY"
} else if (opt$gender == "female") {
  opt$gender <- "XX"
}

theme_set(theme_pubclean())


fai <- read.table(opt$`genome-fai`, header = FALSE, sep = "\t", 
                  col.names = c("chr", "length", "start", "value1", "value2")) %>%
  mutate(chr = substr(chr,4,30))%>%
  mutate(middle = start + (length / 2))
if (opt$`gender` == "XY") {
  fai <- fai[1:24,]
} else {
  fai <- fai[1:23,]
  }

## Segments
# Read the segments table
seg_df <- read.table(opt$`segments`, header = TRUE, sep = "\t")
# Merge with fai to adjust positions
seg_df <- merge(seg_df, fai, by = "chr") %>%
  mutate(
    startpos = startpos + start,
    endpos = endpos + start,
    segment_size = endpos - startpos
  )%>%
  dplyr::rename(minor_allele = nMinor, major_allele = nMajor)%>%
  pivot_longer(cols = c(major_allele, minor_allele), names_to = "allele", values_to = "ascat_ploidy")

# If there are any segments with more than 3 times the median segment ploidy, cap them to 3 times the median
max_ploidy <- median(seg_df$ascat_ploidy, na.rm = TRUE) * 3
if (any(seg_df$ascat_ploidy > max_ploidy)) {
  seg_df <- seg_df %>%
    mutate(ascat_ploidy = ifelse(
      ascat_ploidy > max_ploidy,
      # If the ploidy is greater than max_ploidy, subtract the rounded difference from ascat_ploidy
      ascat_ploidy - (round(ascat_ploidy - max_ploidy, 0)),
      ascat_ploidy
    ))
  breaks <- seq(0, max_ploidy, by = 1)
  labels <- c(as.character(breaks[-length(breaks)]), paste0(">", max_ploidy))
} else {
  breaks <- seq(0, max(seg_df$ascat_ploidy, na.rm = TRUE), by = 1)
  labels <- as.character(breaks)
}

## Load ascat.bc from the Rdata file
load(opt$`Rdata-file`)  # Load the Rdata file containing ascat.bc

## BAF
# Extract BAF from the ascat.bc object
setDT(ascat.bc$Tumor_BAF, keep.rownames = TRUE)[]%>%
  separate_wider_delim(rn, "_", names = c("chr","pos"))%>%
  mutate(pos = as.numeric(pos))%>%
  dplyr::select(1,2,"BAF"= 3)->tumorBAF_df
# Merge with fai to adjust positions
tumorBAF_df_adj <- merge(tumorBAF_df, fai, by = "chr") %>%
  mutate(
    pos = pos + start
  )

## LogR
# Extract LogR from the ascat.bc object
setDT(ascat.bc$Tumor_LogR, keep.rownames = T)[]%>%
  separate_wider_delim(rn, "_", names = c("chr","pos"))%>%
  mutate(pos = as.numeric(pos))%>%
  dplyr::select(1,2,"LogR"= 3)->tumorLogR_df
# Merge with fai to adjust positions
tumorLogR_df_adj <- merge(tumorLogR_df, fai, by = "chr") %>%
  mutate(
    pos = pos + start
  )

## Plot segments, BAF, and LogR
png(opt$`output-plot`, width = 24, height = 24, units = "cm", res = 1200)

Segment_plot <- ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_linerange(data = seg_df, 
                 aes(xmin = startpos, xmax = endpos, y = ascat_ploidy, col = allele), 
                 linewidth = 2.5, position = position_dodge(width = -0.1))+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  scale_y_continuous("Copy number",
                     breaks = breaks,
                     labels = labels,
                     expand = expansion(mult = 0.1)) +
  theme(legend.title=element_blank(), axis.title.x = element_blank())

BAF_plot <-  ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_point(data = tumorBAF_df_adj, aes(pos, BAF), shape = '.', col = "#00000010")+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9), 
                     limits = c(0.05,0.95)) +
  theme(axis.title.x = element_blank())

LogR_plot <-  ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_point(data = tumorLogR_df_adj, aes(pos, LogR), shape = '.', col = "#00000010")+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  scale_y_continuous(limits = c(-2,2))+
  theme(axis.title.x = element_blank())

plot_grid(Segment_plot, BAF_plot, LogR_plot, 
          ncol = 1, align = "v", rel_heights = c(2,1,1))

dev.off()

# Write the segmentation data to a file for IGV visualization
output_ratio_seg <- opt$`output-seg`
if (is.null(output_ratio_seg)) {
  output_ratio_seg <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_copynumber_IGV.seg"))
}

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0", con = output_ratio_seg)

# Process the segments file to output nMajor and nMinor as separate rows
read.table(opt$`segments`, header = TRUE, sep = "\t") %>%
  dplyr::rename(Sample = sample, Chromosome = chr, Start = startpos, End = endpos) %>%
  mutate(
    Chromosome = paste0("chr", Chromosome),
    Copynumber = nMajor + nMinor,
  ) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, Copynumber) %>%
  write.table(file = output_ratio_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)


# Write the BAF data to an IGV-compatible file
output_baf_seg <- opt$`output-baf`
if (is.null(output_baf_seg)) {
  output_baf_seg <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_BAF_IGV.seg"))
}

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=100:100:100 color=0,0,255 altColor=255,0,0", con = output_baf_seg)

# Process the BAF data for IGV visualization
tumorBAF_df %>%
  dplyr::rename(Chromosome = chr, Start = pos) %>%
  mutate(Chromosome = paste0("chr", Chromosome), End = Start, Sample = opt$`tumorname`) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, BAF) %>%
  write.table(file = output_baf_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)

