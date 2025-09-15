# Load required libraries
library(data.table)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(optparse)

# Function to plot ascat panels
plot_ascat_panels <- function(fai, seg_df_adj, seg_df, tumorBAF_df_adj, CNs_adj_plot, breaks, labels, chr = NULL, cytoband = NULL) {
  if (!is.null(chr)) {
    # Single chromosome plots
    # Filter data for the specified chromosome
    fai <- fai %>% filter(chr == !!chr)
    seg_df <- seg_df %>% filter(chr == !!chr)
    tumorBAF_df_adj <- tumorBAF_df_adj %>% filter(chr == !!chr)
    CNs_adj_plot <- CNs_adj_plot %>% filter(chr == !!chr)
    x_limits <- range(
      c(0, CNs_adj_plot$pos, tumorBAF_df_adj$pos, seg_df$endpos),
      na.rm = TRUE
    )
    if (!is.null(cytoband)) {
      cytoband <- cytoband %>% filter(chrom == !!chr)
    }

    # Use unadjusted seg_df for the chromosome-specific plot
    Segment_plot <- ggplot(fai) +
    geom_point(data = dplyr::mutate(CNs_adj_plot, track = "CN_obs"), aes(x = pos, y = CN_smooth, col = track), shape = 20) +
    geom_linerange(data = seg_df, 
                  aes(xmin = startpos, xmax = endpos, y = ascat_ploidy, col = allele), 
                  linewidth = 2.5, position = position_dodge(width = -0.1)) +
    geom_point(data = dplyr::mutate(CNs_adj_plot, track = "CN_call"), aes(x = pos, y = CN_seg, col = track), shape = 20) +
    scale_y_continuous("Copy number",
                         expand = expansion(mult = 0.1)) +
    xlim(x_limits) +
      scale_color_manual(values = c("major_allele" = "#E69F00", "minor_allele" = "#0072B2", "CN_obs" = "#B0B0B0", "CN_call" = "#000000")) +
      theme(legend.title=element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
      labs(title = chr)
  } else {
    # Whole genome plot
    # Use adjusted seg_df for the whole genome plot
    seg_df <- seg_df_adj
    cytoband$chromStart <- cytoband$adjStart
    cytoband$chromEnd <- cytoband$adjEnd
    tumorBAF_df_adj$pos <- tumorBAF_df_adj$adjpos
    CNs_adj_plot$pos <- CNs_adj_plot$adjpos
    x_limits <- range(
      c(0, CNs_adj_plot$pos, tumorBAF_df_adj$pos, seg_df$endpos),
      na.rm = TRUE
    )

    # Use breaks and labels for the y-axis
    Segment_plot <- ggplot(fai) +
    geom_point(data = dplyr::mutate(CNs_adj_plot, track = "CN_obs"), aes(x = pos, y = CN_smooth, col = track), shape = 20) +
    geom_vline(aes(xintercept = start), col = "grey") +
    geom_linerange(data = seg_df, 
                   aes(xmin = adjstartpos, xmax = adjendpos, y = ascat_ploidy, col = allele), 
                   linewidth = 2.5, position = position_dodge(width = -0.1)) +
    geom_point(data = dplyr::mutate(CNs_adj_plot, track = "CN_call"), aes(x = pos, y = CN_seg, col = track), shape = 20) +
    geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5) +
    scale_y_continuous("Copy number",
                       breaks = breaks,
                       labels = labels,
                       expand = expansion(mult = 0.1),
                       limits = c(
                        -0.04*max(2, max(seg_df$ascat_ploidy, na.rm = TRUE)), 
                        max(2, max(seg_df$ascat_ploidy, na.rm = TRUE))
                        )) +
    xlim(x_limits) +
      scale_color_manual(values = c("major_allele" = "#E69F00", "minor_allele" = "#0072B2", "CN_obs" = "#B0B0B0", "CN_call" = "#000000")) +
    theme(legend.title=element_blank(), axis.title.x = element_blank())
  }
  
  # Add cytoband if provided
  if (!is.null(cytoband)) {
    Segment_plot <- Segment_plot +
      geom_rect(
        data = cytoband, 
        aes(
          xmin = chromStart, 
          xmax = chromEnd, 
          ymin = -0.04*max(2, max(seg_df$ascat_ploidy, na.rm = TRUE)), 
          ymax = -0.02*max(2, max(seg_df$ascat_ploidy, na.rm = TRUE)), 
          fill = color
        ), inherit.aes = FALSE) +
      scale_fill_identity()
  }
  
  BAF_plot <- ggplot(fai) +
    geom_point(data = tumorBAF_df_adj, aes(pos, BAF), shape = 20, size = 0.5, col = "#009E73") +
    scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9), limits = c(0.01,0.99)) +
    xlim(x_limits) +
    theme(axis.title.x = element_blank())
  
  if (is.null(chr)) {
    BAF_plot <- BAF_plot +
      geom_vline(aes(xintercept = start), col = "grey")
  }

  plot_grid(Segment_plot, BAF_plot, 
            ncol = 1, align = "v", rel_heights = c(3,1))
}

logR_to_CN <- function(logR_vector, purity = 1, round_to_integer = FALSE) {
  if (purity <= 0 || purity > 1) stop("purity must be in (0,1].")
  CN_obs <- 2 * 2^(logR_vector)
  CN_tumor <- (CN_obs - 2 * (1 - purity)) / purity
  CN_tumor[CN_tumor < 0] <- 0
  if (round_to_integer) CN_tumor <- round(CN_tumor)
  return(CN_tumor)
}

theme_set(theme_pubclean())
show_points = 5E5  # Number of points to show in the BAF/LogR plots

# Define command-line arguments
option_list <- list(
  make_option("--tumorname", type = "character", help = "Tumor sample name", metavar = "character"),
  make_option("--gender", type = "character", help = "Gender of the sample (e.g., XX or XY)", metavar = "character"),
  make_option("--genome-fai", type = "character", help = "Path to the genome FAI file", metavar = "character"),
  make_option("--segments", type = "character", help = "Path to the segments file", metavar = "character"),
  make_option("--Rdata-file", type = "character", help = "Path to the ascat run Rdata file", metavar = "character"),
  make_option("--cytoband", type = "character", help = "Path to the cytoband file", metavar = "character"),
  make_option("--output-plot", type = "character", help = "Path to output plot", metavar = "character"),
  make_option("--output-seg", type = "character", help = "Path to output segments file", metavar = "character"),
  make_option("--output-baf", type = "character", help = "Path to output BAF file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$`genome-fai`) || opt$`genome-fai` == "") {
  stop("Please provide a valid path to the genome FAI file using --genome-fai.")
}
if (is.null(opt$segments) || opt$segments == "") {
  stop("Please provide a valid path to the segments file using --segments.")
}
if (is.null(opt$`Rdata-file`) || opt$`Rdata-file` == "") {
  stop("Please provide a valid path to the ascat run Rdata file using --Rdata-file.")
}
if (is.null(opt$`output-plot`) || opt$`output-plot` == "") {
  stop("Please provide a valid path to the output plot using --output-plot.")
}

# Handle default values for optional arguments
if (is.null(opt$tumorname) || opt$tumorname == "") {
  opt$tumorname <- "Tumor"
}
if (is.null(opt$gender) || opt$gender == "") {
  opt$gender <- "XX"
}
if (is.null(opt$`output-seg`) || opt$`output-seg` == "") {
  opt$`output-seg` <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_copynumber.seg"))
}
if (is.null(opt$`output-baf`) || opt$`output-baf` == "") {
  opt$`output-baf` <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_baf.seg"))
}

## Load ascat.bc from the Rdata file
load(opt$`Rdata-file`)  # Load the Rdata file containing ascat.bc

# Map "male" and "female" to "XY" and "XX"
# Added for compatibility with the calc_sex() function
if (opt$gender == "male") {
  opt$gender <- "XY"
} else if (opt$gender == "female") {
  opt$gender <- "XX"
}

# The fai is used to get the chromosome lengths and plot them sequentially
fai <- read.table(opt$`genome-fai`, header = FALSE, sep = "\t", 
                  col.names = c("chr", "length", "start", "value1", "value2")) %>%
  mutate(chr = substr(chr,4,30))%>%
  mutate(middle = start + (length / 2))
if (opt$gender == "XY") {
  fai <- fai[1:24,]
} else {
  fai <- fai[1:23,]
  }

## Major Minor allele Segments
# Read the segments table
seg_df <- ascat.output$segments
# Merge with fai to adjust positions
seg_df <- merge(seg_df, fai, by = "chr") %>%
  mutate(
    adjstartpos = startpos + start,
    adjendpos = endpos + start,
    segment_size = endpos - startpos
  )%>%
  dplyr::rename(minor_allele = nMinor, major_allele = nMajor)%>%
  pivot_longer(cols = c(major_allele, minor_allele), names_to = "allele", values_to = "ascat_ploidy")

# If there are any segments with more than 3 times the median segment ploidy, cap them to 3 times the median
max_ploidy <- median(seg_df$ascat_ploidy, na.rm = TRUE) * 3
if (any(seg_df$ascat_ploidy > max_ploidy)) {
  seg_df_adj <- seg_df %>%
    mutate(ascat_ploidy = ifelse(
      ascat_ploidy > max_ploidy,
      # If the ploidy is greater than max_ploidy, subtract the rounded difference from ascat_ploidy
      ascat_ploidy - (round(ascat_ploidy - max_ploidy, 0)),
      ascat_ploidy
    ))
  # If the ploidy is capped, make it clear in the plot with the labels
  breaks <- seq(0, max_ploidy, by = 1)
  labels <- c(as.character(breaks[-length(breaks)]), paste0(">", max_ploidy))
} else {
  seg_df_adj <- seg_df
  breaks <- seq(0, max(seg_df_adj$ascat_ploidy, na.rm = TRUE), by = 1)
  labels <- as.character(breaks)
}


## Copynumbers
# Extract LogRs from the ascat.bc object and transform to CN
window <- 51
CNs <- data.frame(
  chr = ascat.bc$SNPpos$Chromosome,
  pos = ascat.bc$SNPpos$Position,
  logR_seg = ascat.bc$Tumor_LogR_segmented[,1],
  logR_raw = ascat.bc$Tumor_LogR[,1])%>%
  mutate(
    logR_smooth = runmed(logR_raw, k = window, endrule = "median"),
    CN_seg = logR_to_CN(logR_seg, purity = ascat.output$purity, round_to_integer = FALSE),
    CN_smooth = logR_to_CN(logR_smooth, purity = ascat.output$purity, round_to_integer = FALSE)
  )
# Merge with fai to adjust positions
CNs_adj <- merge(CNs, fai, by = "chr") %>%
  mutate(
    adjpos = pos + start
  )


## BAF
# Extract BAF from the ascat.bc object
setDT(ascat.bc$Tumor_BAF, keep.rownames = TRUE)[]%>%
  separate_wider_delim(rn, "_", names = c("chr","pos"))%>%
  mutate(pos = as.numeric(pos))%>%
  dplyr::select(1,2,"BAF"= 3)->tumorBAF_df
# Merge with fai to adjust positions
tumorBAF_df_adj <- merge(tumorBAF_df, fai, by = "chr") %>%
  mutate(
    adjpos = pos + start
  )


## Cytoband
# Read the cytoband file if provided
if (!is.null(opt$cytoband) && file.exists(opt$cytoband)) {
  cytoband <- data.frame(read.table(opt$cytoband, header = TRUE)) %>%
    mutate(chrom = substr(chrom, 4, 30)) %>%  # Remove the "chr" prefix from chromosome names
    left_join(fai, by = c("chrom" = "chr")) %>%  # Join with FAI data
    mutate(
      adjStart = start + chromStart,       # Adjust start positions
      adjEnd = start + chromEnd,           # Adjust end positions
      color = case_when(                   # Assign colors based on band type
        gieStain == "acen" ~ "#D55E00",
        gieStain == "gneg" ~ "#E0E0E0",
        gieStain == "gpos25" ~ "#C0C0C0",
        gieStain == "gpos50" ~ "#808080",
        gieStain == "gpos75" ~ "#404040",
        gieStain == "gpos100" ~ "#000000",
        gieStain == "stalk" ~ "#000000",
        gieStain == "gvar" ~ "#56B4E9",
        TRUE ~ "#FFFFFF"
      )
    )
} else {
  cytoband <- NULL
}

## Plot segments, BAF, and LogR
pdf(opt$`output-plot`, width = 18, height = 9)

# Reduce the number of points for the whole genome plot
CNs_adj_plot <- CNs_adj %>%
  slice(which(row_number() %% floor(n() / (show_points/5)) == 1))
tumorBAF_df_adj_plot <- tumorBAF_df_adj %>%
  slice(which(row_number() %% floor(n() / (show_points/5)) == 1))

print(plot_ascat_panels(fai, seg_df_adj, seg_df, tumorBAF_df_adj_plot, CNs_adj_plot, breaks, labels, cytoband = cytoband))

# Reduce the number of points for chromosome-specific plots
CNs_adj_plot <- CNs_adj %>%
  slice(which(row_number() %% floor(n() / show_points) == 1))
tumorBAF_df_adj_plot <- tumorBAF_df_adj %>%
  slice(which(row_number() %% floor(n() / show_points) == 1))

for (chr in unique(fai$chr)) {
  print(plot_ascat_panels(fai, seg_df_adj, seg_df, tumorBAF_df_adj_plot, CNs_adj_plot, breaks, labels, chr, cytoband = cytoband))
}

dev.off()

# Write the segmentation data to a file for IGV visualization
output_ratio_seg <- opt$`output-seg`

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0", con = output_ratio_seg)

# Process the segments file to output nMajor and nMinor as separate rows
CNs %>%
  dplyr::rename(Chromosome = chr, Start = pos) %>%
  mutate(Chromosome = paste0("chr", Chromosome), End = Start, Sample = opt$tumorname) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, CN_smooth) %>%
  write.table(file = output_ratio_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)


# Write the BAF data to an IGV-compatible file
output_baf_seg <- opt$`output-baf`

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=100:100:100 color=0,0,255 altColor=255,0,0", con = output_baf_seg)

# Process the BAF data for IGV visualization
tumorBAF_df %>%
  dplyr::rename(Chromosome = chr, Start = pos) %>%
  mutate(Chromosome = paste0("chr", Chromosome), End = Start, Sample = opt$tumorname) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, BAF) %>%
  write.table(file = output_baf_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)

