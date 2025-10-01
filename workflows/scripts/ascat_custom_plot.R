# Load required libraries
library(data.table)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(optparse)

# Function to plot ascat panels
plot_ascat_panels <- function(fai, seg_df, tumorBAF_df_adj, CNs_df_adj, max_scale = NULL, tumorname = NULL, chr = NULL, cytoband = NULL) {
  # Adjustments for scaling and titles if max_scale is not provided
  title_suffix <- ""
  if (is.null(max_scale)) {
    # The _plot columns are limited to the set y-scale, here we reset them to the actual values for automatic scaling
    seg_df$ascat_ploidy_plot <- seg_df$ascat_ploidy
    CNs_df_adj$CN_call_plot <- CNs_df_adj$CN_call
    CNs_df_adj$CN_smooth_plot <- CNs_df_adj$CN_smooth
    max_scale <- max(c(seg_df$ascat_ploidy, CNs_df_adj$CN_call, CNs_df_adj$CN_smooth), na.rm = TRUE)
    title_suffix <- " (unscaled)"
  }

  if (is.null(chr)) {
    # Whole genome plot
    # Use adjusted seg_df for the whole genome plot
    cytoband$chromStart <- cytoband$adjStart
    cytoband$chromEnd <- cytoband$adjEnd
    tumorBAF_df_adj$pos <- tumorBAF_df_adj$adjpos
    CNs_df_adj$pos <- CNs_df_adj$adjpos
    x_limits <- range(
      c(0, CNs_df_adj$pos, tumorBAF_df_adj$pos, seg_df$endpos),
      na.rm = TRUE
    )

    Segment_plot <- ggplot(fai) +
    geom_point(data = dplyr::mutate(CNs_df_adj, track = "CN_smooth"), aes(x = pos, y = CN_smooth_plot, col = track), shape = 20) +
    geom_vline(aes(xintercept = start), col = "grey") +
    geom_linerange(data = seg_df, aes(xmin = adjstartpos, xmax = adjendpos, y = ascat_ploidy_plot, col = allele), 
                   linewidth = 2.5, position = position_dodge(width = -0.15)) +
    geom_point(data = dplyr::mutate(CNs_df_adj, track = "CN_call"), aes(x = pos, y = CN_call_plot, col = track), shape = 20) +
    geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = c(-0.15, max_scale+0.1)) +
    scale_color_manual(values = c("major_allele" = "#E69F00", "minor_allele" = "#0072B2", "CN_smooth" = "#B0B0B0", "CN_call" = "#000000")) +
    theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(title = if (!is.null(tumorname)) paste0(tumorname, title_suffix) else "", y = "Copy number")
  } else {
    # Single chromosome plots
    # Filter data for the specified chromosome
    fai <- fai %>% filter(chr == !!chr)
    seg_df <- seg_df %>% filter(chr == !!chr)
    tumorBAF_df_adj <- tumorBAF_df_adj %>% filter(chr == !!chr)
    CNs_df_adj <- CNs_df_adj %>% filter(chr == !!chr)
    chr_title <- paste0("chr", chr)
    x_limits <- range(c(0, CNs_df_adj$pos, tumorBAF_df_adj$pos, seg_df$endpos), na.rm = TRUE)
    x_breaks <- seq(x_limits[1], x_limits[2], by = 1e7)
    x_labels <- paste0(x_breaks / 1e6, "Mb")
    if (!is.null(cytoband)) {
      cytoband <- cytoband %>% filter(chrom == !!chr)
    }

    Segment_plot <- ggplot(fai) +
    geom_point(data = dplyr::mutate(CNs_df_adj, track = "CN_smooth"), aes(x = pos, y = CN_smooth_plot, col = track), shape = 20) +
    geom_linerange(data = seg_df, 
                  aes(xmin = startpos, xmax = endpos, y = ascat_ploidy_plot, col = allele), 
                  linewidth = 2.5, position = position_dodge(width = -0.15)) +
    geom_point(data = dplyr::mutate(CNs_df_adj, track = "CN_call"), aes(x = pos, y = CN_call_plot, col = track), shape = 20) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(-0.15, max_scale+0.1)) +
    scale_color_manual(values = c("major_allele" = "#E69F00", "minor_allele" = "#0072B2", "CN_smooth" = "#B0B0B0", "CN_call" = "#000000")) +
    theme(legend.title=element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title = paste0(chr_title, title_suffix), y = "Copy number")
  }
  
  # Add cytoband if provided
  if (!is.null(cytoband)) {
    Segment_plot <- Segment_plot +
      geom_rect(
        data = cytoband, 
        aes(
          xmin = chromStart, 
          xmax = chromEnd, 
          ymin = -0.15, 
          ymax = -0.05, 
          fill = color
        ), inherit.aes = FALSE) +
      scale_fill_identity()
  }
  
  # BAF plot
  BAF_plot <- ggplot(fai) +
    geom_point(data = tumorBAF_df_adj, aes(pos, BAF), shape = 20, size = 0.5, col = "#009E73") +
    scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9), limits = c(0,1)) +
    xlim(x_limits) +
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
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

# Define command-line arguments
option_list <- list(
  make_option("--tumorname", type = "character", default = "Tumor", help = "Tumor sample name [default %default]", metavar = "character"),
  make_option("--gender", type = "character", default = "XY", help = "Gender of the sample (e.g., XX or XY) [default %default]", metavar = "character"),
  make_option("--genome-fai", type = "character", help = "Path to the genome FAI file", metavar = "character"),
  make_option("--Rdata-file", type = "character", help = "Path to the ascat run Rdata file", metavar = "character"),
  make_option("--cytoband", type = "character", help = "Path to the cytoband file", metavar = "character"),
  make_option("--output-plot", type = "character", help = "Path to output plot", metavar = "character"),
  make_option("--output-seg-smooth", type = "character", help = "Path to output segments file with smoothed copynumbers", metavar = "character"),
  make_option("--output-seg-call", type = "character", help = "Path to output segments file with ascat calls", metavar = "character"),
  make_option("--output-baf", type = "character", help = "Path to output BAF file", metavar = "character"),
  make_option("--whole-genome-points", type = "integer", default = 1E5, help = "Number of points to show in the whole genome BAF/LogR plots [default %default]", metavar = "integer"),
  make_option("--chromosome-points", type = "integer", default = 2E4, help = "Number of points to show in the chromosome-specific BAF/LogR plots [default %default]", metavar = "integer"),
  make_option("--smoothing-window", type = "integer", default = 51, help = "Window size for smoothing the LogR values [default %default]", metavar = "integer"),
  make_option("--default-y-scale", type = "integer", default = 6, help = "Default maximum y-scale for the copy number plot [default %default]", metavar = "integer")
  )

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$`genome-fai`) || opt$`genome-fai` == "") {
  stop("Please provide a valid path to the genome FAI file using --genome-fai.")
}
if (is.null(opt$`Rdata-file`) || opt$`Rdata-file` == "") {
  stop("Please provide a valid path to the ascat run Rdata file using --Rdata-file.")
}
if (is.null(opt$`output-plot`) || opt$`output-plot` == "") {
  stop("Please provide a valid path to the output plot using --output-plot.")
}

# Handle default values for optional arguments
if (is.null(opt$`output-seg-smooth`) || opt$`output-seg-smooth` == "") {
  opt$`output-seg-smooth` <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_CN_smooth.seg"))
}
if (is.null(opt$`output-seg-call`) || opt$`output-seg-call` == "") {
  opt$`output-seg-call` <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_CN_call.seg"))
}
if (is.null(opt$`output-baf`) || opt$`output-baf` == "") {
  opt$`output-baf` <- file.path(dirname(opt$`output-plot`), paste0(basename(opt$`output-plot`), "_ascat_BAF_IGV.seg"))
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
  pivot_longer(cols = c(major_allele, minor_allele), names_to = "allele", values_to = "ascat_ploidy")%>%
  mutate(ascat_ploidy_plot = pmin(ascat_ploidy, opt$`default-y-scale`))

## Copynumbers
# Extract LogRs from the ascat.bc object and transform to CN
CNs_df <- data.frame(
  chr = ascat.bc$SNPpos$Chromosome,
  pos = ascat.bc$SNPpos$Position,
  logR_seg = ascat.bc$Tumor_LogR_segmented[,1],
  logR_raw = ascat.bc$Tumor_LogR[,1])%>%
  merge(., fai, by = "chr") %>%
  mutate(
    logR_smooth = runmed(logR_raw, k = opt$`smoothing-window`, endrule = "median"),
    CN_call = logR_to_CN(logR_seg, purity = ascat.output$purity, round_to_integer = FALSE),
    CN_call_plot = pmin(CN_call, opt$`default-y-scale`),
    CN_smooth = logR_to_CN(logR_smooth, purity = ascat.output$purity, round_to_integer = FALSE),
    CN_smooth_plot = pmin(CN_smooth, opt$`default-y-scale`),
    adjpos = pos + start
  )


## BAF
# Extract BAF from the ascat.bc object
tumorBAF_df <- 
  setDT(ascat.bc$Tumor_BAF, keep.rownames = TRUE)[] %>%
  separate_wider_delim(rn, "_", names = c("chr","pos")) %>%
  merge(., fai, by = "chr") %>%
  mutate(pos = as.numeric(pos),
         adjpos = pos + start) %>%
  dplyr::select(1,2,"BAF"= 3, adjpos)
 

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
CNs_df_adj <- CNs_df %>%
  slice(which(row_number() %% ceiling(n() / opt$`whole-genome-points`) == 1))
tumorBAF_df_adj <- tumorBAF_df %>%
  slice(which(row_number() %% ceiling(n() / opt$`whole-genome-points`) == 1))

# Whole genome plot
print(plot_ascat_panels(fai, seg_df, tumorBAF_df_adj, CNs_df_adj, max_scale = opt$`default-y-scale`, tumorname = opt$tumorname, cytoband = cytoband))
if(
    any(seg_df$ascat_ploidy > opt$`default-y-scale`) |
    any(CNs_df$CN_call > opt$`default-y-scale`) |
    any(CNs_df$CN_smooth > opt$`default-y-scale`)
  ){
    print(plot_ascat_panels(fai, seg_df, tumorBAF_df_adj, CNs_df_adj, tumorname = opt$tumorname, cytoband = cytoband))
  }

# Reduce the number of points for chromosome-specific plots
CNs_df_adj <- CNs_df %>%
  group_by(chr) %>%
  slice(which(row_number() %% ceiling(n() / opt$`chromosome-points`) == 1)) %>%
  ungroup()
tumorBAF_df_adj <- tumorBAF_df %>%
  group_by(chr) %>%
  slice(which(row_number() %% ceiling(n() / opt$`chromosome-points`) == 1)) %>%
  ungroup()

# Plot each chromosome separately
for (chr in unique(fai$chr)) {
  print(plot_ascat_panels(fai, seg_df, tumorBAF_df_adj, CNs_df_adj, max_scale = opt$`default-y-scale`, chr = chr, cytoband = cytoband))
  if(
    any(seg_df$ascat_ploidy[seg_df$chr == chr] > opt$`default-y-scale`) |
    any(CNs_df$CN_call[CNs_df$chr == chr] > opt$`default-y-scale`) |
    any(CNs_df$CN_smooth[CNs_df$chr == chr] > opt$`default-y-scale`)
  ){
    print(plot_ascat_panels(fai, seg_df, tumorBAF_df_adj, CNs_df_adj, chr = chr, cytoband = cytoband))
  }
}

dev.off()

## Generate IGV-compatible segmentation and BAF files
# Write the smooth segmentation data to a file for IGV visualization
output_seg_smooth <- opt$`output-seg-smooth`

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0 viewLimits=0:10", con = output_seg_smooth)

# Process the segments file to output nMajor and nMinor as separate rows
CNs_df %>%
  dplyr::rename(Chromosome = chr, Start = pos) %>%
  mutate(Chromosome = paste0("chr", Chromosome), End = Start, Sample = opt$tumorname) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, CN_smooth) %>%
  write.table(file = output_seg_smooth, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)

# Write the call segmentation data to a file for IGV visualization
output_seg_call <- opt$`output-seg-call`

# Add IGV-compatible header
writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0 viewLimits=0:10", con = output_seg_call)

# Process the segments file to output nMajor and nMinor as separate rows
CNs_df %>%
  dplyr::rename(Chromosome = chr, Start = pos) %>%
  mutate(Chromosome = paste0("chr", Chromosome), End = Start, Sample = opt$tumorname) %>%  # Add "chr" prefix for IGV compatibility
  dplyr::select(Sample, Chromosome, Start, End, CN_call) %>%
  write.table(file = output_seg_call, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)


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

