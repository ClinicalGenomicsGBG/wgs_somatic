# Load required libraries
library(ggpubr)
library(dplyr)
library(cowplot)

# Function to calculate the most common increment value in a vector
calculate_increment <- function(vector, num_rows = 5) {
  # Ensure there are enough rows to compare
  if (length(vector) < num_rows) {
    stop("Not enough values to compare.")
  }
  
  # Calculate the differences between the Start values of the first few rows
  increments <- diff(vector[1:num_rows])
  
  # Return the most common increment value
  return(as.numeric(names(sort(table(increments), decreasing = TRUE)[1])))
}

# Constants
show_points = 5E5  # Number of points to show in the BAF plot
max_cutoff = 3     # If the ratio is greater than this value * baseline, it will be capped to this value

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
tumor_id <- args[1]               # Tumor sample ID
ratio_file <- args[2]             # File containing ratio data
BAF_file <- args[3]               # File containing BAF data
fai_file <- args[4]               # FAI file (genome index) used to adjust chromosome positions
cytoband <- args[5]               # Cytoband file (chromosome banding information)
output_ratio_plot <- args[6]      # Output file for the ratio plot
output_ratio_seg <- args[7]       # Output file for the ratio segmentation

# Set the default theme for plots
theme_set(theme_pubclean())

# Check if the FAI file exists
if (!file.exists(fai_file)) {
  stop(paste("FAI file not found:", fai_file))
}

# Check if the ratio file exists and is not empty
if (!file.exists(ratio_file) || file.info(ratio_file)$size == 0) {
  cat("WARNING: Ratio file is missing or empty. Creating placeholder outputs.\n")
  # Create empty placeholder outputs
  png(output_ratio_plot)
  plot.new()
  text(0.5, 0.5, "No data available", cex = 2)
  dev.off()
  
  write.table(data.frame(), output_ratio_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  quit(status = 0)
} else {
  # Read the ratio data into a data frame
  ratio <- data.frame(read.table(ratio_file, header=TRUE))
}

# Read the BAF data into a data frame, skipping header lines starting with '#'
if (!file.exists(BAF_file) || file.info(BAF_file)$size == 0) {
  cat("WARNING: BAF file is missing or empty. Proceeding without BAF data.\n")
  BAF <- data.frame()  # Create an empty data frame
} else {
  BAF <- read.table(
    BAF_file, header = FALSE, comment.char = "#", sep = "\t",
    col.names = c("Chromosome", "Start", "End", "Features", "BAF"))%>%
    mutate(Chromosome = substr(Chromosome, 4, 30)) # Remove the "chr" prefix from chromosome names
}

# Read the FAI file and process it
fai <- read.table(fai_file, header = FALSE, sep = "\t", 
                  col.names = c("chr", "length", "start", "value1", "value2")) %>%
  mutate(chr = substr(chr, 4, 30)) %>%  # Remove the "chr" prefix from chromosome names
  mutate(middle = start + (length / 2)) # Calculate the middle position of each chromosome
fai <- fai[1:24,]  # Keep only the first 24 rows (autosomes and sex chromosomes)

# Read the cytoband file and process it
cytoband <- data.frame(read.table(cytoband, header = TRUE)) %>%
  mutate(chrom = substr(chrom, 4, 30)) %>%  # Remove the "chr" prefix from chromosome names
  left_join(fai, by = c("chrom" = "chr")) %>%  # Join with FAI data
  mutate(adjStart = start + chromStart,       # Adjust start positions
         adjEnd = start + chromEnd,           # Adjust end positions
         color = case_when(                   # Assign colors based on band type
            gieStain == "acen" ~ "#C00000",
            gieStain == "gneg" ~ "#E0E0E0",
            gieStain == "gpos25" ~ "#C0C0C0",
            gieStain == "gpos50" ~ "#808080",
            gieStain == "gpos75" ~ "#404040",
            gieStain == "gpos100" ~ "#000000",
            gieStain == "stalk" ~ "#000000",
            gieStain == "gvar" ~ "#799FC7",
            TRUE ~ "#FFFFFF"
         ))

# Calculate the ploidy based on the ratio data
valid_ratio <- ratio %>% filter(MedianRatio > 0.8 & MedianRatio < 1.2)

if (nrow(valid_ratio) > 0) {
  ploidy <- median(valid_ratio$CopyNumber, na.rm = TRUE)
} else {
  ploidy <- 2  # Default ploidy value
}

cat(c("INFO: Selected ploidy:", ploidy, "\n"))

# Calculate the most common increment value in the ratio data
increment <- calculate_increment(ratio$Start)

# Adjust the ratio data for plotting
ratio_adjusted <- ratio %>%
  left_join(fai, by = c("Chromosome" = "chr")) %>%
  mutate(adjStart = Start + start,  # Adjust start positions
         End = Start + increment,  # Calculate end positions
         corr_ratio = case_when(   # Correct the ratio values
            Ratio * ploidy > max_cutoff * ploidy ~ Ratio * ploidy - (round(Ratio * ploidy) - max_cutoff * ploidy),
            TRUE ~ Ratio * ploidy),
         CopyNumber_adj = case_when(  # Adjust copy number values
            CopyNumber > max_cutoff * ploidy ~ max_cutoff * ploidy,
            TRUE ~ CopyNumber),
         color = case_when(  # Assign colors based on ratio values
            Ratio >= 0.9 & Ratio <= 1.1 ~ "around 1",
            Ratio < 0.9 & Ratio >= 0 ~ "below 0.9",
            Ratio > 1.1 ~ "above 1.1",
            TRUE ~ "other"
         )) %>%
  filter(Ratio > 0)  # Filter out invalid rows

# Calculate breaks and labels for the y-axis
max_corr_ratio <- max(ratio_adjusted$corr_ratio, na.rm = TRUE)
breaks <- seq(0, max_corr_ratio, by = 1)

# If the maximum ratio after correction is the same as before correction, 
# it means the ratio is not adjusted and we can use the original labels
if (max_corr_ratio == max(ratio_adjusted$Ratio * ploidy, na.rm = TRUE)) {
  labels <- as.character(breaks)
} else {
  # Otherwise, we need to adjust the labels to reflect the correction
  labels <- c(as.character(breaks[-length(breaks)]), paste0(">", max_corr_ratio - 1))
}

# Create the ratio plot
Ratio_plot <- ggplot(fai) +                                              
  geom_rect(data = cytoband, aes(xmin = adjStart, xmax = adjEnd, ymin = -(max_cutoff * ploidy * 0.035), ymax = -(max_cutoff * ploidy * 0.01), fill = color)) +
  geom_point(data = ratio_adjusted, aes(adjStart, corr_ratio, color = color), shape = '.') +
  geom_point(data = ratio_adjusted, aes(adjStart, CopyNumber_adj), shape = '.', col = "#00000070") +
  geom_vline(aes(xintercept = start), col = "grey") +
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3) +
  scale_y_continuous("Copy number",
                     breaks = breaks,
                     labels = labels,
                     expand = expansion(mult = 0.1)) +
  scale_fill_identity() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Process the BAF data for plotting
BAF_adjusted <- BAF %>%
  left_join(fai, by = c("Chromosome" = "chr")) %>%
  mutate(adjPosition = Start + start) %>%
  filter(BAF > 0) %>%  # Keep rows with valid BAF values
  select(Chromosome, Start, adjPosition, BAF, Features)

# Reduce the number of points for the BAF plot
BAF_adjusted_red <- BAF_adjusted %>%
  slice(which(row_number() %% floor(n() / show_points) == 1))

# Create the BAF plot
BAF_plot <- ggplot(fai) +                                              
  geom_vline(aes(xintercept = start), col = "grey") +
  geom_point(data = BAF_adjusted_red, aes(adjPosition, BAF), shape = '.', col = "#00000020") +
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3) +
  scale_y_continuous("BAF", limits = c(0, 1), expand = expansion(mult = 0.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save the combined ratio and BAF plots as a PNG file
png(output_ratio_plot, width = 24, height = 12, units = "cm", pointsize = 20, res = 1200)
plot_grid(Ratio_plot, BAF_plot, ncol = 1, align = "v", rel_heights = c(5, 2))
dev.off()

# Write the ratio segmentation data to a file
writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0", con = output_ratio_seg)
ratio_adjusted %>%
  mutate(Sample = tumor_id,
         Chromosome = paste0("chr", Chromosome),
         Copynumber = Ratio * ploidy) %>%
  select(Sample, Chromosome, Start, End, Copynumber) %>%
  write.table(output_ratio_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
