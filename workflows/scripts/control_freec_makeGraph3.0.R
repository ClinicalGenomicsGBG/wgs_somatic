
library(ggpubr)
library(dplyr)
library(cowplot)

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

show_points = 5E5
max_cutoff = 3  # Maximum ratio to plot

args <- commandArgs(trailingOnly = TRUE)
tumor_id <- args[1]
ratio_file <- args[2]
BAF_file <- args[3]
fai_file <- args[4]
cytoband <- args[5]
output_ratio_plot <- args[6]
output_ratio_seg <- args[7]
output_BAF_igv <- args[8]

# Check if the fai file exists
if (!file.exists(fai_file)) {
  stop(paste("FAI file not found:", fai_file))
}

theme_set(theme_pubclean())

ratio <- data.frame(read.table(ratio_file, header=TRUE))

BAF <- data.frame(read.table(BAF_file, header=TRUE))

fai <- read.table(fai_file, header = FALSE, sep = "\t", 
                  col.names = c("chr", "length", "start", "value1", "value2")) %>%
  mutate(chr = substr(chr,4,30))%>%
  mutate(middle = start + (length / 2))
fai <- fai[1:24,]

cytoband <- data.frame(read.table(cytoband, header = TRUE)) %>%
  mutate(chrom = substr(chrom,4,30)) %>%
  left_join(fai, by = c("chrom" = "chr")) %>%
  mutate(adjStart = start + chromStart,
         adjEnd = start + chromEnd,
         color = case_when(
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

ploidy = median(ratio$CopyNumber[which(ratio$MedianRatio>0.8 & ratio$MedianRatio<1.2)], na.rm = T)
cat (c("INFO: Selected ploidy:", ploidy, "\n"))

increment <- calculate_increment(ratio$Start)

ratio_adjusted <- ratio %>%
  left_join(fai, by = c("Chromosome" = "chr")) %>%
  mutate(adjStart = Start + start,
         End = Start + increment,
         corr_ratio = case_when(
            Ratio * ploidy > max_cutoff * ploidy ~ Ratio * ploidy - (round(Ratio * ploidy) - max_cutoff * ploidy),
            TRUE ~ Ratio * ploidy),
         CopyNumber_adj = case_when(
            CopyNumber > max_cutoff * ploidy ~ max_cutoff * ploidy,
            TRUE ~ CopyNumber),
         color = case_when(
            Ratio >= 0.9 & Ratio <= 1.1 ~ "around 1",
            Ratio < 0.9 & Ratio >= 0 ~ "below 0.9",
            Ratio > 1.1 ~ "above 1.1",
            TRUE ~ "other"
         ))%>%
  filter(Ratio > 0, BAF > 0)

# Calculate breaks and labels
max_corr_ratio <- max(ratio_adjusted$corr_ratio, na.rm = TRUE)
breaks <- seq(0, max_corr_ratio, by = 1)
labels <- c(as.character(breaks[-length(breaks)]), paste0(">", max(breaks) - 1))

Ratio_plot <- ggplot(fai) +                                              
  geom_rect(data = cytoband, aes(xmin = adjStart, xmax = adjEnd, ymin = -(max_cutoff*ploidy*0.035), ymax = -(max_cutoff*ploidy*0.01), fill = color)) +
  geom_point(data = ratio_adjusted, aes(adjStart, corr_ratio, color = color), shape = '.') +
  geom_point(data = ratio_adjusted, aes(adjStart, CopyNumber_adj), shape = '.', col = "#00000050") +
  geom_vline(aes(xintercept = start), col = "grey") +
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3) +
  scale_y_continuous("Copy number",
                     breaks = breaks,
                     labels = labels,
                     expand = expansion(mult = 0.1)) +
  scale_fill_identity() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

BAF_adjusted <- BAF %>%
  left_join(fai, by = c("Chromosome" = "chr")) %>%
  mutate(adjPosition = Position + start) %>%
  filter(BAF > 0, uncertainty < 1) %>%
  select(Chromosome, Position, adjPosition, BAF, FittedA, FittedB, A, B, uncertainty)

BAF_adjusted_red <- BAF_adjusted%>%
  slice(which(row_number() %% floor(n()/show_points) == 1))
BAF_plot <-  ggplot(fai)+                                              
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_point(data = BAF_adjusted_red, aes(adjPosition,BAF),shape='.', col = "#00000010")+
  geom_point(data = BAF_adjusted_red[BAF_adjusted_red$uncertainty > 0 & BAF_adjusted_red$uncertainty < 1,], aes(adjPosition,FittedA),shape=15, size = 0.2, col = "#0090FF")+
  geom_point(data = BAF_adjusted_red[BAF_adjusted_red$uncertainty > 0 & BAF_adjusted_red$uncertainty < 1,], aes(adjPosition,FittedB),shape=15, size = 0.2, col = "#FF5050")+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3) +
  scale_y_continuous("BAF", limits = c(0, 1), expand = expansion(mult = 0.1))+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())


png(output_ratio_plot, width = 24, height = 12, units = "cm", pointsize = 20, res = 1200)
plot_grid(Ratio_plot, BAF_plot, ncol = 1, align = "v", rel_heights = c(5, 2))
dev.off()


writeLines("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0", con = output_ratio_seg)
ratio_adjusted %>%
  mutate(Sample = tumor_id,
         Chromosome = paste0("chr", Chromosome)) %>%
  select(Sample, Chromosome, Start, End, corr_ratio) %>%
  write.table(output_ratio_seg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)


writeLines(paste0(
  '#type=GENE_EXPRESSION\n#track graphtype=points name="',
  tumor_id,
  '" color=0,0,255 altColor=255,0,0 maxHeightPixels=160:160:160 viewLimits=0:1',
  "\n#Chromosome\tStart\tEnd\tFeatures\tvalues"
), con = output_BAF_igv)
BAF_adjusted %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Position,
         End = Position,
         Features = "") %>%
  select(Chromosome, Start, End, Features, BAF) %>%
  write.table(output_BAF_igv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
