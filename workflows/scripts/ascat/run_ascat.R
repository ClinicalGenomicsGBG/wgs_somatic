library(ASCAT)
library(data.table)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(plotly)
library(htmlwidgets)

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
genomefai = args[10]

gender = "XX"


fai <- read.table(genomefai, header = FALSE, sep = "\t", 
                  col.names = c("chr", "length", "start", "value1", "value2")) %>%
  mutate(chr = substr(chr,4,30))%>%
  mutate(middle = start + (length / 2))
if (gender == "XY") {
  fai <- fai[1:24,]
} else {
  fai <- fai[1:23,]
  }

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


mySeg <- read.table(file.path(output,paste0(tumorid,".segments_raw.txt")), header = T)
mySeg <- merge(mySeg, fai, by = "chr")
mySeg <- mySeg %>%
  mutate(
    startpos = startpos + start,
    endpos = endpos + start,
    segment_size = endpos - startpos,
    segment_frac = segment_size / sum(segment_size)
  )%>%
  rename(minor_allele = nBraw, major_allele = nAraw)%>%
  pivot_longer(cols = c(major_allele, minor_allele), names_to = "allele", values_to = "ascat_ploidy")

setDT(ascat.bc$Tumor_BAF, keep.rownames = T)[]%>%
  separate_wider_delim(rn, "_", names = c("chr","pos"))%>%
  mutate(pos = as.numeric(pos))%>%
  dplyr::select(1,2,"BAF"= 3)->tumorBAF_df
tumorBAF_df <- merge(tumorBAF_df, fai, by = "chr")
tumorBAF_df <- tumorBAF_df %>%
  mutate(
    pos = pos + start
  )

setDT(ascat.bc$Tumor_LogR, keep.rownames = T)[]%>%
  separate_wider_delim(rn, "_", names = c("chr","pos"))%>%
  mutate(pos = as.numeric(pos))%>%
  dplyr::select(1,2,"LogR"= 3)->tumorLogR_df
tumorLogR_df <- merge(tumorLogR_df, fai, by = "chr")
tumorLogR_df <- tumorLogR_df %>%
  mutate(
    pos = pos + start
  )

theme_set(theme_pubclean())

png(file.path(output,paste0(tumorid,".ascat_out.png")), width = 24, height = 24, units = "cm", res = 1200)

# Set the ylims according to segments longer than 1% of total segments (visible in the plot)
# Smaller fragments can be studied in interactive plot
lim_y <- c(0,0.5+max(mySeg$ascat_ploidy[mySeg$segment_frac>0.001])) 

Segment_plot <- ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_linerange(data = mySeg, 
                 aes(xmin = startpos, xmax = endpos, y = ascat_ploidy, col = allele), 
                 size = 2.5, position = position_dodge(width = -0.1))+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  ylab("Copy number")+
  coord_cartesian(ylim = lim_y)+
  theme(legend.title=element_blank(), axis.title.x = element_blank())

BAF_plot <-  ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_point(data = tumorBAF_df, aes(pos,BAF), shape = '.', col = "#00000010")+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9), 
                     limits = c(0.05,0.95)) +
  theme(axis.title.x = element_blank())

LogR_plot <-  ggplot(fai)+
  geom_vline(aes(xintercept = start), col = "grey")+
  geom_point(data = tumorLogR_df, aes(pos,LogR), shape = '.', col = "#00000010")+
  geom_text(aes(label = chr, x = middle, y = Inf), vjust = 1, size = 3.5)+
  scale_y_continuous(limits = c(-2,2))+
  theme(axis.title.x = element_blank())

plot_grid(Segment_plot, BAF_plot, LogR_plot, 
          ncol = 1, align = "v", rel_heights = c(3,2,2))

dev.off()

vline <- function(x = 0, color = "#00000033") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash = "dot")
  )
}

vlines <- list()
for (i in 1:nrow(fai)) {
  vlines <-  append(vlines, list(vline(fai$start[i])))
}

plot_SEG <- mySeg%>%
  mutate(gap = NA)%>%
  pivot_longer(cols = c(startpos, endpos, gap), names_to = "pos_name", values_to = "pos")%>%
  mutate(ascat_ploidy = ascat_ploidy + ifelse(allele == "major_allele", 0.01, -0.01))

SEG <- plot_ly(plot_SEG,
               x = ~pos,
               y = ~ascat_ploidy,
               color = ~allele,
               colors = c("#FC4E07","#2E9FDF"),
               type = "scatter", mode = "lines")%>%
  layout(shapes = vlines,
         annotations = list(
           x = fai$middle,
           y = 1,
           text = fai$chr,
           xref = "x",
           yref = "paper",
           showarrow = F
         ),
        yaxis = list(title = list(text = "Copy number"),
                     zerolinecolor = "#FFFFFF"),
        xaxis = list(zerolinecolor = "#FFFFFF"),
        legend = list(orientation = "h",
                      y = 100,
                      x = 0)
         )


show_points = 2E4

tumorBAF_df_red <- tumorBAF_df%>%
  dplyr::filter(BAF < 0.99 & BAF > 0.01)%>%
  slice(which(row_number() %% floor(n()/show_points) == 1))

BAF <- plot_ly(tumorBAF_df_red,
               x = ~pos,
               y = ~BAF,
               type = "scatter",
               mode = "markers",
               marker = list(
                 color = "#000000CC",
                 size = 2
               ),
               showlegend = F)%>%
  layout(shapes = vlines,
         annotations = list(
           x = fai$middle,
           y = 1,
           text = fai$chr,
           xref = "x",
           yref = "paper",
           showarrow = F
         ),
         yaxis = list(title = list(text = "BAF"),
                      zerolinecolor = "#FFFFFF"),
         xaxis = list(zerolinecolor = "#FFFFFF"))

tumorLogR_df_red <- tumorLogR_df%>%
  slice(which(row_number() %% floor(n()/show_points) == 1))

LogR <- plot_ly(tumorLogR_df_red,
               x = ~pos,
               y = ~LogR,
               type = "scatter",
               mode = "markers",
               marker = list(
                 color = "#000000CC",
                 size = 2
               ),
               showlegend = F)%>%
  layout(shapes = vlines,
         annotations = list(
           x = fai$middle,
           y = 1,
           text = fai$chr,
           xref = "x",
           yref = "paper",
           showarrow = F
         ),
         yaxis = list(title = list(text = "LogR"),
                      zerolinecolor = "#FFFFFF"),
         xaxis = list(zerolinecolor = "#FFFFFF"))



p <- subplot(SEG,BAF,LogR, nrows = 3, shareX = T, titleY = T, titleX = F)

saveWidget(p, file.path(output,paste0(tumorid,".ascat_interactive.html")), selfcontained = T)
