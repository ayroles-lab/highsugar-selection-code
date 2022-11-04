########## 200727 #########
#Redo the "what happens to the PCs?" analysis excluding the C only signals
load('../data/190916_maf05SNPs_alleleFreqs.RData')
load('../data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')
#Remove C only selection signals
nonSignCut <- 1e-4
p <- binGlmScan_contTime_param1234_200615.gr$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615.gr$emTrend.HS_p > nonSignCut] <- 1

library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(data.table)
library(stringr)
library(ggthemes)
library(cowplot)
library(extrafont)
# p_int <- 10^-c(10:15)
# p_int <- 10^-seq(10.5, 12, by = .2)
{
  #p_int <- sort(c(10^(seq(-9, -14, length.out = 4)), 8*10^-12))
  p_int = c(1*10^-14, 5*10^-13, 8*10^-12, 1*10^-10, 1*10^-9)
  titleCols <- brewer.pal(n = length(p_int), name = 'Spectral')
  col.treat <- c(wes_palette('Darjeeling1', 1), 'deepskyblue2')
  dat = vector("list", length(p_int))
  for (i in 1:length(p_int)) {
    print(i)
    #Get SNPs
    keep <- p > p_int[i]
    keep.snps <- binGlmScan_contTime_param1234_200615.gr@ranges@NAMES[keep]

    #PCA
    allAbove.freqs_keep <- allAbove.freqs[, colnames(allAbove.freqs) %in% keep.snps]
    allAbove.freqs_keep.pc <- prcomp(allAbove.freqs_keep)

    #plot
    dat[[i]] <- data.frame(allAbove.freqs_keep.pc$x)
    dat[[i]]$PC1 <- -dat[[i]]$PC1 #Just flipping PC1 for visualization
    if(i == 5) dat[[i]]$PC2 <- -dat[[i]]$PC2
    samples <- rownames(allAbove.freqs_keep.pc$x)
    t <- as.numeric(sub('G(.*)_.*', '\\1', samples))
    pop <- sub('.*_(.*)', '\\1', samples)
    dat[[i]]$t <- t
    treat <- rep('Control', nrow(dat[[i]]))
    treat[grep(pattern = 'N4|N5|N6', x = pop)] <- 'HS'
    dat[[i]]$trat <- treat
    dat[[i]]$pop <- pop
  }
}
{
  g = vector("list", length(p_int))
  for (i in 1:length(p_int)) {
    g[[i]] <- ggplot(dat[[i]], aes(x = PC1, y = PC2, group = pop, colour = treat)) +
      geom_line(size = 1) + geom_point(size = 3, aes(shape = factor(t))) +
      scale_colour_manual(name = "treat", values = col.treat) +
      theme(plot.title = element_text(hjust = 0.5, size = 25, color = titleCols[i], face="bold.italic"),
            legend.title=element_text(size=24),
            legend.text=element_text(size=20),
            axis.text=element_text(size=25),
            axis.title=element_text(size=28,face="bold"),
            strip.text.y = element_text(size = 20)) +
      guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) +
      ggtitle(bquote(SNPs ~ with ~ p[int] ~ '>' ~ .(p_int[i]))) + theme_tufte()
      if(i < length(p_int)) g[[i]] = g[[i]] + theme(legend.position = "none")
  }
}

######## New Manhatan plot 221103 ############
# Fig 2D
source("utils.R")
load('../data/201001_alleFreqPC_allSNPs.RData')
summary(allSNPs.pc)

library(wesanderson)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(grid)

load(
    "../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData"
)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in%
    chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-04
signCut <- 8e-12
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
binGlmScan_contTime_param1234_200615$p = p

chrs <- c("2L", "2R", "3L", "3R", "4", "X")
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in%
    chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-04
signCut <- 8e-12
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
binGlmScan_contTime_param1234_200615$p = p

nonSignCut <- 1e-04
signCut <- 8e-12

gwas_data_load = binGlmScan_contTime_param1234_200615
gwas_data_load$chr = factor(binGlmScan_contTime_param1234_200615$CHROM, levels = chrs)
gwas_data_load$bp = binGlmScan_contTime_param1234_200615$POS
sig_data <- gwas_data_load %>% 
  subset(p < signCut)
notsig_data <- gwas_data_load %>% 
  subset(p >= signCut) %>%
  group_by(chr) %>% 
  sample_frac(0.2)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(p), 
                                  color = chr)) +
  geom_point(alpha = 0.75, size = 0.25) +
  geom_hline(data = data.frame(p_int = p_int), aes(yintercept = -log10(p_int)), color = "black") + 
  geom_hline(yintercept = -log10(8*10^-12), color = "tomato3") + 
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks = seq(0, 60, 10)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromossome", 
       y = expression("-log10(" ~ p[int] ~ ")")) + 
  theme_tufte() +
  theme( 
    legend.position = "none",
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #axis.title.y = element_markdown(),
    axis.text.x = element_text(size = 8, vjust = 0.5),
    axis.line = element_line(color = 'black')
  )

  g[[6]] = manhplot
  panel_plot = plot_grid(plotlist = g, ncol = 2, nrow = 3, labels = 1:5)
  save_plot("../figures/sequential_PCA.png", panel_plot, base_height = 3.2, base_asp = 1.2, ncol = 2, nrow = 3)
