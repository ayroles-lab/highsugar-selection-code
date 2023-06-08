library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
withinChrR2_HS <- fread('/Genomics/ayroleslab2/simon/minjasSelection/results/2020-10-03_ldEtc/phasedHS_indFilter05_SNPfilter08_merged.geno.ld.gz')
withinChrR2_C <- fread('/Genomics/ayroleslab2/simon/minjasSelection/results/2020-10-09_phasingLDetc_C/phasedC_indFilter05_SNPfilter08_merged.ld.gz')
bins <- c(seq(from = 0, to = 500, by = 50), 3e7)

#HS
withinChrR2_HS$dist <- abs(withinChrR2_HS$POS1 - withinChrR2_HS$POS2)
withinChrR2_HS$dist_bin <- cut(withinChrR2_HS$dist, breaks = bins)

#C
withinChrR2_C$dist <- abs(withinChrR2_C$POS1 - withinChrR2_C$POS2)
withinChrR2_C$dist_bin <- cut(withinChrR2_C$dist, breaks = bins)

#Combine HS and C in the same fig instead
library(wesanderson)
keep <- which(withinChrR2_HS$dist_bin != '(500,3e+07]')
tmp.hs <- data.frame(withinChrR2_HS[keep, 5:7], treat = 'High Sugar')
keep <- which(withinChrR2_C$dist_bin != '(500,3e+07]')
tmp.c <- data.frame(withinChrR2_C[keep, 5:7], treat = 'Control')
dat <- rbind(tmp.hs, tmp.c)

#col.treat <- c(wes_palette('Darjeeling1', 1), 'deepskyblue2')
col.treat <- c("#B3B3FF", "#FFD9B3")

p = ggplot(data = dat, aes(x = dist_bin, y = R.2, fill = treat)) +
  geom_boxplot(width = .5, outlier.shape = NA) +
  ylab(expression('Linkage disequilibrium (' ~ R^2 ~ ')')) + xlab('Maximum Distance (bp)') +
  scale_x_discrete(labels= seq(50, 500, 50)) +
  theme_classic() + scale_fill_manual(name = "treat", values = col.treat) +
  theme(plot.title = element_text(hjust = 0., size = 30, face = 'bold'),
        legend.title=element_text(size=7, face = 'bold'),
        legend.text=element_text(size=8),
        legend.position = c(.778, .8), 
        axis.text=element_text(size=6),
        axis.text.x = element_text(angle = 0, size = 6),
        axis.title=element_text(size=8)) +
  guides(fill=guide_legend(title="Selection Regime"))
save_plot('figures/ldDecayG100_HS_C_samePanel.png', p, base_width = 7.2/3, base_height = 2 )



#Stuff for fig 1
#MAF distribution
load('data/190427_mergedFreqsANDnrSamples_unionSNPs.RData')
library(tidyverse)
library(ggridges)

#G1
g1.freqs <- freqs.merge_union_190426[,grep('G1_', colnames(freqs.merge_union_190426)), with = F]
#Don't know if this makes the most sense
g1.freqs_mean <- rowMeans(g1.freqs[, c(1,2,4,5,6)]) #Skipping N5 g1
g1.freqs_mean.maf <- g1.freqs_mean
g1.freqs_mean.maf[g1.freqs_mean.maf > .5] <- 1 - g1.freqs_mean.maf[g1.freqs_mean.maf > .5]

#G11
g11.freqs <- freqs.merge_union_190426[,grep('G11_', colnames(freqs.merge_union_190426)), with = F]
#Don't know if this makes the most sense
g11.freqs_mean <- rowMeans(g11.freqs)
g11.freqs_mean.maf <- g11.freqs_mean
g11.freqs_mean.maf[g11.freqs_mean.maf > .5] <- 1 - g11.freqs_mean.maf[g11.freqs_mean.maf > .5]

#G25
g25.freqs <- freqs.merge_union_190426[,grep('G25_', colnames(freqs.merge_union_190426)), with = F]
#Don't know if this makes the most sense
g25.freqs_mean <- rowMeans(g25.freqs)
g25.freqs_mean.maf <- g25.freqs_mean
g25.freqs_mean.maf[g25.freqs_mean.maf > .5] <- 1 - g25.freqs_mean.maf[g25.freqs_mean.maf > .5]

#G100
g100.freqs <- freqs.merge_union_190426[,grep('G100_N', colnames(freqs.merge_union_190426)), with = F]
#Don't know if this makes the most sense
g100.freqs_mean <- rowMeans(g100.freqs)
g100.freqs_mean.maf <- g100.freqs_mean
g100.freqs_mean.maf[g100.freqs_mean.maf > .5] <- 1 - g100.freqs_mean.maf[g100.freqs_mean.maf > .5]

tmp <- data.frame(g1 = g1.freqs_mean.maf, g11 = g11.freqs_mean.maf,
                  g25 = g25.freqs_mean.maf, g100 = g100.freqs_mean.maf)
dat <- pivot_longer(tmp, cols = c('g1', 'g11', 'g25', 'g100'))
dat$name <- factor(dat$name, levels = c('g1', 'g11', 'g25', 'g100'))

# png('../results/2020-11-25_msFigs/MAFdistr_v1.png', 1200, 1000)
p = ggplot(dat, aes(x = value, ill = name, color = name)) +
  geom_density(alpha = .5, key_glyph = "abline") + theme_classic() +
  scale_color_discrete(labels = c(1, 11, 25, 100)) +
  scale_x_continuous(limits = c(0, 0.5)) + scale_y_continuous(limits = c(0, 4.5)) + 
  theme(plot.title = element_text(hjust = 0., size = 25),
        legend.title=element_text(size=8, face = 'bold'),
        legend.text=element_text(size=6),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 0),
        legend.key.size = unit(0.5,"line"),
        legend.position = c(0.8, 0.75),
        axis.title=element_text(size=8)) + ylab('Density') + xlab('Minor Allele Frequency') +
        guides(color=guide_legend(title="Generation"))
save_plot('figures/MAFdistr_ridges.png', p, base_width = 7.2/3, base_height = 2 )

# dev.off()

#How about this
p = ggplot(dat, aes(x = value, y = name)) + #fill = ..x..
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  theme_minimal(base_size = 32) +
  # scale_fill_viridis(name = "value", option = "C") +
  xlab('Minor Allele Frequency') +
  ylab('Generation') +
  scale_y_discrete(labels = c(1, 11, 25, 100)) +
  theme(legend.position = 'none')
save_plot('figures/MAFdistr_ridges.png', p, base_width = 3, base_height = 2.5 )
