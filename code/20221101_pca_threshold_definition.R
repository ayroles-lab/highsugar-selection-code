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
  p_int = c(1*10^-14, 5*10^-13, 8*10^-12, 2.15*10^-11, 1*10^-9)
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
  panel_plot = plot_grid(plotlist = g, ncol = 2, nrow = 3, labels = 1:5)
  save_plot("../figures/sequential_PCA.png", panel_plot, base_height = 3.2, base_asp = 1.2, ncol = 2, nrow = 3)
}
#New fig2
library(grid)
load('../data/190916_maf05SNPs_freqPCs.RData')
load('../data/190427_mergedFreqsANDnrSamples_unionSNPs.RData')
samples <- rownames(allAbove.pc$x)
t <- as.numeric(sub('G(.*)_.*', '\\1', samples))
pop <- sub('.*_(.*)', '\\1', samples)
dat <- data.frame(allAbove.pc$x)
dat$pop <- pop
dat$t <- t
treat <- rep('Control', nrow(dat))
treat[grep(pattern = 'N4|N5|N6', x = pop)] <- 'High Sugar'
dat$treat <- factor(treat, levels = c('High Sugar', 'Control')) #To set the colors
dat$PC1 <- -dat$PC1 #Just flipping PC1 for visualization

#Trying with the var expl as inset
col.treat <- c(wes_palette('Darjeeling1', 1), 'deepskyblue2')
vp <- viewport(width = 0.3, height = 0.33, x = 0.205, y = 0.83)

pc <- ggplot(dat, aes(x = PC1, y = PC2, group = pop, colour = treat)) +
  geom_line(size = 2) + geom_point(size = 8, aes(shape = factor(t))) +
  scale_colour_manual(name = "treat", values = col.treat) +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=35,face="bold"),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20)) + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation"))

eig <- allAbove.pc$sdev^2
varExpl <- eig/sum(eig)
dat <- data.frame(pc = 1:length(varExpl), varExpl = varExpl)
pc.var <- ggplot(dat, aes(x = pc, y = varExpl, group = 1)) + geom_bar(stat = "identity") +
  xlab('PC') + ylab('Variance Explained') +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20))

png('../figures/PCs12_inset_v2.png', width = 1600, height = 1000)
print(pc)
print(pc.var, vp = vp)
dev.off()

#Manhattan
source('utils.R')
load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
chrs <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in% chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-4
signCut <- 8e-12
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
col <- rep('black', nrow(binGlmScan_contTime_param1234_200615))
col[binGlmScan_contTime_param1234_200615$emTrend.C_p > nonSignCut] <- wes_palette('Darjeeling1', 1)  #C_p > nonSignCut => HS only signals
col[binGlmScan_contTime_param1234_200615$`generation:treatment_p` > signCut] <- 'grey'
png('../figures/manhattan_pInt_colByTreatSlopeP_HSonly_v3.png', 1600, 1200)
par(mar = c(5,5,4,2))
plot.manhattan3(pos = binGlmScan_contTime_param1234_200615$POS,
                chr = binGlmScan_contTime_param1234_200615$CHROM,
                y = p, col = col, ylim = ylim)
legend('topleft', c('High Sugar', 'Both Treatments'), col = c(pal[1], 'black'), pch = 19, cex = 3,
       title = expression(bold('Selection Response')))
abline(h = -log10(signCut), lty = 2, lwd = 2)
mtext(expression('-log10(' ~ p[int] ~ ')'), side = 2, cex = 2, line = 2.7)
dev.off()

#Examples
binGlmScan_contTime_param1234_200615[order(p), ] |> head()
f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3R' & freqs.merge_union_190426$POS == 24624973, ]
# f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3R' & freqs.merge_union_190426$POS == 24632146, ]
png(paste0('../figures/freqEx_', f$CHROM, f$POS, '_v3.png'), 1200, 1000)
plot.freqVStime_v2(f[,9:32], lines = T, cex.legend = 3, pal = col.treat, cex.points = 5, legend = F, ylim = c(0,1))
dev.off()

#only HS
binGlmScan_contTime_param1234_200615[p < 1e-30 & binGlmScan_contTime_param1234_200615$CHROM == '2R', ]
f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2R' & freqs.merge_union_190426$POS == 12422213, ]
png(paste0('../figures/freqEx_', f$CHROM, f$POS, '.png'), 1200, 1000)
plot.freqVStime_v2(f[,9:32], lines = T, cex.legend = 3, pal = col.treat, cex.points = 5, legend = F)
dev.off()

#only HS ex 2
binGlmScan_contTime_param1234_200615[p < 1e-21
                                     & binGlmScan_contTime_param1234_200615$CHROM == '2L'
                                     & binGlmScan_contTime_param1234_200615$POS > .85e7
                                     & binGlmScan_contTime_param1234_200615$POS < .95e7, ]
# f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2L' & freqs.merge_union_190426$POS == 8712764, ]
f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2L' & freqs.merge_union_190426$POS == 8712101, ]
png(paste0('../figures/freqEx_', f$CHROM, f$POS, '.png'), 1200, 1000)
plot.freqVStime_v2(f[,9:32], lines = T, cex.legend = 3, pal = col.treat, cex.points = 5, legend = F)
dev.off()

#only HS ex 3
binGlmScan_contTime_param1234_200615[p < 1e-15
                                     & binGlmScan_contTime_param1234_200615$CHROM == '3R'
                                     & binGlmScan_contTime_param1234_200615$POS > 2.6e7
                                     & binGlmScan_contTime_param1234_200615$POS < 2.8e7, ]
# f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2L' & freqs.merge_union_190426$POS == 8712764, ]
f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3R' & freqs.merge_union_190426$POS == 27119027, ]
png(paste0('../figures/freqEx_', f$CHROM, f$POS, '_v3.png'), 1200, 1000)
plot.freqVStime_v2(f[,9:32], lines = T, cex.legend = 3, pal = col.treat, cex.points = 5, legend = F, ylim = c(0,1))
dev.off()
