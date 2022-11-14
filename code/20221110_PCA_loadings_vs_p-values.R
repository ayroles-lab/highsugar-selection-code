source("utils.R")
# 2020-11-25
# Fig 1a
load('../data/201001_alleFreqPC_allSNPs.RData')
summary(allSNPs.pc)

library(ggplot2)
library(grid)
library(wesanderson)

{
    samples <- rownames(allSNPs.pc$x)
    t <- as.numeric(sub('G(.*)_.*', '\\1', samples))
    pop <- sub('.*_(.*)', '\\1', samples)
    dat <- data.frame(allSNPs.pc$x)
    dat$pop <- pop
    dat$t <- t
    treat <- rep('Control', nrow(dat))
    treat[grep(pattern = 'N4|N5|N6', x = pop)] <- 'High Sugar'
    dat$treat <- factor(treat, levels = c('High Sugar', 'Control')) #To set the colors
    dat$PC1 <- -dat$PC1 #Just flipping PC1 for visualization
}


### PCs SI figure
#Supplementary loadings manhattan

load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
chrs <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
allSNPs.pc_load <- data.frame(CHROM = binGlmScan_contTime_param1234_200615$CHROM,
                              POS = binGlmScan_contTime_param1234_200615$POS,
                              allSNPs.pc$rotation)
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in% chrs, ]
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM != 'Y', ]
allSNPs.pc_load <- allSNPs.pc_load[allSNPs.pc_load$CHROM %in% chrs, ]
allSNPs.pc_load <- allSNPs.pc_load[allSNPs.pc_load$CHROM != 'Y', ]

cor.test(-allSNPs.pc_load[,"PC1"], binGlmScan_contTime_param1234_200615$`generation_p`)
cor.test(allSNPs.pc_load[,"PC2"], binGlmScan_contTime_param1234_200615$`generation:treatment_p`)
cor.test(allSNPs.pc_load[,"PC2"], allSNPs.pc_load[,"PC1"])

load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
allSNPs.pc_load <- data.frame(CHROM = binGlmScan_contTime_param1234_200615$CHROM,
                              POS = binGlmScan_contTime_param1234_200615$POS,
                              allSNPs.pc$rotation)
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in% chrs, ]
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM != 'Y', ]
allSNPs.pc_load <- allSNPs.pc_load[allSNPs.pc_load$CHROM %in% chrs, ]
allSNPs.pc_load <- allSNPs.pc_load[allSNPs.pc_load$CHROM != 'Y', ]
col <- brewer.pal(6, 'GnBu')[c(4,6)]
png('../figures/manhattan_allSNPsPC12_loadings.png', 1600, 1200)
par(mar = c(5,5,4,2))
plot.manhattan2(pos = allSNPs.pc_load$POS,
                chr = allSNPs.pc_load$CHROM,
                y = list(abs(allSNPs.pc_load$PC1), abs(allSNPs.pc_load$PC2)), col = col, log10 = F)
legend('topleft', c('PC1', 'PC2'), col = col, pch = 19, cex = 3)
mtext('abs(PC loading)', side = 2, cex = 2, line = 2.7)
dev.off()
#more pvals
png('../figures/manhattan_pInt_pTime.png', 1600, 1200)
par(mar = c(5,5,4,2))
plot.manhattan2(pos = binGlmScan_contTime_param1234_200615$POS,
                chr = binGlmScan_contTime_param1234_200615$CHROM,
                y = list(binGlmScan_contTime_param1234_200615$generation_p, binGlmScan_contTime_param1234_200615$`generation:treatment_p`),
                col = col, log10 = T)
legend('topleft', c(expression(p[time]), expression(p[int])), col = col, pch = 19, cex = 3)
mtext('-log10(p)', side = 2, cex = 2, line = 2.7)
dev.off()
#combined
pTime <- log10(binGlmScan_contTime_param1234_200615$generation_p)
pInt <- log10(binGlmScan_contTime_param1234_200615$`generation:treatment_p`)
pc1 <- abs(allSNPs.pc_load$PC1)
pc1 <- pc1*max(abs(pInt))/max(pc1)
pc2 <- abs(allSNPs.pc_load$PC2)
pc2 <- pc2*max(abs(pInt))/max(pc2)
ylim <- c(-60, 60)
png('../figures/manhattan_allSNPsPC12_loadings_pInt_pTime_v2.png', 1600, 1200)
par(mar = c(5,5,4,2))
plot.manhattan2(pos = allSNPs.pc_load$POS,
                chr = allSNPs.pc_load$CHROM,
                y = list(pc1, pc2,
                         pTime, pInt),
                col = c(col, col), log10 = F, ylim = ylim)
legend('topleft', c('PC1', 'PC2'), col = col, pch = 19, cex = 3)
legend('bottomleft', c(expression(p[time]), expression(p[int])), col = col, pch = 19, cex = 3)
mtext('abs(PC loading)', side = 2, cex = 2, line = 2.7, adj = .7)
mtext('-log10(p)', side = 2, cex = 2, line = 2.7, adj = .3)
abline(h = 0, lwd = 10, lty = 1, col = 'white')
dev.off()


#Hm
binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$emTrend.C_p < 1e-45
                                     & binGlmScan_contTime_param1234_200615$CHROM == '2L',]
p <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2L' & freqs.merge_union_190426$POS == 17559090, ]
plot.freqVStime_v2(p[, 9:32], lines = T, ylim = 0:1)
p <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2L' & freqs.merge_union_190426$POS == 20944347, ]
plot.freqVStime_v2(p[, 9:32], lines = T, ylim = 0:1)
plot.manhattan2(pos = binGlmScan_contTime_param1234_200615$POS,
                chr = binGlmScan_contTime_param1234_200615$CHROM,
                y = list(p_c, p_hs, binGlmScan_contTime_param1234_200615$`generation:treatment_p`),
                col = c('blue', 'red', 'black'), ylim = ylim, zoom.chr = '3L', main = '3L')
binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$emTrend.C_p < 1e-35
                                     & binGlmScan_contTime_param1234_200615$CHROM == '3L'
                                     & binGlmScan_contTime_param1234_200615$POS < 6e6,]
p <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3L' & freqs.merge_union_190426$POS == 5647345, ]
plot.freqVStime_v2(p[, 9:32], lines = T, ylim = 0:1)
binGlmScan_contTime_param1234_200615[(binGlmScan_contTime_param1234_200615$emTrend.C_p < 1e-35 | binGlmScan_contTime_param1234_200615$`generation:treatment_p` < 1e-30)
                                     & binGlmScan_contTime_param1234_200615$CHROM == '3L'
                                     & binGlmScan_contTime_param1234_200615$POS < 6e6,]
p <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3L' & freqs.merge_union_190426$POS == 2288829, ]
plot.freqVStime_v2(p[, 9:32], lines = T, ylim = 0:1)