#Check overlap between IHS and regression results
library(wesanderson)
source("utils.R")
load("../data/191008_ihsScans_bothTreat.RData")
load("../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData")

#Combine
identical(wgscan.ihs_HS$ihs$POSITION, wgscan.ihs_C$ihs$POSITION)
snps_HS <- paste(wgscan.ihs_HS$ihs$CHR, wgscan.ihs_HS$ihs$POSITION, sep = '_')
snps_C <- paste(wgscan.ihs_C$ihs$CHR, wgscan.ihs_C$ihs$POSITION, sep = '_')
ihs_HS <- wgscan.ihs_HS$ihs[match(snps_C, snps_HS), ]
sum(wgscan.ihs_C$ihs$POSITION - ihs_HS$POSITION, na.rm = T)

col <- wes_palette('Darjeeling1', n = 2)
png('../figures/ihsAbs_manhattan_comb.png', width = 1200, height = 1000)
plot.manhattan2(pos = wgscan.ihs_C$ihs$POSITION,
                chr = wgscan.ihs_C$ihs$CHR,
                y = list(abs(wgscan.ihs_C$ihs$IHS), abs(ihs_HS$IHS)), log10 = F, col = col)
mtext('abs(IHS)', side = 2, line = 1.8, cex = 2)
legend('topleft', c('Control', 'HS'), col = col, pch = 19, cex = 2)
dev.off()

png('../figures/ihsLogP_manhattan_comb.png', width = 1200, height = 1000)
plot.manhattan2(pos = wgscan.ihs_C$ihs$POSITION,
                chr = wgscan.ihs_C$ihs$CHR,
                y = list(wgscan.ihs_C$ihs$LOGPVALUE, ihs_HS$LOGPVALUE), log10 = F)
mtext('-log10(p)', side = 2, line = 1.8, cex = 1.5)
legend('topleft', c('Control', 'HS'), col = col, pch = 19, cex = 2)
dev.off()
# wgscan.ihs$ihs[order(abs(wgscan.ihs$ihs$IHS), decreasing = T)[1:10], ]

# ptm <- Sys.time()
# HS.3L <- data2haplohh(hap_file = "../data/2019-09-30_currentGen_vcfSplit/filteredSNPset_190419_3R_HS.vcf.bgz",
#                          polarize_vcf = FALSE, min_perc_geno.mrk = 70, min_perc_geno.hap = 50)
# Sys.time() - ptm

# snp <- which(HS.3L@positions == 24624973)
# hh_top3L.ehs <- calc_ehh(HS.3L, mrk = snp, phased = F, polarized = F, include_nhaplo = T)
# plot(hh_top3L.ehs)


#Check overlap between IHS and regression results
# load('190429_binGlmScan_contTime_param1234_SNPunion_update.RData')
snps_HS <- paste(wgscan.ihs_HS$ihs$CHR, wgscan.ihs_HS$ihs$POSITION, sep = '_')
snps_C <- paste(wgscan.ihs_C$ihs$CHR, wgscan.ihs_C$ihs$POSITION, sep = '_')
snps_binGlmScan <- paste(binGlmScan_contTime_param1234_190429$CHROM, binGlmScan_contTime_param1234_190429$POS, sep = '_')
#Control
binGlmScan.match <- binGlmScan_contTime_param1234_190429[match(snps_C, snps_binGlmScan), ]
sum(binGlmScan.match$POS - wgscan.ihs_C$ihs$POSITION, na.rm = T)
png('../figures/ihsControl_VS_inP.png', width = 1000, height = 1000)
plot(-log10(binGlmScan.match$`generation:treatment_p`), abs(wgscan.ihs_C$ihs$IHS),
     xlab = '-log10(p_int)', ylab = 'IHS Control', cex.lab = 1.5, xlim = c(0, 61), ylim = c(0, 7.8))
dev.off()
logP <- -log10(binGlmScan.match$`generation:treatment_p`)
summary(lm(abs(wgscan.ihs_C$ihs$IHS) ~ logP))
#HS
binGlmScan.match <- binGlmScan_contTime_param1234_190429[match(snps_HS, snps_binGlmScan), ]
sum(binGlmScan.match$POS - wgscan.ihs_HS$ihs$POSITION, na.rm = T)
png('../figures/ihsHS_VS_inP.png', width = 1000, height = 1000)
plot(-log10(binGlmScan.match$`generation:treatment_p`), abs(wgscan.ihs_HS$ihs$IHS),
     xlab = '-log10(p_int)', ylab = 'IHS HS', cex.lab = 1.5, xlim = c(0, 61), ylim = c(0, 7.8))
dev.off()
logP <- -log10(binGlmScan.match$`generation:treatment_p`)
summary(lm(abs(wgscan.ihs_HS$ihs$IHS) ~ logP))
