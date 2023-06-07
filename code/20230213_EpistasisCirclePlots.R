########## 210126 #########
#More stringent set of candidate pairs
# load('210126_chisqLinkageClumps_distantChisqP5.7e-08_corP0.01_clumpDist250000_clumpChisqP1e-05.RData')
# load('210126_chisqLinkageClumps_distantChisqP5.7e-08_corP1e-05_clumpDist250000_clumpChisqP1e-05.RData')
load('../archive/R/210126_chisqLinkageClumps_distantChisqP5.7e-08_corP0.001_clumpDist250000_clumpChisqP1e-05.RData')
load('../data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')

#Pick the clumps with > n clumped SNPs at both ends
n <- 2
clumpsOfInterest <- list()
clumpsOfInterest_bySNP1 <- list()
for (i in 1:length(chisqLinkage_phasedGenos_snp1ANDsnp2Clumps)) {
  snp1Clumps <- chisqLinkage_phasedGenos_snp1ANDsnp2Clumps[[i]]
  snp1.name <- names(chisqLinkage_phasedGenos_snp1ANDsnp2Clumps)[i]

  tmp <- list()
  if(length(snp1Clumps) > 0){
    for (j in 1:length(snp1Clumps)) {
      clump <- snp1Clumps[[j]]
      snp2.name <- names(snp1Clumps)[j]

      snps1 <- unique(as.character(clump$snp1))
      snps2 <- unique(as.character(clump$snp2))
      if(length(snps1) > n & length(snps2) > n){
        clumpsOfInterest[[paste(snp1.name, snp2.name, sep = '__')]] <- clump
        tmp[[snp2.name]] <- clump
      }
    }
  }
  if(length(tmp) > 0)
    clumpsOfInterest_bySNP1[[snp1.name]] <- tmp
}
#save(clumpsOfInterest_bySNP1, file = paste0('210126_chisqLinkageClumps_moreThan', n, 'SNPsPerEnd.RData'))

#Circos
library(RCircos)
library(GenomicRanges)
library(RColorBrewer)
library(tidyverse)

chrs <- c('2L', '2R', '3L', '3R')
binGlm.chrs <- as.character(rep.int(binGlmScan_contTime_param1234_200615.gr@seqnames@values, binGlmScan_contTime_param1234_200615.gr@seqnames@lengths))
binGlmScan_contTime_param1234_200615.gr <- binGlmScan_contTime_param1234_200615.gr[binGlm.chrs %in% chrs]
binGlmScan_contTime_param1234_200615.gr$chr <- as.character(rep.int(binGlmScan_contTime_param1234_200615.gr@seqnames@values, binGlmScan_contTime_param1234_200615.gr@seqnames@lengths))
binGlmScan_chrs.gr <- unlist(reduce(split(binGlmScan_contTime_param1234_200615.gr, ~chr), min.gapwidth = 1e6))
myGenome.rcircos <- data.frame(Chromosome = binGlmScan_chrs.gr@seqnames@values,
                               ChromStart = 0, #binGlmScan_chrs.gr@ranges@start,
                               ChromEnd = binGlmScan_chrs.gr@ranges@start + binGlmScan_chrs.gr@ranges@width,
                               Band = c('p36.33', 'p36.32', 'p36.31', 'p36.23'), #I don't know what this is. Just took it from data(UCSC.HG19.Human.CytoBandIdeogram)
                               Stain = c('gneg', 'gpos25', 'gneg', 'gpos25'))

# pairsPerClump <- sapply(chisqLinkage_phasedGenos_snp1Clumps, nrow)
# clumps <- chisqLinkage_phasedGenos_snp1Clumps[pairsPerClump > 6]
# clumps <- clumpsOfInterest
clumps <- clumpsOfInterest_bySNP1

# col <- colorRampPalette(brewer.pal(n = 9, 'Set3'))(length(clumps))
col <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(clumps))
clump <- bind_rows(clumps[[1]])
linkData <- data.frame(Chromosome = clump$snp1.chr,
                       chromStart = clump$snp1.pos,
                       chromEnd = clump$snp1.pos + 1,
                       Chromosome.1 = clump$snp2.chr,
                       chromStart.1 = clump$snp2.pos,
                       chromEnd.1 = clump$snp2.pos + 1,
                       PlotColor =  col[1])
for (i in 2:length(clumps)) {
  clump <- bind_rows(clumps[[i]])

  tmp <- data.frame(Chromosome = clump$snp1.chr,
                    chromStart = clump$snp1.pos,
                    chromEnd = clump$snp1.pos + 1,
                    Chromosome.1 = clump$snp2.chr,
                    chromStart.1 = clump$snp2.pos,
                    chromEnd.1 = clump$snp2.pos + 1,
                    PlotColor = col[i])
  linkData <- rbind(linkData, tmp)
}

# png(paste('../results/2021-01-26_newClumps/circos_distantChisqP5.7e-08_corP0.001_clumpDist250000_clumpChisqP1e-05_atleast3snps.png'), 1200, 1200)
# png(paste('../results/2021-01-26_newClumps/circos_distantChisqP5.7e-08_corP0.01_clumpDist250000_clumpChisqP1e-05_atleast5snps.png'), 1200, 1200)
# png(paste('../results/2021-01-26_newClumps/circos_distantChisqP5.7e-08_corP1e-05_clumpDist250000_clumpChisqP1e-05_atleast2snps.png'), 1200, 1200)

png(paste('../figures/2023-02-13_newClumps-circos_distantChisqP5.7e-08_corP0.001_clumpDist250000_clumpChisqP1e-05_atleast3snps_v3.png'), 1000, 1000)
par(mar=c(0, 0, 0, 0))
RCircos.Set.Core.Components(myGenome.rcircos, tracks.inside = 1, tracks.outside = 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(link.data = linkData, track.num = 1, by.chromosome = F, lineWidth = rep(3, nrow(linkData)))
dev.off()
