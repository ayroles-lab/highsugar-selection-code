
########## 200904 #########
#Gene lists for Julien
load('data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')
load('data/200608_DEres_GRanges.RData')
chrs <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
binGlm.chrs <- as.character(rep.int(binGlmScan_contTime_param1234_200615.gr@seqnames@values, binGlmScan_contTime_param1234_200615.gr@seqnames@lengths))
binGlmScan_contTime_param1234_200615.gr <- binGlmScan_contTime_param1234_200615.gr[binGlm.chrs %in% chrs]
nonSignCut <- 1e-4
p <- binGlmScan_contTime_param1234_200615.gr$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615.gr$emTrend.HS_p > nonSignCut] <- 1
selSNPs <- binGlmScan_contTime_param1234_200615.gr[p < 8e-12]

#Geno or Diet DE
DE.body <- dmel.gff_swapGenes_body[dmel.gff_swapGenes_body$padj_geno < .05 | dmel.gff_swapGenes_body$padj_diet < .05]
overlap.body <- findOverlaps(selSNPs, DE.body)
DE.body_overlap <- DE.body[unique(overlap.body@to)]
dat <- as.data.frame(DE.body_overlap[,c(1:5, 61:74)])
#write.table(dat, file = '../figures/genes_DEgenoORdietFDR05_selP8e12_body', quote = F, sep = '\t', row.names = F)

DE.head <- dmel.gff_swapGenes_head[dmel.gff_swapGenes_head$padj_geno < .05 | dmel.gff_swapGenes_head$padj_diet < .05]
overlap.head <- findOverlaps(selSNPs, DE.head)
DE.head_overlap <- DE.head[unique(overlap.head@to)]
dat <- as.data.frame(DE.head_overlap[,c(1:5, 61:74)])
#write.table(dat, file = '../figures/genes_DEgenoORdietFDR05_selP8e12_head', quote = F, sep = '\t', row.names = F)

#Figure 4 (or 5 depending on the paper structure)
#Count DE genes & overlaps
body.DEgeno <- dmel.gff_swapGenes_body$ID[dmel.gff_swapGenes_body$padj_geno < .01]
body.DEdiet <- dmel.gff_swapGenes_body$ID[dmel.gff_swapGenes_body$padj_diet < .01]
length(body.DEgeno)
length(body.DEdiet)
sum(body.DEgeno %in% body.DEdiet)

head.DEgeno <- dmel.gff_swapGenes_head$ID[dmel.gff_swapGenes_head$padj_geno < .01]
head.DEdiet <- dmel.gff_swapGenes_head$ID[dmel.gff_swapGenes_head$padj_diet < .01]
length(head.DEgeno)
length(head.DEdiet)
sum(head.DEgeno %in% head.DEdiet)
