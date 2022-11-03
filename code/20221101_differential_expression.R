
########## 200904 #########
#Gene lists for Julien
load('../data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')
load('../data/200608_DEres_GRanges.RData')
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
write.table(dat, file = '../figures/genes_DEgenoORdietFDR05_selP8e12_body', quote = F, sep = '\t', row.names = F)

DE.head <- dmel.gff_swapGenes_head[dmel.gff_swapGenes_head$padj_geno < .05 | dmel.gff_swapGenes_head$padj_diet < .05]
overlap.head <- findOverlaps(selSNPs, DE.head)
DE.head_overlap <- DE.head[unique(overlap.head@to)]
dat <- as.data.frame(DE.head_overlap[,c(1:5, 61:74)])
write.table(dat, file = '../figures/genes_DEgenoORdietFDR05_selP8e12_head', quote = F, sep = '\t', row.names = F)

#Figure 4 (or 5 depending on the paper structure)
#Count DE genes & overlaps
body.DEgeno <- dmel.gff_swapGenes_body$ID[dmel.gff_swapGenes_body$padj_geno < .05]
body.DEdiet <- dmel.gff_swapGenes_body$ID[dmel.gff_swapGenes_body$padj_diet < .05]
length(body.DEgeno)
length(body.DEdiet)
sum(body.DEgeno %in% body.DEdiet)

head.DEgeno <- dmel.gff_swapGenes_head$ID[dmel.gff_swapGenes_head$padj_geno < .05]
head.DEdiet <- dmel.gff_swapGenes_head$ID[dmel.gff_swapGenes_head$padj_diet < .05]
length(head.DEgeno)
length(head.DEdiet)
sum(head.DEgeno %in% head.DEdiet)

load('../data/200904_swapExp_marginalDE_varyingPint_overlap_cOnlyExcluded.RData')
load('../data/200904_swapExp_marginalDE_varyingPint_overlap_cOnlyExcluded_permBody.RData')
load('../data/200904_swapExp_marginalDE_varyingPint_overlap_cOnlyExcluded_permHead.RData')
library(ggplot2)
library(ggpubr)
#Fraction overlapping SNPs
ylim <- range(c(overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nOverlapSNPs_geno.perm/overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nOverlapSNPs_geno/overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nOverlapSNPs_diet.perm/overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nOverlapSNPs_diet/overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm$nSNPs))
#Plot body
g.geno <- ggplot(data = overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm, aes(y = nOverlapSNPs_geno.perm/nSNPs, x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr05_varyingSelP_body, mapping = aes(y = nOverlapSNPs_geno/nSNPs, x = factor(-log10(p)), group = 1, size = 1)) + #x = p_int, y = overlap_HS,
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 15)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Selection History Effect')

g.diet <- ggplot(data = overlap_swapExpSelection_DEfdr05_varyingSelP_body.perm, aes(y = nOverlapSNPs_diet.perm/nSNPs, x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr05_varyingSelP_body, mapping = aes(y = nOverlapSNPs_diet/nSNPs, x = factor(-log10(p)), group = 1, size = 1)) + #x = p_int, y = overlap_HS,
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 15)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Plastic Effect')

png('../figures/fig_4_swapMarginalDE_pInt_overlap_body.png', width = 1200, height = 800)
ggarrange(g.geno, g.diet)
dev.off()

#Plot head
ylim <- range(c(overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nOverlapSNPs_geno.perm/overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nOverlapSNPs_geno/overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nOverlapSNPs_diet.perm/overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nOverlapSNPs_diet/overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm$nSNPs))

g.geno <- ggplot(data = overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm, aes(y = nOverlapSNPs_geno.perm/nSNPs, x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr05_varyingSelP_head, mapping = aes(y = nOverlapSNPs_geno/nSNPs, x = factor(-log10(p)), group = 1, size = 1)) + #x = p_int, y = overlap_HS,
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 15)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Selection History Effect')

g.diet <- ggplot(data = overlap_swapExpSelection_DEfdr05_varyingSelP_head.perm, aes(y = nOverlapSNPs_diet.perm/nSNPs, x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr05_varyingSelP_head, mapping = aes(y = nOverlapSNPs_diet/nSNPs, x = factor(-log10(p)), group = 1, size = 1)) + #x = p_int, y = overlap_HS,
  theme(legend.position = 'none',
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 15)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Plastic Effect')

png('../figures/fig_4_swapMarginalDE_pInt_overlap_head.png', width = 1200, height = 800)
ggarrange(g.geno, g.diet)
dev.off()

#Plastic vs Evolved DE effects. Like fig2 here https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006336
load('../data/200608_DEres_GRanges.RData')

keep <- dmel.gff_swapGenes_body$padj_diet < .05 & dmel.gff_swapGenes_body$padj_geno < .05
png('../figures/fig_4_plastic_vs_evol_DE_body.png', 1200, 1200)
par(mar = c(5,6,4,2) + .1)
plot(dmel.gff_swapGenes_body$lfc_diet[keep], dmel.gff_swapGenes_body$lfc_geno[keep],
     xlab = 'Plastic Effect', ylab = 'Selection History Effect', cex.lab = 3, cex.axis = 2, main = '', cex.main = 3)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
grid()
dev.off()

keep <- dmel.gff_swapGenes_head$padj_diet < .05 & dmel.gff_swapGenes_head$padj_geno < .05
png('../figures/fig_4_plastic_vs_evol_DE_head.png', 1200, 1200)
par(mar = c(5,6,4,2) + .1)
plot(dmel.gff_swapGenes_head$lfc_diet[keep], dmel.gff_swapGenes_head$lfc_geno[keep],
     xlab = 'Plastic Effect', ylab = 'Selection History Effect', cex.lab = 3, cex.axis = 2, main = '', cex.main = 3)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
grid()
dev.off()
