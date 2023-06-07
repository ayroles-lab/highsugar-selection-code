
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

load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_w5000.RData')
load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_permBody_w5000.RData')
load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_permHead_w5000.RData')
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(dplyr)
#Fraction overlapping SNPs
ylim_body <- range(c(overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nOverlapSNPs_geno.perm/      
                overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nSNPs,

                overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nOverlapSNPs_geno/overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nSNPs,

                overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nOverlapSNPs_diet.perm/overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nSNPs,

                overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nOverlapSNPs_diet/overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm$nSNPs))
#Plot body
g.geno_b <- ggplot(data = overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm, 
                 aes(y = nOverlapSNPs_geno.perm/nSNPs, 
                     x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr01_varyingSelP_body, 
            mapping = aes(y = nOverlapSNPs_geno/nSNPs, 
                          x = factor(-log10(p))), 
            group = 1, linewidth = 1) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5, size = 8)) + ylim(ylim_body) +
  ylab('Fraction of selected SNPs in DE genes') + 
  xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Selection History Effect')

#Plot head
ylim_head <- range(c(overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nOverlapSNPs_geno.perm/overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nOverlapSNPs_geno/overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nOverlapSNPs_diet.perm/overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nSNPs,
                overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nOverlapSNPs_diet/overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm$nSNPs))

g.geno_h <- ggplot(data = overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm, 
                 aes(y = nOverlapSNPs_geno.perm/nSNPs, x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = overlap_swapExpSelection_DEfdr01_varyingSelP_head,
            mapping = aes(y = nOverlapSNPs_geno/nSNPs, x = factor(-log10(p))), 
            group = 1, size = 1) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5, size = 8)) + ylim(ylim_head) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) +
  ggtitle('Selection History Effect - HEAD')

save_plot("test.png", g.geno_h + g.geno_b)

ylim = c(0.06, 0.22)

dat.perm = bind_rows(Head = overlap_swapExpSelection_DEfdr01_varyingSelP_head.perm, 
                     Body = overlap_swapExpSelection_DEfdr01_varyingSelP_body.perm, 
                     .id = "tissue")
dat.obs = bind_rows(Head = overlap_swapExpSelection_DEfdr01_varyingSelP_head, 
                     Body = overlap_swapExpSelection_DEfdr01_varyingSelP_body, 
                     .id = "tissue")
g.geno_facet <- ggplot(data = dat.perm, 
                       aes(y = nOverlapSNPs_geno.perm/nSNPs,
                           x = factor(-log10(p)))) +
  geom_boxplot(key_glyph = 'rect') +
  geom_line(data = dat.obs,
            mapping = aes(y = nOverlapSNPs_geno/nSNPs, x = factor(-log10(p))), 
            group = 1, size = 1) + 
  facet_grid(cols = vars(tissue)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5, size = 8)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) 
  
save_plot("test.png", g.geno_facet )


