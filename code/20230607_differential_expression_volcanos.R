#Volcano plots
options(stringsAsFactors = F)
swapExp_DESeqRes_body <- read.table('data/Results_DESeq_SWAPexperiment_body_May12.txt', header = T)
swapExp_DESeqRes_head <- read.table('data/Results_DESeq_SWAPexperiment_head_May12.txt', header = T)

#body
DE.fdr <- .01

#ggplot facet version
dat <- data.frame(tissue = rep(c('Body', "Head"), 
                               c(nrow(swapExp_DESeqRes_body), 
                                 nrow(swapExp_DESeqRes_head))),
                  lfc_geno = c(swapExp_DESeqRes_body$lfc_geno, 
                               swapExp_DESeqRes_head$lfc_geno),
                  p_geno = c(swapExp_DESeqRes_body$padj_geno,                
                             swapExp_DESeqRes_head$padj_geno))
dat$sign <- 0
dat$sign[dat$p_geno < DE.fdr] <- 1


g <- ggplot(data = dat, aes(x = lfc_geno, y = -log10(p_geno), color = sign)) +
  geom_point(size = 0.8) + theme_classic() +
  facet_grid(cols = vars(tissue)) +
  # scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8),
        strip.text.x = element_text(size = 8)) +
  xlab('Log2 fold change') + ylab(expression('-log10(' ~ p[DE] ~ ')')) +
  geom_hline(yintercept = -log10(DE.fdr), linetype = 'dashed')

png('figures/DEgeno_bodyANDhead_volcano_ggplot_v2.png', 
    width = 3, height = 3, unit = "in", res = 300)
print(g)
dev.off()

load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_w5000.RData')
load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_permBody_w5000.RData')
load('archive/R/210209_swapExp_marginalDEfdr01_varyingPint_overlap_cOnlyExcluded_permHead_w5000.RData')
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(dplyr)

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
            group = 1, size = 1, color = 2) + 
  facet_grid(cols = vars(tissue)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) + ylim(ylim) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) 
  
save_plot("figures/DE.png", g + ggtitle("A.") + g.geno_facet + ggtitle("B.") ,
          base_width = 5.5, base_height = 3 )
