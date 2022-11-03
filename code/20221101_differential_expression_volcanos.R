#Volcano plots
options(stringsAsFactors = F)
swapExp_DESeqRes_body <- read.table('../data/Results_DESeq_SWAPexperiment_body_May12.txt', header = T)
swapExp_DESeqRes_head <- read.table('../data/Results_DESeq_SWAPexperiment_head_May12.txt', header = T)

#body
DE.fdr <- .01
col <- rep('gray', nrow(swapExp_DESeqRes_body))
col[swapExp_DESeqRes_body$padj_geno < DE.fdr] <- 'red'
png('../figures/DEgeno_body_volcano.png', 1200, 1200)
par(mar = c(5,5,4,2))
plot(swapExp_DESeqRes_body$lfc_geno, -log10(swapExp_DESeqRes_body$padj_geno), col = col, pch = 19, cex = 1,
     xlab = 'Log2 fold change', ylab = '-log10(p)', cex.lab = 2.5, cex.axis = 2)
abline(h = -log10(DE.fdr), lty = 2, lwd = 2)
dev.off()

#head
DE.fdr <- .01
col <- rep('gray', nrow(swapExp_DESeqRes_head))
col[swapExp_DESeqRes_head$padj_geno < DE.fdr] <- 'red'
png('../figures/DEgeno_head_volcano.png', 1200, 1200)
par(mar = c(5,5,4,2))
plot(swapExp_DESeqRes_head$lfc_geno, -log10(swapExp_DESeqRes_head$padj_geno), col = col, pch = 19, cex = 1,
     xlab = 'Log2 fold change', ylab = '-log10(p)', cex.lab = 2.5, cex.axis = 2)
abline(h = -log10(DE.fdr), lty = 2, lwd = 2)
dev.off()

#Combined
xlim <- range(c(swapExp_DESeqRes_body$lfc_geno, swapExp_DESeqRes_head$lfc_geno))
ylim <- range(-log10(c(swapExp_DESeqRes_body$padj_geno, swapExp_DESeqRes_head$padj_geno)))

png('../figures/DEgeno_bodyANDhead_volcano.png', width = 2*1200, height = 1200)
par(mar = c(5,5,4,2), mfrow = 1:2)

col <- rep('gray', nrow(swapExp_DESeqRes_body))
col[swapExp_DESeqRes_body$padj_geno < DE.fdr] <- 'red'
plot(swapExp_DESeqRes_body$lfc_geno, -log10(swapExp_DESeqRes_body$padj_geno), col = col, pch = 19, cex = 1,
     xlab = 'Log2 fold change', ylab = '-log10(p)', cex.lab = 2.5, cex.axis = 2, main = 'Body', cex.main = 3,
     xlim = xlim, ylim = ylim)
abline(h = -log10(DE.fdr), lty = 2, lwd = 2)

col <- rep('gray', nrow(swapExp_DESeqRes_head))
col[swapExp_DESeqRes_head$padj_geno < DE.fdr] <- 'red'
plot(swapExp_DESeqRes_head$lfc_geno, -log10(swapExp_DESeqRes_head$padj_geno), col = col, pch = 19, cex = 1,
     xlab = 'Log2 fold change', ylab = '-log10(p)', cex.lab = 2.5, cex.axis = 2, main = 'Head', cex.main = 3,
     xlim = xlim, ylim = ylim)
abline(h = -log10(DE.fdr), lty = 2, lwd = 2)
dev.off()

#ggplot facet version
dat <- data.frame(tissue = rep(c('Body', "Head"), c(nrow(swapExp_DESeqRes_body), nrow(swapExp_DESeqRes_head))),
                  lfc_geno = c(swapExp_DESeqRes_body$lfc_geno, swapExp_DESeqRes_head$lfc_geno),
                  p_geno = c(swapExp_DESeqRes_body$padj_geno, swapExp_DESeqRes_head$padj_geno))
dat$sign <- 0
dat$sign[dat$p_geno < DE.fdr] <- 1

g <- ggplot(data = dat, aes(x = lfc_geno, y = -log10(p_geno), color = sign)) +
  geom_point(size = 4) + theme_classic() +
  facet_grid(cols = vars(tissue)) +
  # scale_color_manual(values = c('0' = 'black', '1' = 'red')) +
  theme(legend.position = 'none',
        axis.text=element_text(size=35),
        axis.title=element_text(size=40),
        plot.title = element_text(hjust = 0.5, size = 25),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        strip.text.x = element_text(size = 40, face = 'bold')) +
  xlab('Log2 fold change') + ylab(expression('-log10(' ~ p[DE] ~ ')')) +
  geom_hline(yintercept = -log10(DE.fdr), linetype = 'dashed')

png('../figures/DEgeno_bodyANDhead_volcano_ggplot_v2.png', width = 1200, height = 1200)
print(g)
dev.off()
