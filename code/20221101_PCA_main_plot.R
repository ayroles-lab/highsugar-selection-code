source("code/utils.R")
# 2020-11-25
# Fig 1a
load('data/201001_alleFreqPC_allSNPs.RData')
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

eig <- allSNPs.pc$sdev^2
varExpl <- eig/sum(eig)
dat_var <- data.frame(pc = 1:length(varExpl), varExpl = varExpl)

#Trying with the var expl as inset
# col.treat <- c(wes_palette('Darjeeling1', 1), 'deepskyblue2')
#col.treat <- c("#B3B3FF", "#FFD9B3")
col.treat <- c(2, 1)

vp <- viewport(width = 0.25, height = 0.33, x = 0.25, y = 0.8)

pc <- ggplot(dat, aes(x = PC1, y = PC2, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("First principal component (", round(dat_var$varExpl[1], 3)*100, "% variance explained)"), 
       y = paste0("Second principal component (", round(dat_var$varExpl[2], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = c(0.9, 0.5)) + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 


pc.var <- ggplot(dat_var, aes(x = pc, y = varExpl, group = 1)) + geom_bar(stat = "identity") +
  xlab('PC') + ylab('Variance Explained') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        strip.text.y = element_text(size = 8), 
        axis.line = element_line(color = 'black'))

png('figures/PCs12_inset_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc)
print(pc.var, vp = vp)
dev.off()

pc23 <- ggplot(dat, aes(x = PC2, y = PC3, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("Second principal component (", round(dat_var$varExpl[2], 3)*100, "% variance explained)"), 
       y = paste0("Third principal component (", round(dat_var$varExpl[3], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = "bottom") + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 
png('figures/PCs23_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc23)
dev.off()

pc13 <- ggplot(dat, aes(x = PC1, y = PC3, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("First principal component (", round(dat_var$varExpl[1], 3)*100, "% variance explained)"), 
       y = paste0("Third principal component (", round(dat_var$varExpl[3], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = "bottom") + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 
png('figures/PCs13_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc13)
dev.off()
pc13 <- ggplot(dat, aes(x = PC1, y = PC3, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("First principal component (", round(dat_var$varExpl[1], 3)*100, "% variance explained)"), 
       y = paste0("Third principal component (", round(dat_var$varExpl[3], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = "bottom") + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 
png('figures/PCs13_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc13)
dev.off()

pc14 <- ggplot(dat, aes(x = PC1, y = PC4, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("First principal component (", round(dat_var$varExpl[1], 3)*100, "% variance explained)"), 
       y = paste0("Fourth principal component (", round(dat_var$varExpl[4], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = "bottom") + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 
png('figures/PCs14_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc14)
dev.off()

pc15 <- ggplot(dat, aes(x = PC1, y = PC5, group = pop, color = treat)) +
  geom_line(linewidth = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x =  paste0("First principal component (", round(dat_var$varExpl[1], 3)*100, "% variance explained)"), 
       y = paste0("Fifth principal component (", round(dat_var$varExpl[5], 3)*100, "% variance explained)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = "bottom") + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 
png('figures/PCs15_allSNPs.png', width = 6, height = 4.5, units = "in", res = 300)
print(pc15)
dev.off()

}

pc.poster <- ggplot(dat, aes(x = PC1, y = PC2, group = pop, color = treat)) +
  geom_line(size = 1, key_glyph = draw_key_rect) + 
  geom_point(size = 4, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  labs(x = "First principal component", y = "Second principal component") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text.y = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        legend.position = c(0.9, 0.5)) + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 

save_plot("figures/PCs-poster.pdf", pc.poster ,
          base_width = 7.8, base_height = 5 )
save_plot("figures/PCs-poster.png", pc.poster ,
          base_width = 7.8, base_height = 5 )
### PCs SI figure

load('200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
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
