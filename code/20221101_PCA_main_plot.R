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

#Trying with the var expl as inset
col.treat <- c(wes_palette('Darjeeling1', 1), 'deepskyblue2')
vp <- viewport(width = 0.25, height = 0.33, x = 0.205, y = 0.83)

pc <- ggplot(dat, aes(x = PC1, y = PC2, group = pop, colour = treat)) +
  geom_line(size = 2, key_glyph = 'crossbar') + geom_point(size = 8, aes(shape = factor(t))) +
  scale_shape_manual(values = 15:18) +
  scale_colour_manual(name = "treat", values = col.treat) +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=35,face="bold"),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20),
        axis.line = element_line(color = 'black'),
        legend.position = c(0.8, 0.5)) + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation")) 

eig <- allSNPs.pc$sdev^2
varExpl <- eig/sum(eig)
dat <- data.frame(pc = 1:length(varExpl), varExpl = varExpl)
pc.var <- ggplot(dat, aes(x = pc, y = varExpl, group = 1)) + geom_bar(stat = "identity") +
  xlab('PC') + ylab('Variance Explained') +
  theme_tufte() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20), 
        axis.line = element_line(color = 'black'))

png('../figures/PCs12_inset_allSNPs.png', width = 1600, height = 1000)
print(pc)
print(pc.var, vp = vp)
dev.off()
}