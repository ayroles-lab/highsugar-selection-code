source("utils.R")
# 2020-11-25
# Fig 1a
load('../data/201001_alleFreqPC_allSNPs.RData')
summary(allSNPs.pc)

library(ggplot2)
library(grid)
library(wesanderson)
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
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_text(size=35,face="bold"),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20)) + #legend.justification = c(1,.2)
  guides(colour=guide_legend(title="Treatment"), shape = guide_legend(title="Generation"))

eig <- allSNPs.pc$sdev^2
varExpl <- eig/sum(eig)
dat <- data.frame(pc = 1:length(varExpl), varExpl = varExpl)
pc.var <- ggplot(dat, aes(x = pc, y = varExpl, group = 1)) + geom_bar(stat = "identity") +
  xlab('PC') + ylab('Variance Explained') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"),
        strip.text.y = element_text(size = 20))

png('../figures/PCs12_inset_allSNPs.png', width = 1600, height = 1000)
print(pc)
print(pc.var, vp = vp)
dev.off()

########## 210128 #########
# Fig 1b
library(wesanderson)
load(
    "../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData"
)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in%
    chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-04
signCut <- 1e-09
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
col <- rep(
    wes_palette("Darjeeling2")[2],
    nrow(binGlmScan_contTime_param1234_200615)
)
col[binGlmScan_contTime_param1234_200615$`generation:treatment_p` > signCut] <- "grey"
png("../figures/manhattan_pInt_HSonly_cut1e-9.png", 1600, 1200)
par(mar = c(5, 5, 4, 2))
plot.manhattan3(
    pos = binGlmScan_contTime_param1234_200615$POS, chr = binGlmScan_contTime_param1234_200615$CHROM,
    y = p, col = col, ylim = ylim
)
abline(
    h = -log10(signCut),
    lty = 2, lwd = 2
)
mtext(
    expression("-log10(" ~ p[int] ~ ")"),
    side = 2, cex = 2, line = 2.7
)
dev.off()
