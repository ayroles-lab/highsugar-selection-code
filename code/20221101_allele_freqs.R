library(wesanderson)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(grid)

#Fig 3 with these SNPs
w <- 800
h <- 1000
load('../data/200917_meanFreqDiffs.RData')
nonSignCut <- 1e-4
signCut <- 8e-12 #1e-9
p_int <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p_HS <- binGlmScan_contTime_param1234_200615$emTrend.HS_p

#Distribution
#Selected SNPs
idx <- which(p_int < signCut & p_HS < nonSignCut)
tmp <- abs(freqDiff_200917[idx, 9:11])
colnames(tmp) <- c('G1 → G11', 'G1 → G25', 'G1 → G100')
dat <- pivot_longer(tmp, c('G1 → G11', 'G1 → G25', 'G1 → G100'))
dat$name <- factor(dat$name, levels = c('G1 → G11', 'G1 → G25', 'G1 → G100'))
library(ggthemes)
library(viridis)

pal = c("#FEE0D2", "#E69F00", "#DE2D26")
#plot distribution
p2 <- ggplot(dat, aes(x=value, fill=name)) + geom_histogram(alpha=1, closed="left",boundary = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))  + xlab('') + xlim(0:1) +
        scale_fill_manual(values = pal, name = "Comparison") +
        # scale_fill_viridis(discrete=TRUE,option="B") +
        ylab('Number of SNPs') +
        xlab('Mean change in allele frequency') +
        theme_tufte() +
        theme(legend.position = "none") +
        ggtitle("B. Selected SNPs") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0), limits = c(0, 0.7)) +
        theme(axis.line = element_line(color = 'black'))

#All SNPs
tmp <- abs(freqDiff_200917[, 9:11])
colnames(tmp) <- c('G1 → G11', 'G1 → G25', 'G1 → G100')
dat <- pivot_longer(tmp, c('G1 → G11', 'G1 → G25', 'G1 → G100'))
dat$name <- factor(dat$name, levels = c('G1 → G11', 'G1 → G25', 'G1 → G100'))

#plot distribution
p1 <- ggplot(dat, aes(x=value, fill=name)) + geom_histogram(alpha=1, closed="left",boundary = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))  + xlab('') + xlim(0:1) +
        scale_fill_manual(values = pal, name = "Comparison") +
        # scale_fill_viridis(discrete=TRUE,option="B") +
        ylab('Number of SNPs') +
        xlab('Mean change in allele frequency') +
        theme_tufte() +
        xlim(0, 0.75) +
        theme(legend.position = c(0.7,.8)) +
        ggtitle('A. All SNPs') +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0), limits = c(0, 0.7)) +
        theme(axis.line = element_line(color = 'black'))


library(patchwork)
panel = (p1/p2)
library(cowplot)
save_plot("../figures/delta_af.png", panel, base_height = 6)
