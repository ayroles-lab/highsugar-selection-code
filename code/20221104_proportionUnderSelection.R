load('data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
load('data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')
source("code/utils.R")

signCut <- 8e-12
#signCut <- 10^-9
nonSignCut <- 1e-4
p <- binGlmScan_contTime_param1234_200615.gr$`generation:treatment_p`
names(p) <- binGlmScan_contTime_param1234_200615.gr@ranges@NAMES
p[binGlmScan_contTime_param1234_200615.gr$emTrend.HS_p > nonSignCut] <- 1
sum(p < signCut)
sum(p < signCut & binGlmScan_contTime_param1234_200615$CHROM %in% chrs)

width <- 200 #Look this far around every SNP
merge <- 1 #Merge regions closer than this
binGlmScan_contTime_param1234_200615.gr_keep <- binGlmScan_contTime_param1234_200615.gr[p < signCut]
binGlmScan_contTime_param1234_200615.gr_keep <- binGlmScan_contTime_param1234_200615.gr_keep + width/2
binGlmScan_contTime_param1234_200615.gr_keep_regions <- reduce(binGlmScan_contTime_param1234_200615.gr_keep, min.gapwidth = merge)
ranges <- binGlmScan_contTime_param1234_200615.gr_keep_regions@ranges
sum(ranges@width)/1e6
sum(ranges@width)/140e6


load('data/190916_maf05SNPs_alleleFreqs.RData')
chrs <- c("2L", "2R", "3L", "3R", "4", "X")

keep <- p < signCut & binGlmScan_contTime_param1234_200615$CHROM %in% chrs
selected.snps <- binGlmScan_contTime_param1234_200615.gr@ranges@NAMES[keep]

colnames(allAbove.freqs) %in% keep.snps

library(tidyverse) 
library(rethinking)

G1_HS <- t(allAbove.freqs[paste0("G1_", c("N4", "N5", "N6")),])
dat = tibble(id = rownames(G1_HS), HS1 = G1_HS[,1], HS2 = G1_HS[,2], HS3 = G1_HS[,3]) |>
    mutate(selected = id %in% selected.snps, pval = p[id]) |>
    filter(pval != 1)



names(dat)

png("HS-selected-SFS.png")
par(mfrow=c(1, 1))
hist(1-dat$HS1[dat$selected], breaks = 20, main = "Observed Allele frequency distributio\nSelected sites", xlab = "Allele Frequency", prob = TRUE)
dev.off()

png("figures/af_selected_plot.png")
ggplot(pivot_longer(dat, HS1:HS3, names_to = "pop", values_to = "af"), aes(selected, af, group = interaction(pop, selected))) + geom_boxplot() + facet_wrap(~pop) + labs(y = "Allele frequency at Generation 1", x = "Selected") +  theme_half_open(16) +   panel_border() +
  background_grid()
dev.off()

dat$selected = as.numeric(dat$selected)
m1 = glm(selected ~ HS1 + HS2 + HS3 - 1, data = dat, family = "binomial") 
summary(m1)
plot(precis(m1, pars = c("HS1", "HS2", "HS3")))
dev.off()

m2 = lm(log(pval) ~ HS1 + HS2 + HS3 - 1, data = dat) 
summary(m2)
plot(precis(m2, pars = c("HS1", "HS2", "HS3")))
dev.off()

logit

png("figures/af_pval_lm.png")
ggplot(pivot_longer(dat, HS1:HS3, names_to = "pop", values_to = "af"), aes(logit(af), -log10(pval))) + geom_point() + geom_smooth() + facet_wrap(~pop)
dev.off()


library(devtools)
install_github("rmcelreath/rethinking")
library(rethinking)

