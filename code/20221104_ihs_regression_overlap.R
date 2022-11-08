#iHS
load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
# load('200925_ihsScan_HS_SNPmissing09_indMissing09.RData')
load('../archive/R/200925_ihsScan_HS_SNPmissing08_indMissing05.RData')
chrs <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
binGlmScan_contTime_param1234_200615.chrs <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in% chrs, ]
#HS
#Match binGlm and IHS
binGlmScan.names <- paste(binGlmScan_contTime_param1234_200615.chrs$CHROM, binGlmScan_contTime_param1234_200615.chrs$POS, sep = '_')
wgscan.ihs.names <- paste(wgscan.ihs$ihs$CHR, wgscan.ihs$ihs$POSITION, sep = '_')
binGlmScan.reOrder <- binGlmScan_contTime_param1234_200615.chrs[match(wgscan.ihs.names, binGlmScan.names), ]
sum(binGlmScan.reOrder$POS - wgscan.ihs$ihs$POSITION, na.rm = T) #Check
sum(wgscan.ihs.names %in% binGlmScan.names) #Not all ihs SNPs are a subset of the binGlm SNPs. That's because a the binglm ones are the intersect with the ones that where typed in g1-25
sum(binGlmScan.reOrder$CHROM == wgscan.ihs$ihs$CHR, na.rm = T) #Check

mod <- lm(-log10(binGlmScan.reOrder$`generation:treatment_p`) ~ abs(wgscan.ihs$ihs$IHS))
summary(mod)

nonSignCut <- 1e-4
p <- binGlmScan.reOrder$`generation:treatment_p`
p[binGlmScan.reOrder$emTrend.HS_p > nonSignCut] <- 1

#Prettier version
library(ggplot2)
library(ggpmisc)
dat <- data.frame(x = abs(wgscan.ihs$ihs$IHS), y = -log10(p))
my.formula <- y~x
png('../figures/ihs_vs_pInt_exclConly_SNPmissing08_indMissing05_v2.png', 1000, 1000)
ggplot(dat, aes(x=x, y=y)) +
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed") + theme_classic() +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x = .86, label.y = 0.35,
               formula = my.formula, parse = TRUE, size = 15) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = my.formula),
                  geom = 'text',
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 2), sep = "")),
                  label.x = 'right', label.y = 10, size = 15) +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.title=element_blank(),
        legend.text=element_text(size=30),
        axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold")) +
  xlab('abs(iHS)') + ylab('-log10(p)')
dev.off()

###### Filter significant SNPs, check if their iHS is significant

dat <- data.frame(ihs = abs(wgscan.ihs$ihs$IHS), 
                  ihs_logpvalue = wgscan.ihs$ihs$LOGPVALUE , 
                  regression_logpvalue = -log10(p)) %>% na.omit() %>%   
                  filter(regression_logpvalue > -log10(8*10^-12)) %>%
                  mutate(ihs_sig = ihs_logpvalue > -log10(0.01)) %>%
                  group_by(ihs_sig) %>% count()
