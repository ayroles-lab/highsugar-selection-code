library(ggplot2)
library(ggpmisc)
#iHS
load('data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
# load('200925_ihsScan_HS_SNPmissing09_indMissing09.RData')
load('archive/R/200925_ihsScan_HS_SNPmissing08_indMissing05.RData')
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
{
dat <- data.frame(x = abs(wgscan.ihs$ihs$IHS), y = -log10(p))
my.formula <- y~x
cor.test(dat$x, dat$y)
high_sugar = ggplot(dat, aes(x=x, y=y)) +
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x = 'right', label.y = 0.7, rr.digits = 4,
               formula = my.formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = my.formula),
                  geom = 'text',
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 2), sep = "")),
                  label.x = 'right', label.y = 25, size = 3) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 45)) +
  theme(plot.title = element_text(size = 12),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        axis.line = element_line(color = 'black')) +
  xlab('abs(iHS)') + ylab('-log10(p)') + ggtitle(expression('A. ' ~ p[int] ~ ' by iHS'))

#Time p
mod <- lm(-log10(binGlmScan.reOrder$generation_p) ~ abs(wgscan.ihs$ihs$IHS))
summary(mod)
#C p
mod <- lm(-log10(binGlmScan.reOrder$emTrend.C_p) ~ abs(wgscan.ihs$ihs$IHS))
summary(mod)
dat <- data.frame(x = abs(wgscan.ihs$ihs$IHS), y = -log10(binGlmScan.reOrder$emTrend.C_p))
my.formula <- y~x
control = ggplot(dat, aes(x=x, y=y)) +
  geom_point() +
  geom_smooth(method=lm,  linetype="dashed") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               label.x = 'right', label.y = 0.7, rr.digits = 4,
               formula = my.formula, parse = TRUE, size = 3) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = my.formula),
                  geom = 'text',
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 2), sep = "")),
                  label.x = 'right', label.y = 25, size = 3) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 45)) +
  theme(plot.title = element_text(size = 12),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        axis.line = element_line(color = 'black'),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8)) +
  xlab('abs(iHS)') + ylab('-log10(p)') + ggtitle(expression('B. ' ~ p[C] ~ ' by iHS'))

#HS p
# mod <- lm(-log10(binGlmScan.reOrder$emTrend.HS_p) ~ abs(wgscan.ihs$ihs$IHS))
# summary(mod) #Huh
# cor.test(-log10(binGlmScan.reOrder$emTrend.HS_p), abs(wgscan.ihs$ihs$IHS))

# dat <- data.frame(x = abs(wgscan.ihs$ihs$IHS), y = -log10(binGlmScan.reOrder$emTrend.HS_p))
# my.formula <- y~x
# high_sugar = ggplot(dat, aes(x=x, y=y)) +
#   geom_point() +
#   geom_smooth(method=lm,  linetype="dashed") +
#   stat_poly_eq(aes(label = paste(..rr.label..)),
#                label.x = 'right', label.y = 0.7,
#                formula = my.formula, parse = TRUE, size = 10) +
#   stat_fit_glance(method = 'lm',
#                   method.args = list(formula = my.formula),
#                   geom = 'text',
#                   aes(label = paste("p-value = ", signif(..p.value.., digits = 2), sep = "")),
#                   label.x = 'right', label.y = 25, size = 10) +
#   theme_tufte() + 
#   theme(plot.title = element_text(hjust = 0.5, size = 25),
#         legend.title=element_blank(),
#         legend.text=element_text(size=30),
#         axis.text=element_text(size=25),
#         axis.title=element_text(size=28,face="bold"),
#         axis.line = element_line(color = 'black')) +
#   xlab('abs(iHS)') + ylab('-log10(p)') + ggtitle("A. Selection signal in HS")
panel = high_sugar + control
save_plot("figures/ihs_vs_pInt_exclConly_SNPmissing08_indMissing05_HScontrasts_Ccontrast.png", panel, base_height = 3.5, ncol = 2, base_asp = 1)
}
###### Filter significant SNPs, check if their iHS is significant

dat <- data.frame(ihs = abs(wgscan.ihs$ihs$IHS), 
                  ihs_logpvalue = wgscan.ihs$ihs$LOGPVALUE , 
                  regression_logpvalue = -log10(p)) %>% na.omit() %>%   
                  filter(regression_logpvalue > -log10(8*10^-12)) %>%
                  mutate(ihs_sig = ihs_logpvalue > -log10(0.05)) %>%
                  group_by(ihs_sig) %>% count()
dat
