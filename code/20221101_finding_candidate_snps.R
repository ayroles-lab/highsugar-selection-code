load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
library(tidyverse)

chrs <- c('2L', '2R', '3L', '3R', '4', 'X')


signCut <- 1e-09
nonSignCut <- 1e-4

# Sig interaction and sig in HS
binGlmScan_contTime_param1234_200615 %>% filter(`generation:treatment_p` < signCut, CHROM %in% chrs, emTrend.HS_p < nonSignCut) %>% count


signCut <- 8e-12
nonSignCut <- 1e-4

candidate_snps <- binGlmScan_contTime_param1234_200615 %>% filter(`generation:treatment_p` < signCut, CHROM %in% chrs, emTrend.HS_p < nonSignCut)
candidate_snps %>% count

candidate_snps %>% group_by(CHROM) %>% count
