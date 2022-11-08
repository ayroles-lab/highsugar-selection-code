load('../data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
load('../data/200618_binGlmScan_slopePerTreat_GRangesFormat.RData')
signCut <- 8e-12
nonSignCut <- 1e-4
p <- binGlmScan_contTime_param1234_200615.gr$`generation:treatment_p`
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

