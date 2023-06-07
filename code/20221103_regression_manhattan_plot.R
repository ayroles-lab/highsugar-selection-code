######## New Manhatan plot 221103 ############
# Fig 2D
source("code/utils.R")
load('data/201001_alleFreqPC_allSNPs.RData')
summary(allSNPs.pc)

library(wesanderson)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(grid)

load(
    "data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData"
)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in%
    chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-04
signCut <- 8e-12
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
binGlmScan_contTime_param1234_200615$p = p

chrs <- c("2L", "2R", "3L", "3R", "4", "X")
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in%
    chrs, ]

ylim <- c(0, 60)
nonSignCut <- 1e-04
signCut <- 8e-12
p <- binGlmScan_contTime_param1234_200615$`generation:treatment_p`
p[binGlmScan_contTime_param1234_200615$emTrend.HS_p > nonSignCut] <- 1
binGlmScan_contTime_param1234_200615$p = p

nonSignCut <- 1e-04
signCut <- 8e-12

gwas_data_load = binGlmScan_contTime_param1234_200615
gwas_data_load$chr = factor(binGlmScan_contTime_param1234_200615$CHROM, levels = chrs)
gwas_data_load$bp = binGlmScan_contTime_param1234_200615$POS
sig_data <- gwas_data_load %>% 
  subset(p < signCut)
notsig_data <- gwas_data_load %>% 
  subset(p >= signCut) %>%
  group_by(chr) %>% 
  sample_frac(0.5)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

manhplot <- ggplot(filter(gwas_data, p < signCut), aes(x = bp_cum, y = -log10(p), 
                                  color = chr)) +
  geom_hline(yintercept = -log10(signCut), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, size = 0.25) +
  geom_point(data = filter(gwas_data, p > signCut), alpha = 0.75, size = 0.25, color = "grey") +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks = seq(0, 60, 10)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromossome", 
       y = expression("-log10(" ~ p[int] ~ ")")) + 
  theme_classic() +
  theme( 
    legend.position = "none",
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #axis.title.y = element_markdown(),axis.line = element_line(color = 'black'), 
    plot.title = element_text(size = 10), 
    axis.text = element_text(size = 7), 
    axis.text.x = element_text(size = 7, vjust = 0.5),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )
save_plot("figures/manhattan_pInt_HSonly_cut8e-12_ggplot.png", manhplot, base_height = 5)

######## Regression plots 221103 #################
load('data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
load('data/190427_mergedFreqsANDnrSamples_unionSNPs.RData')

chrs <- c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
binGlmScan_contTime_param1234_200615 <- binGlmScan_contTime_param1234_200615[binGlmScan_contTime_param1234_200615$CHROM %in% chrs, ]

plotFreqVsTime = function(f){
  col.treat <- c('DarkOrange', '#8A2BE2')
  lf = f %>% 
    select(contains("_N")) %>% 
    pivot_longer(G25_N5:G100_N6) %>% 
    separate(name, c("gen", "pop")) %>%
    mutate(treat = ifelse(pop %in% c('N4','N5','N6'), "HS", "Control"),
           gen = factor(gen, levels = c("G1", "G11", "G25", "G100")))
  p = ggplot(lf, aes(gen, value, color = treat, group = pop, shape = treat)) + 
    geom_point() + 
    geom_line() + 
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values = col.treat, name = "") +
    scale_shape_discrete(name = "") +
    theme_classic() + 
    theme(axis.line = element_line(color = 'black'), 
          plot.title = element_text(size = 10), 
          axis.text = element_text(size = 7), 
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    labs(y = "Allele frequency", x = "Generation")
  p
}


binGlmScan_contTime_param1234_200615 = binGlmScan_contTime_param1234_200615 %>% 
  mutate(p_int = `generation:treatment_p`,
         p_c = emTrend.C_p,
         p_HS = emTrend.HS_p,
         p_t = generation_p)

#Both
arrange(binGlmScan_contTime_param1234_200615, p_int) |> head()

binGlmScan_contTime_param1234_200615 %>%
  filter(p_int < 10e-20, p_c < 10e-20, p_HS < 10e-20) %>% arrange(p_int, p_c, p_HS) %>% head(10)
f <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == 'X' & freqs.merge_union_190426$POS == 16255696, ]

# Both same direction
binGlmScan_contTime_param1234_200615 %>%
  filter(p_t < 10e-20, p_int > 0.05) %>% arrange(p_int, p_c, p_HS) %>% head(10)
panelA <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '2R' & freqs.merge_union_190426$POS == 11258278, ]
save_plot("test.png", plotFreqVsTime(f1))

#only HS
binGlmScan_contTime_param1234_200615 %>%
  filter(p_int < 10e-20, p_c > 0.05, p_HS < 10e-20) %>% arrange(p_int, p_c, p_HS) %>% head(10)
panelC <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3R' & freqs.merge_union_190426$POS == 24851020, ]

#only C
binGlmScan_contTime_param1234_200615 %>%
  filter(p_int < 10e-20, p_HS > 0.05, p_c < 10e-20) %>% arrange(p_int, p_c, p_HS) %>% head(10)
panelB <- freqs.merge_union_190426[freqs.merge_union_190426$CHROM == '3R' & freqs.merge_union_190426$POS == 7581540, ]


p1 = plotFreqVsTime(panelA) + ggtitle("A. HS and Control") 
p2 = plotFreqVsTime(panelB) + ggtitle("B. Control only") 
p3 = plotFreqVsTime(panelC) + ggtitle("C. HS only") 

library(patchwork)
panel = (p1 + p2 + p3 + plot_layout(guides = 'collect')) / manhplot + ggtitle("D.") 

save_plot("figures/regressionExamples_manhattan_pInt_HSonly_cut8e-12.png", panel, 
          base_height = 4, base_width = 7.2)



