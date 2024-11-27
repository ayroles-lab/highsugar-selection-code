library(here)
setwd(here::here("archive/R"))
load('200120_binGlmScan_GRangesFormat.RData')
library(GenomicRanges)

binGlmScan_contTime_param1234_190429.gr_p1e13 <- binGlmScan_contTime_param1234_190429.gr[binGlmScan_contTime_param1234_190429.gr$`generation:treatment_p` < 1e-13]

load(here('DE_tablea/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData'))
load(here('DE_tablea/200618_binGlmScan_slopePerTreat_GRangesFormat.RData'))
source(here("code/utils.R"))

window <- 200

signCut <- 8e-12
#signCut <- 10^-9
nonSignCut <- 1e-4
p <- binGlmScan_contTime_param1234_200615.gr$`generation:treatment_p`
names(p) <- binGlmScan_contTime_param1234_200615.gr@ranges@NAMES
p[binGlmScan_contTime_param1234_200615.gr$emTrend.HS_p > nonSignCut] <- 1

binGlmScan_contTime_param1234_200615.gr[p < signCut] 

seqlevelsStyle(binGlmScan_contTime_param1234_200615.gr) <- "UCSC"

options(stringsAsFactors = F)
swapExp_DESeqRes_body <- read.table(here('data/Results_DESeq_SWAPexperiment_body_May12.txt'), header = T)
swapExp_DESeqRes_head <- read.table(here('data/Results_DESeq_SWAPexperiment_head_May12.txt'), header = T)

#body
DE.fdr <- .01

#ggplot facet version
DE_table <- data.frame(tissue = rep(c('Body', "Head"), 
                               c(nrow(swapExp_DESeqRes_body), 
                                 nrow(swapExp_DESeqRes_head))),
                  gene_id = c(rownames(swapExp_DESeqRes_body), 
                              rownames(swapExp_DESeqRes_head)),
                  lfc_geno = c(swapExp_DESeqRes_body$lfc_geno, 
                               swapExp_DESeqRes_head$lfc_geno),
                  p_geno = c(swapExp_DESeqRes_body$padj_geno,                
                             swapExp_DESeqRes_head$padj_geno))
DE_table$sign <- 0
DE_table$sign[DE_table$p_geno < DE.fdr] <- 1

pak::pkg_install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
gene_locations_GR <- transcripts(txdb, columns = "gene_id")
gene_locations_GR$gene_id = unlist(gene_locations_GR$gene_id)
chrs <- c("2L", "2R", "3L", "3R", "4", "X")

gene_locations_GR <- gene_locations_GR[gene_locations_GR@seqnames %in% paste0("chr", chrs)]

seqlevelsStyle(gene_locations_GR)

tissue = "Body"
signCut <- 8e-12
signCut <- 1e-5
DE.fdr <- 1e-3
runDEoverlap = function(signCut, tissue, DE.fdr){
    DE_table_t = DE_table[DE_table$tissue == tissue,]
    DE_gene_ids = DE_table_t[DE_table_t$p_geno < DE.fdr, "gene_id"]
    print(length(DE_gene_ids))
    DE_gene_locations = gene_locations_GR[gene_locations_GR$gene_id %in% DE_gene_ids]

    selected_snps = binGlmScan_contTime_param1234_200615.gr[p < signCut] 
    #selected_snps <- selected_snps[selected_snps@seqnames %in% chrs]
    n_snps = length(selected_snps)

    gene_location = DE_gene_locations + 200

    tiles <- tileGenome(seqinfo(txdb), ntile=1000000)
    gene_tiles = findOverlaps(tiles, gene_location)
    is_gene = rep(0, length(tiles))
    is_gene[unique(gene_tiles@from)] = 1

    overlap = findOverlaps(selected_snps, tiles[is_gene == 1])
    n_overlap_snps = overlap@from |> unique() |> length()
    observed_overlap = n_overlap_snps/n_snps * 100

    permute_tiles = function(){
        perm_is_gene = sample(is_gene)
        overlap = findOverlaps(selected_snps, tiles[perm_is_gene == 1])
        n_overlap_snps = overlap@from |> unique() |> length()
        permuted_overlap = n_overlap_snps/n_snps * 100
        return(permuted_overlap)
    }
    null = replicate(1000, permute_tiles())
    return(list(obs = observed_overlap, null = null))
}
x = runDEoverlap(1e-5, "Head")
x$obs
summary(x$null)

pval_thresh <- c(1e-5, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30)

library(purrr)
library(furrr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
head_DE_overlap_perm = map(pval_thresh, runDEoverlap, "Head", 1e-2, .progress = TRUE)
body_DE_overlap_perm = map(pval_thresh, runDEoverlap, "Body", 1e-2, .progress = TRUE)

obs_df = bind_rows("Head" = data.frame(variable = pval_thresh, value = map_dbl(head_DE_overlap_perm, pluck, "obs")), 
                   "Body" = data.frame(variable = pval_thresh, value = map_dbl(body_DE_overlap_perm, pluck, "obs")), 
                   .id = "tissue")

null_df_head = map(head_DE_overlap_perm, pluck, "null") |> as.data.frame()
colnames(null_df_head) = pval_thresh
null_df_head = reshape2::melt(null_df_head)
head(null_df_head)

null_df_body = map(body_DE_overlap_perm, pluck, "null") |> as.data.frame()
colnames(null_df_body) = pval_thresh
null_df_body = reshape2::melt(null_df_body)
head(null_df_body)
null_df = bind_rows(Head = null_df_head, Body = null_df_body, .id = "tissue")

g.geno_facet <- ggplot(data = null_df, 
                       aes(y = value,
                           x = factor(-log10(as.numeric(as.character(variable)))))) +
    geom_boxplot(key_glyph = 'rect', linewidth = 0.3, outlier.size = 0.1) +
    geom_point(data = obs_df,
            mapping = aes(y = value, 
                          x = factor(-log10(as.numeric(as.character(variable))))), 
            group = 1, size = 1, color = 2) +
    geom_line(data = obs_df,
              aes(y = value, 
                  x = factor(-log10(as.numeric(as.character(variable))))), 
            group = 1, linewidth = 0.7, color = 2) + 
            theme_classic() +
    facet_grid(cols = vars(tissue)) +
   theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) + 
  ylab('Percent of selected SNPs in DE genes (%)') + xlab(expression('-log10(' ~ p[int] ~ ')')) 
  
save_plot(here::here("test.png"), g.geno_facet)  
  
  facet_grid(cols = vars(tissue)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)) + ylim(c(0.06, 0.22)) +
  ylab('Fraction of selected SNPs in DE genes') + xlab(expression('-log10(' ~ p[int] ~ ')')) 

### GO enrichment DE genes

library(clusterProfiler)
library(enrichplot)
library(org.Dm.eg.db)

perform_enrichment <- function(genes, background) {
  genes = na.omit(genes)
  go_enrichment_result <- enrichGO(
    gene = genes,
    universe = background,
    OrgDb = org.Dm.eg.db,
    keyType = "FLYBASE",
    ont = "BP",  # Biological Process ontology
    pAdjustMethod = "BH",
    qvalueCutoff = 0.01
  )
 return(go_enrichment_result)
}

tissue = "Body"
DE_table_t = DE_table[DE_table$tissue == tissue,]
DE_gene_ids = DE_table_t[DE_table_t$p_geno < 0.01, "gene_id"]
print(length(DE_gene_ids))
GO_body = perform_enrichment(DE_gene_ids, DE_table_t$gene_id)

tissue = "Head"
DE_table_t = DE_table[DE_table$tissue == tissue,]
DE_gene_ids = DE_table_t[DE_table_t$p_geno < 0.01, "gene_id"]
print(length(DE_gene_ids))
GO_head = perform_enrichment(DE_gene_ids, DE_table_t$gene_id)

library(rio)
export(GO_body@result, here("GO_DEgenes_body.csv"))
export(GO_head@result, here("GO_DEgenes_head.csv"))

edo <- pairwise_termsim(GO_head)
p1 <- emapplot(edo, cex_category=1.5,layout="kk") + ggtitle("Head")

edo <- pairwise_termsim(GO_body)
p2 <- emapplot(edo, cex_category=1.5,layout="kk") + ggtitle("Body")
save_plot(here("DEtest-Head.png"), p1, base_width = 15, base_height=10)
save_plot(here("DEtest-Body.png"), p2, base_width = 15, base_height=10)
