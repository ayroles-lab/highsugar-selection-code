#col scales for a draft MS fig
load('archive/R/210119_pvals_chisqANDcor_HS_matForm.RData')

tmp <- round(-log10(corP_chisqP_HS.mat))
mat1 <- mat2 <- tmp
#mat1[mat1 > 70] <- 70
mat1[lower.tri(mat1)] <- NA
mat2[upper.tri(mat2)] <- NA
col <- colorRampPalette(c('white', '#ad0303'))
dim(mat1)
library(corrplot)
splits = table(strsplit(colnames(mat1), "_") |> sapply(\(x) x[[2]]))
splits = cumsum(c(1, splits))
png('figures/20230216_corrplotEpistaticChisqCorr.png', width = 1200, height = 1200)
par(family = "Serif")
corrplot(mat1, diag = F, order = "original", col = col(200),
         tl.pos = "n", tl.cex = 0.5, method = "color", type = "upper", is.corr = F,
         cl.cex = 2, mar = c(4,5,1,4)) |>
corrRect(index = splits, lwd = 0.8)
corrplot(mat2, diag = F, order = "original", col = col(200),
         tl.pos = "n", tl.cex = 0.5, method = "color", type = "lower", is.corr = F,
         cl.cex = 2, add = T, mar = c(4,5,1,4)) |>
corrRect(index = splits, lwd = 0.8)
dev.off()


#A manhattan plot with the same SNPs
#Ideally, I'd like to make the SNPs line up exactly in the heatmap + manhattan on the top/side. That seems like a pain though
load('archive/R/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData')
snps <- colnames(corP_chisqP_HS.mat)
selection.snps <- paste(binGlmScan_contTime_param1234_200615$POS, binGlmScan_contTime_param1234_200615$CHROM, sep = '_')
selection <- binGlmScan_contTime_param1234_200615[match(snps, selection.snps), ]

nonSignCut <- 1e-4
signCut <- 8e-12
p <- selection$`generation:treatment_p`
p[selection$emTrend.HS_p > nonSignCut] <- 1
png('test.png', width = 1200, height = 600)
par(mar = c(5,5,4,2))
plot.manhattan2(pos = selection$POS,
                chr = selection$CHROM,
                y = p, col = wes_palette('Darjeeling2')[2], cex = .8)
mtext(expression('-log10(' ~ p[int] ~ ')'), side = 2, cex = 2, line = 2.7)
dev.off()

pos <- list()
chrs <- unique(selection$CHROM)
for (i in 1:length(chrs)) {
  pos[[i]] <- order(selection$POS[selection$CHROM == chrs[i]])
}
pos <- unlist(pos)

png('test.png', width = 1200, height = 600)
par(mar = c(5,5,4,2))
plot.manhattan2(pos = pos,
                chr = selection$CHROM,
                y = p, col = wes_palette('Darjeeling2')[2], cex = .8)
mtext(expression('-log10(' ~ p[int] ~ ')'), side = 2, cex = 2, line = 2.7)
dev.off()


########## 210120 #########
#Can I make the heatmap "spaced" by physical positions?
pos <- as.numeric(sub(pattern = '(.*)_(.*)', replacement = '\\1', rownames(corP_chisqP_HS.mat)))
chr <- sub(pattern = '(.*)_(.*)', '\\2', rownames(corP_chisqP_HS.mat))
tmp <- sortmap(chr, pos)
pos.cum <- tmp$cummap

# mat <- apply(corP_chisqP_HS.mat, 2, rev)
# pos <- as.numeric(sub(pattern = '(.*)_(.*)', replacement = '\\1', rownames(mat)))
# chr <- sub(pattern = '(.*)_(.*)', '\\2', rownames(mat))
# tmp <- sortmap(chr, pos)
# pos.cum <- tmp$cummap

upper <- lower <- -log10(corP_chisqP_HS.mat)
upper[upper > 70] <- 70
upper[lower.tri(upper)] <- NA
lower[upper.tri(lower)] <- NA
png('test.png', 1200, 1200)
par(mfrow = c(2,1))
#manhattan
plot(pos.cum, -log10(p), col = wes_palette('Darjeeling2')[2], cex = .8, pch = 19)
# mtext(expression('-log10(' ~ p[int] ~ ')'), side = 2, cex = 2, line = 2.7)

#heatmap
# image(x = rev(pos.cum), y = rev(pos.cum), z = -log10(t(mat)))
# image(z = -log10(t(mat)))
image(x = pos.cum, y = pos.cum, z = upper)
image(x = pos.cum, y = pos.cum, z = lower, add = T)
dev.off()

#Trying a ggplot approach
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
#Manhattan
nonSignCut <- 1e-4
signCut <- 8e-12
p <- selection$`generation:treatment_p`
p[selection$emTrend.HS_p > nonSignCut] <- 1
df <- data.frame(pos = pos.cum, p = -log10(p))
p1 <- ggplot(df, aes(x = pos, y = p)) +
  geom_point(pch = 19) + theme_bw() +
  # Change the expand values
  scale_x_discrete(expand = c(0.05, 0.05)) +
  #scale_y_discrete(breaks = NULL) +
  theme(legend.position = "none") + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())

#Heatmap
mat <- -log10(corP_chisqP_HS.mat)
rownames(mat) <- colnames(mat) <- pos.cum
# mat.melt <- melt(mat)
mat.melt <- melt(mat[2:5, 2:5])

w <- diff(mat.melt$Var1)
w[w < 0] <- sort(unique(w))[2]
w <- c(w, min(w))

p2 <- ggplot(mat.melt,
             aes(x = Var1, y = Var2, fill = value)) +
  # geom_bar(stat="identity")
  geom_tile(aes(width = w, height = 1e4)) +

  scale_fill_gradientn(colours = c("blue", "white" , "#ad0303")) +
  theme(legend.position = "none") + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) +
  scale_x_discrete(expand = c(0, 0), breaks = letters[1:10]) +
  scale_y_discrete(expand = c(0, 0))

png('test.png', 1200, 1200)
p2
dev.off()



#Trying to create scaled positions
pos <- as.numeric(sub(pattern = '(.*)_(.*)', replacement = '\\1', rownames(corP_chisqP_HS.mat)))
chr <- sub(pattern = '(.*)_(.*)', '\\2', rownames(corP_chisqP_HS.mat))
nonSignCut <- 1e-4
signCut <- 8e-12
p <- selection$`generation:treatment_p`
p[selection$emTrend.HS_p > nonSignCut] <- 1

tmp <- sortmap(chr, pos, delta = 5e6)
pos.cum <- tmp$cummap
pos.dist <- diff(pos.cum)
#This determines the amount of "scaling".
#   Large => closer to true physical pos (good manhattan, bad heatmap)
#   Small => further from physical pos (bad manhattan, good heatmap)
breaks <- 50
pos.dist_cuts <- cut(pos.dist, breaks = breaks, labels = F)
# pos.cum_scaled <- c(1, pos.dist_cuts)
pos.cum_scaled <- rep(NA, length(pos.cum))
pos.cum_scaled[1] <- 1
for (i in 1:length(pos.dist)) {
  pos.cum_scaled[i+1] <- pos.cum_scaled[i] + pos.dist_cuts[i]
}
#Get midpoints per chr for x-axis
chrMids <- 1:4
for(i in 1:length(unique(chr)))
  chrMids[i] <- mean(pos.cum_scaled[chr == unique(chr)[i]])
#Put in NA entries at chr boundaries. This is a workaround to draw them in the heatmap
bounds <- cumsum(table(chr))
pos.cum_scaled <- c(pos.cum_scaled[1:bounds[1]],
                    pos.cum_scaled[bounds[1]] + breaks/2,
                    pos.cum_scaled[(bounds[1] + 1):bounds[2]],
                    pos.cum_scaled[bounds[2]] + breaks/2,
                    pos.cum_scaled[(bounds[2] + 1):bounds[3]],
                    pos.cum_scaled[bounds[3]] + breaks/2,
                    pos.cum_scaled[(bounds[3] + 1):bounds[4]])
p <- c(p[1:bounds[1]], NA,
       p[(bounds[1] + 1):bounds[2]], NA,
       p[(bounds[2] + 1):bounds[3]], NA,
       p[(bounds[3] + 1):bounds[4]])

mat <- matrix(NA, length(pos.cum_scaled), length(pos.cum_scaled)) #New larger matrix, with extra cols/rows corresponding to chr bounds
mat[-which(is.na(p)), -which(is.na(p))] <- -log10(corP_chisqP_HS.mat) #Insert p-vals everwhere but at the bounds
diag(mat) <- 0
rownames(mat) <- colnames(mat) <- pos.cum_scaled
#test
# plot(pos.cum, -log10(p), col = wes_palette('Darjeeling2')[2], cex = .8, pch = 19)
# plot(pos.cum_scaled, -log10(p), col = wes_palette('Darjeeling2')[2], cex = .8, pch = 19)

#format data for ggplot
df.manhattan <- data.frame(pos = pos.cum_scaled, p = -log10(p))
upper <- lower <- mat
upper[upper > 70] <- 70
upper[lower.tri(upper) & !is.na(upper)] <- 0
lower[upper.tri(lower) & !is.na(lower)] <- 0
upper.melt <- melt(upper)
lower.melt <- melt(lower)
df.heatmap <- data.frame(upper.melt[,1:2], upper = upper.melt$value, lower = lower.melt$value)

#test
# idx <- 1:20
# df.manhattan <- df.manhattan[idx, ]
# upper.melt <- melt(upper[idx,idx])
# lower.melt <- melt(lower[idx,idx])
# df.heatmap <- data.frame(upper.melt[,1:2], upper = upper.melt$value, lower = lower.melt$value)

bounds.idx <- which(is.na(df.manhattan$p))
shadeCol <- 'gray80'
p1 <- ggplot(df.manhattan, aes(x = pos, y = p)) +
  geom_rect(aes(xmin=0, xmax=pos[bounds.idx[1]], ymin=-Inf, ymax=Inf), fill = shadeCol) +
  geom_rect(aes(xmin=pos[bounds.idx[2]], xmax=pos[bounds.idx[3]], ymin=-Inf, ymax=Inf), fill = shadeCol) +
  geom_point(pch = 19) + theme_bw() +
  # Change the expand values
  scale_x_continuous(breaks = chrMids, labels = c('2L', '2R', '3L', '3R')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete() +
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")


p.heat_chisq <- ggplot(dplyr::filter(df.heatmap, upper > 0)) +  theme_bw() +
  # geom_bar(stat="identity")
  geom_tile(aes(x = Var1, y = Var2, fill = upper)) +
  # geom_tile(aes(x = Var1, y = Var2, fill = lower)) +
  scale_fill_gradientn(colours = c("white" , "#ad0303"), na.value = 'black') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = df.manhattan$pos[bounds.idx]) +
  geom_hline(yintercept = df.manhattan$pos[bounds.idx]) + 
  theme_minimal() +   theme(legend.position = "none") + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) 
p.heat_cor <- ggplot(dplyr::filter(df.heatmap, lower > 0)) +  theme_bw() +
  # geom_bar(stat="identity")
  # geom_tile(aes(x = Var1, y = Var2, fill = upper)) +
  geom_tile(aes(x = Var1, y = Var2, fill = lower)) +
  scale_fill_gradientn(colours = c("white" , "#ad0303"), na.value = 'black') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = df.manhattan$pos[bounds.idx]) +
  geom_hline(yintercept = df.manhattan$pos[bounds.idx]) + 
  theme_minimal() + theme(legend.position = "none") + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())

p_both <- arrangeGrob(p1, p.heat_chisq, heights = c(2,8))
png('testX2.png', 1200, 1200/.8, bg = "transparent")
plot(p_both)
dev.off()
p_both <- arrangeGrob(p1, p.heat_cor, heights = c(2,8))
png('testCor.png', 1200, 1200/.8)
plot(p_both)
dev.off()

source('https://gist.githubusercontent.com/eliocamp/eabafab2825779b88905954d84c82b32/raw/8af496a49e72361e305b580c8701e2337b37fd77/new_aes.R')

library(patchwork)
chr_boxes = data.frame(start = c(1, df.manhattan$pos[bounds.idx]),
                         end = c(df.manhattan$pos[bounds.idx], df.heatmap$Var1 |> max())) 
p.heat_chisqCor <- ggplot() + 
  # geom_bar(stat="identity")
  geom_tile(data= dplyr::filter(df.heatmap, upper > 0), aes(x = Var1, y = Var2, fill = upper)) +
  scale_fill_gradientn(colours = c("white" , "#ad0303"), na.value = 'black') +
  new_scale("fill") +
  geom_tile(data= dplyr::filter(df.heatmap, lower > 0), aes(x = Var1, y = Var2, fill = lower)) +
  # geom_tile(aes(x = Var1, y = Var2, fill = lower)) +
  scale_fill_gradientn(colours = c("white" , "#ad0303"), na.value = 'black') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_rect(data = chr_boxes, color = "black", alpha = 0, linewidth = 1,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end)) +
  #geom_vline(xintercept = df.manhattan$pos[bounds.idx]) +
  #geom_hline(yintercept = df.manhattan$pos[bounds.idx]) + 
  theme_bw() +   theme(legend.position = "none") + 
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) 
layout = 
"A
B
B
B
B"
p_both <- arrangeGrob(p1, p.heat_chisqCor, heights = c(2,8))
png('test.png', 1200, 1200/.8)
plot(p_both)
dev.off()
