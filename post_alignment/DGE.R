#BiocManager::install('edgeR')
library(edgeR)
library(RColorBrewer)
library(rtracklayer)
library(dplyr)
rm(list = ls())
setwd('C:/Users/zeno/Dropbox/Hemp/post_alignment/')
# ------------------- 
DGE_norm <- readRDS('../2024Jun-RNASeq/2024Jun-RNASeq/Gene counts/rRNA-filtered_and_normalized_DGEList.rds')
mdt <- DGE_norm$samples

files <- paste0('../2024Jun-RNASeq/2024Jun-RNASeq/Gene counts/raw gene counts/',
                 dir('../2024Jun-RNASeq/2024Jun-RNASeq/Gene counts/raw gene counts/'))
DGE_raw <- readDGE(files, sep = ' ')

my_obj <- import("../2024Jun-RNASeq/2024Jun-RNASeq/GCF_029168945.1_ASM2916894v1_genomic.gtf.gz")
class(my_obj)
dt_GTF <- data.table::as.data.table(my_obj@elementMetadata@listData)

# ----------- PCA plot ------------ #
lcpm <- cpm(DGE_norm, log=TRUE)
par(mfrow=c(1,2))
group <- DGE_norm$samples$group
cultivar <- DGE_norm$samples$cultivar
col.group <- as.factor(cultivar)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=cultivar, col=col.group)

sampling_stage <- DGE_norm$samples$growth.stage
col.group <- as.factor(sampling_stage)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=sampling_stage, col=col.group)

#plotMDS(lcpm, labels=cultivar, col=col.group)

# d <- dist(lcpm)
# mds <- as.data.frame(as.matrix(cmdscale(d)))
# library(ggrepel)
# dt_mds_plot <- data.frame(Sample=rownames(mds), mds)
# ggplot(df, aes(x=V1, y=V2, label=Sample)) +
#   geom_point() +
#   geom_label_repel(min.segment.length = 0)


# ------------ filter out extremely skewed data ------------- #

cutoff <- 1
drop <- which(apply(cpm(DGE_norm), 1, max) < cutoff)
DGE_norm_filtered <- DGE_norm[-drop,] 
dim(DGE_norm_filtered) # number of genes left

# ------------ contrast individual groups ---------- #
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)

par(mfrow=c(1,1))
v <- voom(DGE_norm, design, plot=TRUE)
vfit <- lmFit(v, design)
colnames(coef(vfit))

v_filter <- voom(DGE_norm_filtered, design, plot=TRUE)
vfit_filtered <- lmFit(v_filter, design)
colnames(coef(vfit_filtered))

# -------- test contrasts ------- #
design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }

contr.matrix <- design.pairs(colnames(coef(vfit_filtered)))
efit <- eBayes(contrasts.fit(vfit_filtered, contrasts=contr.matrix))
plotSA(efit, main="Final model: Mean-variance trend")

tfit <- treat(contrasts.fit(vfit_filtered, contrasts=contr.matrix), lfc=1)

comb_treats <- colnames(contr.matrix)

ls.top.treat <- list()
for (i in comb_treats) {
  ls.top.treat[[i]] <- topTreat(tfit, coef=i, n=Inf)
}

ls.top.table <- list()
for (i in comb_treats) {
  ls.top.table[[i]] <- topTable(efit, coef=1, sort.by = "P", n = Inf)
}

head(ls.top.table[['EAV-EAV_F']], 5)
sum(ls.top.table[['EAV-EAV_F']]$adj.P.Val < 0.05)
head(ls.top.treat[['EAV-EAV_F']], 5)
sum(ls.top.treat[['EAV-EAV_F']]$logFC > 2)

# ----------- select a pair to find DE gene candidates -------------- #
dt <- decideTests(tfit)
summary(dt)
comb_treats
sel_pair_1 <- which(comb_treats == 'EAV-EAV_F')
sel_pair_2 <- which(comb_treats == 'SG-SG_F')

vennDiagram(dt[,c(sel_pair_1, sel_pair_2)],
            circle.col=c("turquoise", "salmon"))

de.common <- which(dt[,sel_pair_1]!=0 & dt[,sel_pair_2]!=0)
length(de.common)

de.unique <- which((dt[,sel_pair_1]!=0 & dt[,sel_pair_2]==0) | (dt[,sel_pair_1]==0 & dt[,sel_pair_2]!=0))
length(de.unique)

dimnames(dt)[[1]][de.unique][1:10]
# ----- join annotation data ----- #

dt_gene_cand <- data.table::data.table(gene_id = dimnames(dt)[[1]][de.unique])
dt_gene_cand <- dplyr::left_join(dt_gene_cand, dt_GTF, by = 'gene_id')
dt_gene_cand <- filter(dt_gene_cand, gbkey == 'Gene')


# ------- heatmap -------- #
library(gplots)
mycol <- colorpanel(1000,"blue","white","red")
dt_plot <- lcpm[dt_gene_cand$gene_id, ]
colnames(dt_plot) <- paste0(group, paste0('_',seq(1,3)))
dt_plot_f <- dt_plot[, grepl('EAV|SG', colnames(dt_plot))]

heatmap.2(dt_plot_f, scale="row",
          labRow=dt_gene_cand$description, labCol=colnames(dt_plot_f),
          col=mycol, trace="none", density.info="none",
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

png(paste0('DGE_hm_', comb_treats[sel_pair_1], '_VS_',
           comb_treats[sel_pair_2], '.png'), plot, width = 14, height = 8, units = 'in', res = 200)

heatmap.2(dt_plot_f, scale="row",
          labRow=dt_gene_cand$description, labCol=colnames(dt_plot_f),
          col=mycol, trace="none", density.info="none",
          margin=c(12,32), 
          lhei=c(2,10), 
          cexRow=1.2,
          cexCol=1.5,
          adjRow = c(0,NA),
          #srtRow=45,
          dendrogram="column")

dev.off()
