
module load R/3.5.0
module load hdf5/1.8.18_gcc-6.2.0  # for hdf5r required by Seurat
module load gcc/7.2.0  # for new installs
R

require(tidyverse)
require(gmodels)
require(Seurat)
require(rms)
require(scales)
require(reshape)
library(Matrix)
options(tibble.print_max = 150)

ifelse(dir.exists("/mnt/scc/project_workspace/"),
       cremwd <- "/mnt/scc/",
       cremwd <- "/restricted/projectnb/crem-bioinfo/")
calculations_path <- paste0(cremwd, "/project_workspace/18_09_25_kostas/calculations/analysis_B2.C11/")
data_path <- paste0(cremwd, "/project_workspace/18_09_25_kostas/calculations/cellranger_aggr/combined_data/outs/filtered_gene_bc_matrices_mex/hg19_tdtomato_genome_cellranger")
dir.create(calculations_path)

sparse.mat <-  Read10X(data_path)
summary(apply(sparse.mat, 2, sum)) # total counts per cell
summary(apply(sparse.mat, 1, sum)) # total counts per gene
sc <- CreateSeuratObject(raw.data = sparse.mat, min.cells = 0, min.genes = 200, normalization.method = "LogNormalize", scale.factor = 1e4, is.expr = 0, do.scale = F, do.center = F, save.raw = T, project = "sc")
sc

# add metadata
sc <- AddMetaData(sc, data.frame(row.names = sc@cell.names,
                                 orig.ident = factor(sub("^.*-", "", sc@cell.names), labels = c("B2", "C11"))),
                  "orig.ident")



# mito QC
mito.genes <- grep("^MT-", rownames(sc@data), value = T)
percent.mito <- Matrix::colSums(expm1(sc@data[mito.genes, ]))/Matrix::colSums(expm1(sc@data))
summary(percent.mito)
sc <- AddMetaData(sc, percent.mito, "percent.mito")
VlnPlot(sc, c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(paste0(calculations_path, "qc.prefilt.pdf"), width = 12, height = 4)
threshold <- 4*mad(percent.mito)
sc <- SubsetData(sc, subset.name = "percent.mito", accept.high = threshold)
VlnPlot(sc, c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(paste0(calculations_path, "qc.postfilt.",threshold,".pdf"), width = 12, height = 4)

# variable genes
pdf(file.path(calculations_path, "findVarGenes.pdf"), width = 8, height = 8)
sc <- FindVariableGenes(sc, 
                        mean.function = ExpMean,
                        dispersion.function = LogVMR,
                        x.low.cutoff = 0.2, x.high.cutoff = 4, 
                        y.cutoff = 1, do.contour = F, num.bin = 20)
dev.off()
length(sc@var.genes)
# scaling
sc <- ScaleData(sc, vars.to.regress = c("nUMI", "percent.mito"), scale.max = 40)
# dim. reduction
sc <- RunPCA(sc, 
             pc.genes = rownames(sc@scale.data),
             do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5)
# sc <- ProjectPCA(sc, do.center = T, replace.pc = T)
sc <- SetIdent(sc, ident.use = sc@meta.data$orig.ident)
PCAPlot(sc, 1, 2)
ggsave(paste0(calculations_path, "PCA.pdf"))


# clustering/modularity optimization
sc <- FindClusters(sc , reduction.type = "pca", resolution = seq(0.5, 2, 0.25), prune.SNN = 1/15, save.SNN = T, force.recalc = T)
colnames(sc@meta.data)
table(sc@meta.data$res.2)
table(sc@meta.data$res.1.75)
table(sc@meta.data$res.1.5)
table(sc@meta.data$res.1.25)
table(sc@meta.data$res.1)
table(sc@meta.data$res.0.75)
table(sc@meta.data$res.0.5)
# CrossTable(sc@meta.data$cell.type, sc@meta.data$res.0.75, format = "SPSS", prop.r=T, prop.c=T, prop.t=F, prop.chisq=F)
saveRDS(sc, paste0(calculations_path, "sc.Rds"))

# manifold learning
pdf(file.path(calculations_path, "elbow.pdf"), width = 8, height = 8)
PCElbowPlot(sc)
dev.off()
sc <- RunTSNE(sc, reduction.use = "pca",
              dims.use = 1:11, do.fast = T,
              add.iter = 0, perplexity = 30)
sc <- SetIdent(sc, ident.use = sc@meta.data$orig.ident)
TSNEPlot(sc, do.label = TRUE, pt.size = 0.5)
ggsave(paste0(calculations_path, "tSNE.pdf"))
plot.tsne.feat <- function(feat) {
  if(!(feat %in% rownames(sc@scale.data))) {return(cat(paste0(feat, " not in dataset\n")))}
  else {
    FeaturePlot(sc, cols.use = c("grey91", "red"), features.plot = feat, pt.size = 1.5)
    ggsave(paste0(calculations_path, "feat_tSNE_", feat, "_bigger1.5.lightgray.pdf"))
  }
}

interesting_genes <- c(
  "NKX2-1", "EGFP", "SFTA3", "SFTPB", "SOX2",
  "AFP", "ALB",
  "TFF1", # gut/gastric
  "CDX2", # intestinal
  "SCGB3A2", # type 2 cells
  "TOP2A", "AURKB", "MKI67", "BIRC5", # cell cycle
  "TP63", "KRT5",
  "SCGB1A1", 
  "SFTPC", "TFF2",
  "SOX9", 
  "HP",
  "FOXJ1",
  "BMP3", "FOXP2", "CPM", "FGF9",
  "NKD1",
  "GIF", "CLDN18", "MUC6",
  "TF", "APOA1", "CYP4F3", "CYP3A7",
  "PAX9", "KRT14",
  "SFTA2", "AQP4",
  "AXIN2",
  "TDTOMATO_EXTRA_IVS")
sapply(interesting_genes, plot.tsne.feat)
interesting_genes <- c(
  "SCGB3A2", "SCGB1A1", "NKX2-1", "SFTA3", "SOX9", "SOX2", "AXIN2", "TP63", "KRT5", "HP",
  "SFTPC", "SFTPB", "AFP", "ALB", "TFF1", "TFF2", "AURKB", "KIF11", "KIF14", "MKI67", "ETV4", "ETV5", "CYP4B1"
)
sapply(interesting_genes, plot.tsne.feat)

plot.tsne.clust <- function(clust) {
  sc <- SetAllIdent(sc, id = clust)
  TSNEPlot(sc, do.label = TRUE, pt.size = 1.5)
  ggsave(paste0(calculations_path, "clust_tSNE_", clust, ".pdf"))
}
plot.tsne.clust("orig.ident")
plot.tsne.clust("res.1.5")
plot.tsne.clust("res.1.25")
plot.tsne.clust("res.1")
plot.tsne.clust("res.0.75")
plot.tsne.clust("res.0.5")

saveRDS(sc, paste0(calculations_path, "sc.Rds"))


# cell cycle regression
cc.genes <- readLines(con = paste0(cremwd, "/reference_data/gene_sets/regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
sc <- CellCycleScoring(sc, s.genes = s.genes, g2m.genes = g2m.genes, 
                       set.ident = TRUE)
sc <- SetIdent(sc, ident.use = sc@meta.data$res.0.5)
RidgePlot(sc, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
          nCol = 2) # CC regression not needed
ggsave(paste0(calculations_path, "ridge.cell.cycle.pdf"))
# sc <- RunPCA(sc, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
# PCAPlot(sc)
# # correct differences among differentiating cells (G2M and S) but don't correct separation with non-differenciating cells (G1)
# sc@meta.data$CC.Difference <- sc@meta.data$S.Score - sc@meta.data$G2M.Score
# sc <- ScaleData(object = sc, vars.to.regress = "CC.Difference", display.progress = FALSE)
# # cell cycle effects strongly mitigated in PCA
# sc <- RunPCA(object = sc, pc.genes = rownames(sc@scale.data), genes.print = 10)
# # when running a PCA on cell cycle genes, actively proliferating cells
# # remain distinct from G1 cells however, within actively proliferating
# # cells, G2M and S phase cells group together
# sc <- RunPCA(object = sc, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
# PCAPlot(object = sc)
# PCElbowPlot(sc)
# sc <- RunTSNE(sc, reduction.use = "pca",
#                dims.use = 1:18, do.fast = T,
#                add.iter = 0, perplexity = 30)
# TSNEPlot(sc, do.label = TRUE, pt.size = 0.5)
# sc <- SetAllIdent(sc, id = "res.1.5")

fmwil <- function(sc, res) {
  file_out <- paste0(calculations_path, res, ".top30.clusterMarkers.wilcox.csv")
  sc <- SetAllIdent(sc, id = res)
  table(sc@ident)
  seurat.markers <- FindAllMarkers(sc, test.use = "wilcox", 
                                   only.pos = F,
                                   logfc.threshold = 0.25, do.print = T)
  seurat.markers.summary <- as.data.frame(seurat.markers %>% group_by(cluster) %>% top_n(30, abs(avg_logFC)))
  write.table(seurat.markers.summary, file = file_out, quote = F, row.names = F, col.names = T, sep = ",")
  return(list(all = seurat.markers, top = seurat.markers.summary))
}
wil.markers.0.5 <- fmwil(sc, "res.0.5")
wil.markers.0.5 <- fmwil(sc, "orig.ident")

plot.marker.heat <- function(clust) {
  file_out <- paste0(calculations_path, clust, ".top30.clusterMarkers.xlsx")
  sc <- SetAllIdent(sc, id = clust)
  seurat.markers <- FindAllMarkers(sc, only.pos = F, thresh.use = 0.2, do.print = T)
  sc.top10 <- seurat.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% top_n(10, avg_logFC)
  openxlsx::write.xlsx(seurat.markers, file_out, zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)
  tmp <- AverageExpression(object = sc, return.seurat = TRUE) # Average of all cells within cluster
  DoHeatmap(tmp, genes.use = c(interesting_genes, sc.top10$gene), col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#8B1A1A", group.label.rot = TRUE, group.cex = 0)
  ggsave(paste0(calculations_path, clust,  ".sc.avg.heatmap.pdf"))
  return(list(all = seurat.markers, top = sc.top10))
  
}
marker.orig.ident <- plot.marker.heat("orig.ident")
marker.res.1 <- plot.marker.heat("res.1")
# marker.res.0.75 <- plot.marker.heat("res.0.75")
# marker.res.0.5 <- plot.marker.heat("res.0.5")

sc <- SetAllIdent(sc, id = "orig.ident")
genes.use <- marker.orig.ident[[1]]  %>% group_by(cluster) %>% filter(p_val_adj < 0.05 & cluster == "B2" & abs(pct.1 - pct.2) > 0.1) %>% arrange(avg_logFC)
DoHeatmap(sc, genes.use = genes.use$gene, col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#8B1A1A", slim.col.label = T, group.label.rot = T, group.cex = 0)
ggsave(paste0(calculations_path, "orig.ident.sc.cells.heatmap.pdf"))


plot.marker.pair <- function(clust, ident.1, ident.2) {
  file_out <- paste0(calculations_path, clust, ".group_", ident.1, "vs", ident.2, "_pairwise.top50.clusterMarkers.xlsx")
  sc <- SetAllIdent(sc, id = clust)
  seurat.markers <- FindMarkers(sc, ident.1= ident.1, ident.2=ident.2, only.pos = F, thresh.use = 0.25, do.print = T)
  # seurat.markers$gene <- rownames(seurat.markers)
  # sc.top10 <- seurat.markers %>%  filter(p_val_adj < 0.05) %>% top_n(50, pct.1)
  openxlsx::write.xlsx(seurat.markers, file_out, zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)
  tmp <- AverageExpression(object = sc, return.seurat = TRUE) # Average of all cells within cluster
  DoHeatmap(tmp, genes.use = c(interesting_genes, sc.top10$gene), col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#8B1A1A", group.label.rot = TRUE, group.cex = 0)
  ggsave(paste0(calculations_path, clust, ".group_", ident.1, "vs", ident.2, "_pairwise.sc.avg.heatmap.pdf"))
}
plot.marker.pair("orig.ident", "B2", "C11")
# 0 is lung AT2 cells
# 3 might be myeloid progenitors
plot.marker.pair("res.1", 5, NULL)
# 5 is fetal liver

sc.markers.wilcox.cell.type <- fmwil(sc, "cell.type")
sc.markers.wilcox.res.0.5 <- fmwil(sc, "res.0.5")
sc.markers.wilcox.res.0.75 <- fmwil(sc, "res.0.75")

openxlsx::write.xlsx(
  seurat.markers,
  file.path(calculations_path, "seurat.markers.wilcox.xlsx"), 
  zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T
)


dim(sc@raw.data)
dim(sc@data)
dim(sc@scale.data)
export.data <- sc@raw.data[-grep("^MT-", rownames(sc@raw.data), value = F), colnames(sc@data)]
export.norm.data <- sc@data[-grep("^MT-", rownames(sc@data), value = F), colnames(sc@data)]
export.genes <- rownames(sc@raw.data)[-grep("^MT-", rownames(sc@raw.data), value = F)]
export.norm.genes <- rownames(sc@raw.data)[-grep("^MT-", rownames(sc@raw.data), value = F)]
export.barcodes <- colnames(sc@data)
all.genes <- read.table(paste0(data_path, "/genes.tsv"))
ensembl.genes <- all.genes[match(export.genes, all.genes$V2), 1]
write.table(as.matrix(export.data), paste0(calculations_path, "spring.data.tsv"), sep = "\t", quote = F, col.names = T, row.names = T )
write.table(as.matrix(export.norm.data), paste0(calculations_path, "spring.norm.data.tsv"), sep = "\t", quote = F, col.names = T, row.names = T )
write.table(export.genes, paste0(calculations_path, "spring.genes.tsv"), sep = "\t", quote = F, col.names = F, row.names = F )
write.table(ensembl.genes, paste0(calculations_path, "loupe_genes.csv"), sep = "\t", quote = F, col.names = "Gene", row.names = F )
write.table(export.barcodes, paste0(calculations_path, "barcodes.csv"), sep = "\t", quote = F, col.names = "Barcode", row.names = F )


groupings <- t(sc@meta.data)[c(3,5:11,14),] # select metadata
write.table(groupings, paste0(calculations_path, "spring.groupings.csv"), sep = ",", quote = F, col.names = F, row.names = T )
save.image(paste0(calculations_path, ".RData"))

# for Heatmap of GSEA core set genes:
# grep -w Yes HALL*xls > core.set.ordered.txt
# sed "s/HALLMARK_//" core.set.ordered.txt > test
# sed "s/.xls:row_/\t/" test > core.set.ordered.txt

gsea.core.set <- scan(paste0(calculations_path, "core.set.uniq.txt"), what = "character") # too many ~ 1000
gsea.core.set <- scan(paste0(calculations_path, "core.set.top10.uniq.txt"), what = "character") # still too many and most are not expressed in most cells
gsea.core.set <- read.table(paste0(calculations_path, "core.set.ordered.txt"), sep = "\t", stringsAsFactors = F)
gsea.ord <- read.table(paste0(calculations_path, "gsea_report_for_na_pos_1537282701433.xls"), header = T, sep = "\t", stringsAsFactors = F)  # get set NES
gsea.ord$NAME <- sub("HALLMARK_", "", gsea.ord$NAME)
gsea.ord <- gsea.ord[gsea.ord$FDR.q.val < 0.05, ] # only those significant
gsea.core.set <- gsea.core.set[which(gsea.core.set$V1 %in% gsea.ord$NAME),]   # subset core set genes only from gene sets that were significant

gsea.core.set <- gsea.core.set %>%
  mutate(set =  factor(V1, levels = gsea.ord$NAME)) %>%
  arrange(set)

seurat.markers.pos <- FindMarkers(sc, ident.1= "E", ident.2="F", only.pos = T, min.pct=0.05, do.print = T) # subset those that are upreg and detected in a min.pct 
gsea.core.subset <- gsea.core.set[gsea.core.set$V3 %in% rownames(seurat.markers.pos), ]  # subset the initial GSEA core set
write.table(gsea.core.subset, paste0(calculations_path, "gsea.core.subset.txt"), col.names=F, row.names= F, sep= "\t", quote = F)
sc <- SetAllIdent(sc, id = 'orig.ident')
DoHeatmap(sc, genes.use = gsea.core.subset$V3, use.scaled = F, slim.col.label = TRUE, remove.key = TRUE, col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#CD3333", cex.row=1)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.norm_notScaled.pdf"), width = 10, height = 10)

sc.average <- AverageExpression(sc, return.seurat = TRUE) # Average of all cells within cluster
DoHeatmap(sc.average, genes.use = gsea.core.subset$V3, use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#CD3333", cex.row=3)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.average.pdf"), width = 1, height = 10)

DoHeatmap(sc, genes.use = gsea.core.subset$V3, cells.use = sc@cell.names[which(sc@cell.names %in%  sample(sc@cell.names, 1000, replace=F))], use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, cex.row=1)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.col_orig.subsample.pdf"), width = 10, height = 10)


library(pheatmap)
data.use <- sc@scale.data
cells.use = sc@cell.names[which(sc@cell.names %in%  sample(sc@cell.names, 1000, replace=F))]
genes.use <- gsea.core.subset$V3
cells.ident <- sc@ident[cells.use]
data.use <- data.use[genes.use, cells.use, drop = FALSE]
data.use <- MinMax(data = data.use, min = -2.5, max = 2.5)
row.names(data.use) <- make.unique(rownames(data.use))
pdf(file.path(calculations_path, "heat.p.3.pdf"), width=10, height=14)
pheatmap(data.use,
         show_colnames = F,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(50),
         cluster_cols = F,
         cluster_rows = F,
         fontsize_row = 8,
         annotation_col = data.frame(row.names = cells.use, sample = cells.ident),
         annotation_row = data.frame(row.names = make.unique(genes.use), hallmark = factor(gsea.core.subset$V1, levels = unique(gsea.core.subset$V1))),
         labels_row = genes.use
)
dev.off()

library(gplots)
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(50) # divergent red blue
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # gradient

heatmap.2(mat, trace="none", dendrogram = "none", col = rev(hmcol), margin=c(5, 5), main="Title")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), main="Title")

distsRL <- dist(t(rlogMat))
mat <- as.matrix(distsRL)


sc <- SetAllIdent(sc, id = "orig.ident")
seurat.markers <- FindMarkers(sc, ident.1 = "E", logfc.threshold = 0, only.pos = F)
seurat.markers <- seurat.markers[-grep("^MT-", rownames(seurat.markers), value = F), ]
openxlsx::write.xlsx(seurat.markers, paste0(calculations_path, "allgenes.EvsF_pairwise.xlsx"), zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)

ggplot(data = data.use, mapping = aes(x = cell, 
                                      y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low, 
                                                                                                         mid = col.mid, high = col.high, name = "Expression", 
                                                                                                         guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
  scale_y_discrete(position = "right", labels = rev(genes.use)) + 
  theme(axis.line = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
        axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
        axis.title.x = element_blank()) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank())




############################################
############################################
############################################
############################################
sc <- readRDS(paste0(calculations_path, "sc.Rds"))
sc <- UpdateSeuratObject(sc)
table(sc$orig.ident)
table(sc$annot)
msig <- scan(paste0(cremwd, "/reference_data/gene_sets/lung.matur.differ.MSig.txt"), what="", sep="\n") # Read in the data
msig <- sapply(msig, function(k) gsub(",$", "", k))
msig <- strsplit(msig, ",[[:space:]]+")  # Separate elements by a comma and one or more whitepace
names(msig) <- sapply(msig, `[[`, 1)  # subset the first vector element and set it as the list element name
msig <- lapply(msig, `[`, -1) # Remove the first vector element from each list element
msig <- msig[1:2]
msig <- unlist(msig)
genes.plot <- unique(c( "SFTPC", msig, "SFTPB", "SFTPD", "SFTPA1", "SFTPA2", "NKX2-1", "MKI67", "TM4SF1" ,"SQSTM1", "KRT5", "SCGB1A1", "FOXJ1","PTPRC", "CDH5", "EPCAM", "COL1A1","ATP6V0B"  ))
DefaultAssay(sc) <- "RNA"
Idents(sc) <- "orig.ident"
project <- sub(".*workspace/(.*tas).*", "\\1", calculations_path)
p <- DotPlot(sc, features = rev(genes.plot),  assay = "RNA") 
p <- p & ggmin::theme_powerpoint()  & RotatedAxis() & theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(paste0(calculations_path, "Dot.JCL.", project, ".pdf"), p, width = 10, height = 4.5)
p <- VlnPlot(sc, features = genes.plot, group.by = "orig.ident", pt.size = 0.1, assay = "RNA")
p <- p & ggmin::theme_powerpoint()  & RotatedAxis() & theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(paste0(calculations_path, "Vln.JCL.", project, ".pdf"), p, width = 15, height = 15)

library(future)
plan("multiprocess", workers = 16)
options(future.globals.maxSize = 4000 * 1024^2)
Idents(sc) <- "orig.ident"
fmas.pair <- function(sc, res, id1, id2, subsets, suffix) {
    file_out <- paste0(calculations_path, "DEG.pairwise.", suffix, ".", res, ".xlsx")
    seurat.markers <- FindMarkers(sc, test.use = "MAST", only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0, verbose = F, group.by = res,  ident.1 = id1, ident.2 = id2, subset.ident = subsets)
    seurat.markers.summary <- as.data.frame(seurat.markers) %>%  rownames_to_column('gene') %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_logFC)) %>% dplyr::filter(row_number() %in% 1:20 | row_number() %in% (n()-19):n())
    openxlsx::write.xlsx(seurat.markers, file_out, zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)
    return(list(all = seurat.markers, top = seurat.markers.summary))
}
mas.pair1 <- fmas.pair(sc, "orig.ident", id1 = "C11", id2 = "B2", subsets = NULL, suffix = "C11_vs_B2")
mas.pair1[[1]] %>% rownames_to_column("row.names") %>% filter(p_val_adj < 0.05) %>% group_by(avg_logFC > 0) %>% count(`avg_logFC > 0`)

saveRDS(mas.pair1 , paste0(calculations_path, "sc.markers.mas.", resDEG, ".Rds"))

mas.pair1 <- readxl::read_excel(paste0(calculations_path, "DEG.pairwise.C11_vs_B2.orig.ident.xlsx")) 
mas.pair1 %>%  filter(avg_logFC < 0.15) %>%  filter(p_val_adj < 0.15) %>%
  summarise(upDEG = sum(avg_logFC>0),
         dnDEG = sum(avg_logFC<0))

fmas.pair.final <- function(sc, res, id1, id2, subsets, suffix) {
    file_out <- paste0(calculations_path, "DEG.pairwise.", suffix, ".", res, ".pct.0.10.logfc.0.15.FDR.0.05.xlsx")
    seurat.markers <- FindMarkers(sc, test.use = "MAST", only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.15, verbose = F, group.by = res,  ident.1 = id1, ident.2 = id2, subset.ident = subsets)
    seurat.markers.summary <- as.data.frame(seurat.markers) %>%  rownames_to_column('gene') %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(desc(avg_logFC)) 
    openxlsx::write.xlsx(seurat.markers.summary, file_out, zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)
    return(list(all = seurat.markers, top = seurat.markers.summary))
}


mas.pair1.final <- fmas.pair.final(sc, "orig.ident", id1 = "C11", id2 = "B2", subsets = NULL, suffix = "C11_vs_B2")
mas.pair1.final[[2]] %>%   
  summarise(upDEG = sum(avg_logFC>0),
         dnDEG = sum(avg_logFC<0))
# upDEG dnDEG
# 1    73   151