
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
calculations_path <- paste0(cremwd, "/project_workspace/18_09_25_kostas/calculations/r_analysis.C11/")
data_path <- paste0(cremwd, "/project_workspace/18_09_25_kostas/calculations/cellranger_count_with_tdTomato/C11/outs/filtered_gene_bc_matrices/hg19_tdtomato_genome_cellranger")
dir.create(calculations_path)

sparse.mat <-  Read10X(data_path)
summary(apply(sparse.mat, 2, sum)) # total counts per cell
summary(apply(sparse.mat, 1, sum)) # total counts per gene
sc <- CreateSeuratObject(raw.data = sparse.mat, min.cells = 0, min.genes = 200, normalization.method = "LogNormalize", scale.factor = 1e4, is.expr = 0, do.scale = F, do.center = F, save.raw = T, project = "sc")
sc

# add metadata
# sc <- AddMetaData(sc, data.frame(row.names = sc@cell.names,
#                                  orig.ident = factor(sub("^.*-", "", sc@cell.names), labels = c("E", "F"))),
#                   "orig.ident")



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

sc <- FindVariableGenes(sc, 
                        mean.function = ExpMean,
                        dispersion.function = LogVMR,
                        x.low.cutoff = 0.2, x.high.cutoff = 4, 
                        y.cutoff = 2, do.contour = F, num.bin = 30)
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
    FeaturePlot(sc, features.plot = feat, pt.size = 0.5)
    ggsave(paste0(calculations_path, "feat_tSNE_", feat, ".pdf"))
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

plot.tsne.clust <- function(clust) {
  sc <- SetAllIdent(sc, id = clust)
  TSNEPlot(sc, do.label = TRUE, pt.size = 0.5)
  ggsave(paste0(calculations_path, "clust_tSNE_", clust, ".pdf"))
}
plot.tsne.clust("orig.ident")
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
  sc.top30 <- seurat.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% top_n(30, avg_logFC)
  openxlsx::write.xlsx(seurat.markers, file_out, zoom = 100, asTable = TRUE, tableStyle = "none", rowNames = T)
  tmp <- AverageExpression(object = sc, return.seurat = TRUE) # Average of all cells within cluster
  DoHeatmap(tmp, genes.use = c(interesting_genes, sc.top10$gene), col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#8B1A1A", group.label.rot = TRUE, group.cex = 0)
  ggsave(paste0(calculations_path, clust,  ".sc.avg.heatmap.pdf"))
  return(list(all = seurat.markers, top = sc.top30))
  
}
marker.orig.ident <- plot.marker.heat("orig.ident")
marker.res.1 <- plot.marker.heat("res.1")
marker.res.0.75 <- plot.marker.heat("res.0.75")
marker.res.0.5 <- plot.marker.heat("res.0.5")

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
plot.marker.pair("res.1", 0, 3)
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


groupings <- t(sc@meta.data)[c(3,5:8,14),] # select metadata
write.table(groupings, paste0(calculations_path, "spring.groupings.csv"), sep = ",", quote = F, col.names = F, row.names = T )
save.image(paste0(calculations_path, ".RData"))

# spring server: 
gsea.core.set <- scan(paste0(calculations_path, "core.set.uniq.txt"), what = "character") # too many ~ 1000
gsea.core.set <- scan(paste0(calculations_path, "core.set.top10.uniq.txt"), what = "character") # still too many and most are not expressed in most cells
gsea.core.set <- read.table(paste0(calculations_path, "core.set.ordered.txt"), sep = "\t", stringsAsFactors = F)

# grep -w Yes HALL*xls > core.set.ordered.txt
# sed "s/HALLMARK_//" core.set.ordered.txt > test
# sed "s/.xls:row_/\t/" test > core.set.ordered.txt

seurat.markers.pos <- FindMarkers(sc, ident.1= "E", ident.2="F", only.pos = T, min.pct=0.05, do.print = T) # subset those that are upreg and detected in a min.pct 
gsea.core.subset <- gsea.core.set[gsea.core.set$V3 %in% rownames(seurat.markers.pos), ]  # subset the initial GSEA core set
write.table(gsea.core.subset, paste0(calculations_path, "gsea.core.subset.txt"), col.names=F, row.names= F, sep= "\t", quotes = F)
sc <- SetAllIdent(sc, id = 'orig.ident')
DoHeatmap(sc, genes.use = gsea.core.subset$V3, use.scaled = F, slim.col.label = TRUE, remove.key = TRUE, col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#CD3333", cex.row=1)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.norm_notScaled.pdf"), width = 10, height = 10)

sc.average <- AverageExpression(sc, return.seurat = TRUE) # Average of all cells within cluster
DoHeatmap(sc.average, genes.use = gsea.core.subset$V3, use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, col.low = "#0000EE", col.mid = "#FFFFFF", col.high = "#CD3333", cex.row=3)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.average.pdf"), width = 1, height = 10)

DoHeatmap(sc, genes.use = gsea.core.subset$V3, cells.use = sc@cell.names[which(sc@cell.names %in%  sample(sc@cell.names, 1000, replace=F))], use.scaled = T, slim.col.label = TRUE, remove.key = TRUE, cex.row=1)
ggsave(file.path(calculations_path, "gsea.coreset.heat.subset.col_orig.subsample.pdf"), width = 10, height = 10)

library(pheatmap)