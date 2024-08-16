#Rustenhoven meninges--------------------------------------------------------------------------------------------------
#SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate seurat4

#Seurat
R
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)

library(purrr)
library(tibble)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven")

MO1_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/MO1")
MO1 = CreateSeuratObject(MO1_10x, project="MO1")

MO2_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/MO2")
MO2 = CreateSeuratObject(MO2_10x, project="MO2")

MY1_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/MY1")
MY1 = CreateSeuratObject(MY1_10x, project="MY1")

MY2_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/MY2")
MY2 = CreateSeuratObject(MY2_10x, project="MY2")

MO_merged <- merge(MO1, y = MO2, add.cell.ids = c("MO1", "MO2"), project = "MO")

MY_merged <- merge(MY1, y = MY2, add.cell.ids = c("MY1", "MY2"), project = "MY")

mem_merged <- merge(MO_merged, y = MY_merged)

mem_merged[["percent.mt"]] <- PercentageFeatureSet(mem_merged, pattern = "^mt-")

mean(mem_merged$nFeature_RNA) + 3 * sd(mem_merged$nFeature_RNA)
mean(mem_merged$nCount_RNA) + 3 * sd(mem_merged$nCount_RNA)
mean(mem_merged$percent.mt) + 3 * sd(mem_merged$percent.mt)

table(mem_merged$orig.ident)
mem_merged <- subset(mem_merged, subset = nFeature_RNA > 200 & nCount_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 35000 & percent.mt < 25)
table(mem_merged$orig.ident)

png(file = paste0("VlnPlot_QC_subset_mem_merged_data.png"), width = 15, height = 17, units = "cm", res = 500)
p <- VlnPlot(mem_merged, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()

mem_merged[["RNA"]] <- split(mem_merged[["RNA"]], f = mem_merged$orig.ident)
mem_merged <- NormalizeData(mem_merged)
mem_merged <- FindVariableFeatures(mem_merged)
mem_merged <- ScaleData(mem_merged, verbose = FALSE)
mem_merged <- RunPCA(mem_merged, npcs = 30, verbose = FALSE)
#New seurat5 integration method
options(future.globals.maxSize = 8000 * 1024^2)
mem_merged <- IntegrateLayers(object = mem_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
mem_merged[["RNA"]] <- JoinLayers(mem_merged[["RNA"]])
#Run UMAP with desired dimension
mem_merged <- FindNeighbors(mem_merged, reduction = "integrated.cca", dims = 1:30)
mem_merged <- FindClusters(mem_merged, resolution = c(0.5))
mem_merged <- RunUMAP(mem_merged, reduction = "integrated.cca", dims = 1:30)

save(mem_merged, file = "mem_merged.Robj")
load(file = "mem_merged.Robj", verbose = TRUE)


annotations <- read.csv("annotations_Mus_EnsDb.csv")
#FindAllMarkers for single resolution
Idents(object = mem_merged) <- "seurat_clusters"	#Change to desired resolution
mem_merged.markers <- FindAllMarkers(mem_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.csv(mem_merged.markers,'mem_merged_merged_findallmarkers_0.5_30pc.csv')


#Annotation
Idents(object = mem_merged) <- mem_merged$seurat_clusters
mem_merged <- RenameIdents(mem_merged, 
	`0` = "Monocyte/Macrophages", 
	`1` = "B cells", 
	`2` = "T cells",
    `3` = "Monocyte/Macrophages", 
	`4` = "Monocyte/Macrophages", 
	`5` = "Neutrophils", 
	`6` = "Fibroblasts", 
	`7` = "Neutrophils", 
	`8` = "T cells", 
	`9` = "Erythroid cells", 
	`10` = "Monocyte/Macrophages", 
	`11` = "ILC2s", 
	`12` = "Dendritic cells", 
	`13` = "B cells", 
	`14` = "Proliferating cells", 
	`15` = "T cells", 
	`16` = "B cells", 
	`17` = "NK cells", 
	`18` = "Monocyte/Macrophages", 
	`19` = "Endothelial cells", 
	`20` = "Erythroid cells", 
	`21` = "Dendritic cells", 
	`22` = "Monocyte/Macrophages",
	`23` = "Mast cells",
	`24` = "B cells",
	`25` = "Monocyte/Macrophages",
	`26` = "Neutrophils",
	`27` = "Dendritic cells",
	`28` = "Unknown",
	`29` = "Unknown-2",
	`30` = "Schwann cells",
	`31` = "Mural cells")
mem_merged$cell_ident <- Idents(object = mem_merged)

total_order <- c("B cells","T cells","NK cells","ILC2s", "Dendritic cells","Monocyte/Macrophages","Neutrophils","Mast cells","Erythroid cells", "Fibroblasts","Endothelial cells","Mural cells","Proliferating cells","Schwann cells","Unknown","Unknown-2")

#Reorganize levels 
Idents(object = mem_merged) <- mem_merged$cell_ident
mem_merged$cell_ident <- factor(mem_merged$cell_ident, levels=total_order)


#Annotation markers---------------------------------------------------------
mono <- c("Mrc1", "C1qc", "Aif1", "Csf1r")
b <- c("Cd79a", "Cd19", "Ly6d", "Ighm")
t <- c("Cd3e", "Cd8a", "Ccr5")
ilc <- c("Gata3", "Il1rl1", "Il17rb")
neutro <- c("Csf3r", "S100a8", "Ly6g")
fibro <- c("Col1a1", "Dcn")
eryth <- c("Alas2", "Car2", "Hba-a1")
dendrit <- c("Cd209a", "H2-Aa", "Itgax", "Xcr1")
prolif <- c("Mki67", "Top2a", "Birc5", "Dbi")
nk <- c("Klrb1b", "Ncr1")
endo <- c("Pecam1", "Tie1", "Podxl")
mast <- c("Cpa3", "Mcpt4", "Fcer1a")
unknown2 <- c("Neurod1", "Neurod4", "Slc6a11", "Nrcam", "Nrxn3", "Clu")
schwann <- c("Plp1", "Sox10", "Atp1a2")
mural <- c("Rgs5", "Rgs4", "Pdgfrb")

total_markers <- c(b, t, nk, ilc, dendrit, mono, neutro, mast, eryth, fibro, endo, mural, prolif, schwann, unknown2)

png(file = "Annotations_dot_total.png", width = 37, height = 17, units = "cm", res = 900)
DotPlot(mem_merged, features = total_markers, group.by = "cell_ident", dot.scale = 8, dot.min = 0.01, scale=TRUE) + RotatedAxis()
dev.off()

#UMAP---------------------------------------------------------
Idents(object = mem_merged) <- "seurat_clusters"
png(file = "Cell_ident_mem_0.5_30_small.png", width = 20, height = 17, units = "cm", res = 700)
DimPlot(mem_merged, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.8)
dev.off()
#UMAP Final---------
#Solo cell ident
Idents(object = mem_merged) <- "cell_ident"
png(file = "Cell_ident_mem_0.5_small_legend_cellident.png", width = 30, height = 21, units = "cm", res = 700)
DimPlot(mem_merged, reduction = "umap", label = FALSE, repel = TRUE, label.size = 4, pt.size = 0.05)
dev.off()
png(file = "Cell_ident_mem_0.5_small_label_cellident.png", width = 30, height = 21, units = "cm", res = 700)
DimPlot(mem_merged, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.05)
dev.off()
png(file = "Cell_ident_mem_0.5_small_nolegend_cellident.png", width = 25, height = 25, units = "cm", res = 700)
DimPlot(mem_merged, reduction = "umap", label = FALSE, repel = TRUE, label.size = 4, pt.size = 0.05) + theme(legend.position = "none")
dev.off()



#Markers of interest---------------------------------------------------------
png(file = "Men_glut_asma_genes.png", width = 20, height = 10, units = "cm", res = 500)
VlnPlot(mem_merged, features = c("Slc2a1", "Acta2"))
dev.off()
png(file = "Men_glut_gene.png", width = 12, height = 10, units = "cm", res = 500)
VlnPlot(mem_merged, features = c("Slc2a1")) + theme(legend.position = 'none')
dev.off()
png(file = "Men_asma_gene.png", width = 12, height = 10, units = "cm", res = 500)
VlnPlot(mem_merged, features = c("Acta2")) + theme(legend.position = 'none')
dev.off()

png(file = "Mem_merged_glut_all_genes.png", width = 22, height = 30, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Acta2", "Tagln", "Rgs5", "Pdgfrb"))
dev.off()

png(file = "Mem_merged_clus_5_7_genes.png", width = 20, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Ly6g", "Csf3r"))
dev.off()

#Blend
png(file = "Mem_merged_glut_asma_genes_b0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Acta2"), blend = TRUE, blend.threshold = 0)
dev.off()
png(file = "Mem_merged_glut_tagln_genes_b0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Tagln"), blend = TRUE, blend.threshold = 0)
dev.off()
png(file = "Mem_merged_glut_rgs5_genes_b0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Rgs5"), blend = TRUE, blend.threshold = 0)
dev.off()
png(file = "Mem_merged_glut_pdgfrb_genes_b0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1","Pdgfrb"), blend = TRUE, blend.threshold = 0)
dev.off()


#Check with no QC subsetting----------------------------------------------------------------------
load(file = "mem_merged_noQC.Robj")
Idents(object = mem_merged) <- "seurat_clusters"

png(file = "Cell_ident_mem_0.5_30_small_noQC.png", width = 20, height = 17, units = "cm", res = 700)
DimPlot(mem_merged, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.8)
dev.off()

png(file = "Men_merged_glut_asma_genes_noQC_blend0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Acta2"), blend = TRUE, blend.threshold = 0)
dev.off()

png(file = "Men_merged_glut_tagln_genes_noQC_blend0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Tagln"), blend = TRUE, blend.threshold = 0)
dev.off()

png(file = "Men_merged_glut_rgs5_genes_noQC_blend0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Rgs5"), blend = TRUE, blend.threshold = 0)
dev.off()

png(file = "Men_merged_glut_pdgfrb_genes_noQC_blend0.png", width = 30, height = 10, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Pdgfrb"), blend = TRUE, blend.threshold = 0) 
dev.off()

png(file = "Men_merged_glut_all_genes.png", width = 22, height = 30, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Acta2", "Tagln", "Rgs5", "Pdgfrb"))
dev.off()

#Check other MY MO object (endothelial & mural cells)-------------------------------------------------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/MY_MO")

MO_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/data/MO")
MO = CreateSeuratObject(MO_10x, project="MO")

MY_10x <- Read10X("/work/biernaskie_lab/apun/nilesh_glut1/meninges/rustenhoven/data/MY")
MY = CreateSeuratObject(MY_10x, project="MY")

mem_merged <- merge(MO, y = MY, add.cell.ids = c("MO", "MY"), project = "men")

mem_merged[["percent.mt"]] <- PercentageFeatureSet(mem_merged, pattern = "^mt-")

table(mem_merged$orig.ident)

#mem_merged[["RNA"]] <- split(mem_merged[["RNA"]], f = mem_merged$orig.ident)
mem_merged <- NormalizeData(mem_merged)
mem_merged <- FindVariableFeatures(mem_merged)
mem_merged <- ScaleData(mem_merged, verbose = FALSE)
mem_merged <- RunPCA(mem_merged, npcs = 30, verbose = FALSE)
#New seurat5 integration method
options(future.globals.maxSize = 8000 * 1024^2)
mem_merged <- IntegrateLayers(object = mem_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
mem_merged[["RNA"]] <- JoinLayers(mem_merged[["RNA"]])
#Run UMAP with desired dimension
mem_merged <- FindNeighbors(mem_merged, reduction = "integrated.cca", dims = 1:30)
mem_merged <- FindClusters(mem_merged, resolution = c(0.5))
mem_merged <- RunUMAP(mem_merged, reduction = "integrated.cca", dims = 1:30)

png(file = "Men_merged_glut_all_genes.png", width = 22, height = 30, units = "cm", res = 500)
FeaturePlot(mem_merged, features = c("Slc2a1", "Acta2", "Tagln", "Rgs5", "Pdgfrb"))
dev.off()

save(mem_merged, file = "mem_merged_mymo.Robj")



