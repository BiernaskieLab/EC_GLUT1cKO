#Habib 5xFAD-----------------------------------------------------------------------------------------------------------------------------------

#SLURM--------
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate seurat4

#Seurat----------------------
R
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)

library(purrr)
library(tibble)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad")

counts <- data.table::fread("GSE143758_Admouse_Hippocampus_7m_AllNuclei_UMIcounts.txt.gz", data.table = FALSE)

#Check and fix row names
head(counts[, 1:5], n = 10)
row.names(counts) <- counts[, 1]
counts <- counts[, -1]

fad = CreateSeuratObject(counts, project="5xFAD")

#Fix metadata from row names-------------
extract_text <- function(row_name) {
  last_underscore_index <- max(gregexpr("_", row_name)[[1]])
  extracted_text <- substr(row_name, 1, last_underscore_index)
  return(extracted_text)
}

extracted_text <- character(length = nrow(fad@meta.data))
# Apply the function to each row name
for (i in seq_len(nrow(fad@meta.data))) {
  extracted_text[i] <- extract_text(rownames(fad@meta.data)[i])
}

split_text <- strsplit(extracted_text, "[._]")
split_df <- do.call(rbind, split_text)
colnames(split_df) <- c("lysis", "batch", "sample")

fad@meta.data <- cbind(fad@meta.data, split_df)


#Remove bad orig.ident column
fad@meta.data <- subset(fad@meta.data, select = -orig.ident)

#Find genotype--------------------------------
fad@meta.data$genotype <- NA


wt_indices <- grep("Wt|WT", fad@meta.data$sample)
fad@meta.data$genotype[wt_indices] <- "WT"

ad_indices <- grep("AD|Untreated", fad@meta.data$sample)
fad@meta.data$genotype[ad_indices] <- "AD"

#Add age sex----------------------------

fad@meta.data$age <- "7 mo"
fad@meta.data$sex <- "male"


#Processing-------------------------------------------------------------------------------

fad[["percent.mt"]] <- PercentageFeatureSet(fad, pattern = "^mt-")

mean(fad$nFeature_RNA) + 4 * sd(fad$nFeature_RNA)
mean(fad$nCount_RNA) + 4 * sd(fad$nCount_RNA)
mean(fad$percent.mt) + 4 * sd(fad$percent.mt)

table(fad$genotype)

png(file = paste0("VlnPlot_QC_fad_data.png"), width = 15, height = 17, units = "cm", res = 500)
p <- VlnPlot(fad, group.by = "genotype", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()

fad <- subset(fad, subset = nFeature_RNA > 200 & nCount_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 4500 & percent.mt < 4)
table(fad$genotype)

png(file = paste0("VlnPlot_QC_subset_fad_data.png"), width = 15, height = 17, units = "cm", res = 500)
p <- VlnPlot(fad, group.by = "genotype", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()


fad[["RNA"]] <- split(fad[["RNA"]], f = fad$genotype)
fad <- NormalizeData(fad)
fad <- FindVariableFeatures(fad)
fad <- ScaleData(fad, verbose = FALSE)
fad <- RunPCA(fad, npcs = 30, verbose = FALSE)

#New seurat5 integration method
options(future.globals.maxSize = 8000 * 1024^2)
fad <- IntegrateLayers(object = fad, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
fad[["RNA"]] <- JoinLayers(fad[["RNA"]])
#Run UMAP with desired dimension
fad <- FindNeighbors(fad, reduction = "integrated.cca", dims = 1:30)
fad <- FindClusters(fad, resolution = c(0.5))
fad <- RunUMAP(fad, reduction = "integrated.cca", dims = 1:30)


png(file = "5xfad_cluster_0.5.png", width = 17, height = 15, units = "cm", res = 700)
DimPlot(fad, reduction = "umap")
dev.off()


save(fad, file = "5xfad.Robj")
load(file = "5xfad.Robj", verbose = TRUE)

#Annotation--------------------------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad")

#FindAllMarkers for single resolution (best for single conditions and small clusters with insufficient cells)
Idents(object = fad) <- "seurat_clusters"	#Change to desired resolution
fad.markers <- FindAllMarkers(fad, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.csv(fad.markers,'5xfad_merged_findallmarkers_0.5_30pc.csv')


Idents(object = fad) <- fad$seurat_clusters
fad <- RenameIdents(fad, `1` = "Astrocyte", `15` = "Oligodendrocyte", `17` = "Ependymal", 
`19` =  "Microglia", `20` = "Astrocyte", `22` = "Astrocyte")
#Need this line too for renaming to actually rename the old column
fad$cell_ident <- Idents(object = fad)

total_markers <- c("Ccdc153", "Rarres2", "Foxj1", "Dnah12", "Odf3b", "Slc1a2", "Clu", "Ndrg2", "Gja1", "Gfap", "Mag", "Mog", "Mbp", "Sox10", "Opalin", "Hexb", "Cx3cr1", "Csf1r")
total_order <- c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")
png(file = "Annotations_dot_total_5xFAD.png", width = 20, height = 15, units = "cm", res = 900)
DotPlot(fad, features = total_markers, group.by = "cell_ident", idents = total_order, dot.scale = 8, dot.min = 0.01, scale=TRUE) + RotatedAxis()
dev.off()


#GSEA DEG prep 
for(cell_type in unique(fad@meta.data$cell_ident)){ 
	if (cell_type %in% c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")) {
		Idents(fad) <- "cell_ident"
		ofinterest_cell = subset(fad, idents = cell_type)
		Idents(ofinterest_cell) <- "genotype"
		markers <- FindMarkers(ofinterest_cell, ident.1 = "AD", ident.2 = "WT", min.pct = 0, logfc.threshold = 0) %>%
			rownames_to_column(var = "gene") %>% 
			left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
			dplyr::arrange(p_val_adj)
		write.csv(markers,file = paste("5xFAD_", cell_type,"_markers_fem.gsea.csv", sep = ""))
	}
}