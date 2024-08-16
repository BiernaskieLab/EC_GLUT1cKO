#Lee Apoe----------------------------------------------------------------

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

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/")
apoe <- readRDS("APOE_x.rds")
DefaultAssay(apoe) <- "RNA"

#Annotation--------------------------------------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/")

Idents(object = apoe) <- "seurat_clusters"
apoe.markers <- FindAllMarkers(apoe, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.csv(apoe.markers,'apoe_findallmarkers_RNA.csv')

Idents(object = apoe) <- "seurat_clusters"
png(file = "Apoe_cluster_0.5.png", width = 17, height = 15, units = "cm", res = 700)
DimPlot(apoe, reduction = "umap", label = TRUE)
dev.off()


#Annotation --------------------------------------------
DefaultAssay(apoe) <- "RNA"
Idents(object = apoe) <- apoe$seurat_clusters
apoe <- RenameIdents(apoe, `0` = "Oligodendrocyte", `1` = "Microglia", `2` = "Astrocyte", `3` = "Astrocyte", `5` = "Astrocyte", 
`6` = "Microglia", `9` = "Microglia", `10` = "Oligodendrocyte", `12` = "Ependymal", `14` = "Microglia",
`20` = "Oligodendrocyte", `24` = "Oligodendrocyte")
apoe$cell_ident <- Idents(object = apoe)

total_markers <- c("Ccdc153", "Rarres2", "Foxj1", "Dnah12", "Odf3b", "Slc1a2", "Clu", "Ndrg2", "Gja1", "Gfap", "Mag", "Mog", "Mbp", "Sox10", "Opalin", "Hexb", "Cx3cr1", "Csf1r")
total_order <- c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")
png(file = "Annotations_dot_total_apoe.png", width = 20, height = 15, units = "cm", res = 900)
DotPlot(apoe, features = total_markers, group.by = "cell_ident", idents = total_order, dot.scale = 8, dot.min = 0.01, scale=TRUE) + RotatedAxis()
dev.off()

Idents(object = apoe) <- "age"
apoe_old <- subset(apoe, idents = c("12mo", "24mo"))

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/deg")

#GSEA DEG prep --------------------------------------------
for(cell_type in unique(apoe_old@meta.data$cell_ident)){ 
	if (cell_type %in% c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")) {
	
		ofinterest = apoe_old
	
		Idents(ofinterest) <- "cell_ident"
		ofinterest_cell = subset(ofinterest, idents = cell_type)
		Idents(ofinterest_cell) <- "allele"
		markers <- FindMarkers(ofinterest_cell, ident.1 = "E4", ident.2 = "E3", min.pct = 0, logfc.threshold = 0) %>%
			rownames_to_column(var = "gene") %>% 
			left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
			dplyr::arrange(p_val_adj)
		write.csv(markers,file = paste(cell_type,"_apoe_markers_old.gsea.csv", sep = ""))
	}
}

save(apoe, file = "apoe_updated.Robj")