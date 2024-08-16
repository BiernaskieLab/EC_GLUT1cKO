#Yang human AD---------------------------------------------------------------------------------------------------------
#SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=09:30:00
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate seurat4

#Seurat--------------------------------
R

library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)
library(tidyverse)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human")

hippo <- readRDS("x_hippocampus.rds")
DefaultAssay(hippo) <- "RNA"


#UMAP-----------------------------
Idents(object = hippo) <- "celltype"
png(file = "Cell_ident_hippo_0.5_small.png", width = 25, height = 17, units = "cm", res = 700)
DimPlot(hippo, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.8)
dev.off()
Idents(object = hippo) <- "celltype"
png(file = "Cell_ident_hippo_0.5_small_nolegend.png", width = 17, height = 17, units = "cm", res = 700)
DimPlot(hippo, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.8) + theme(legend.position = "none")
dev.off()
Idents(object = hippo) <- "celltype"
png(file = "Cell_ident_hippo_0.5_small_nolegend_nolabel.png", width = 17, height = 17, units = "cm", res = 700)
DimPlot(hippo, reduction = "umap", label = FALSE, pt.size = 0.8) + theme(legend.position = "none")
dev.off()

#Findallmarkers-----------------------------
Idents(object = hippo) <- "celltype"
hippo.markers <- FindAllMarkers(hippo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(hippo.markers,'hippo_AD_findallmarkers_celltype.csv')

Idents(object = hippo) <- "seurat_clusters"
hippo.markers <- FindAllMarkers(hippo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(hippo.markers,'hippo_AD_findallmarkers_clusters.csv')

#Ependymal markers-----------------------------
png(file = "VlnPlot_ependyal_markers_hippo_human.png", width = 20, height = 10, units = "cm", res = 300)
VlnPlot(hippo, features = c("DNAH12", "AK9"), group.by = "celltype")
dev.off()

png(file = "DotPlot_ependyal_markers_hippo_human.png", width = 15, height = 10, units = "cm", res = 300)
DotPlot(hippo, features = c("DNAH12", "FOXJ1", "PIFO", "DYNLRB2", "AQP4", "GFAP"), group.by = "celltype") + RotatedAxis()
dev.off()

png(file = "DotPlot_ependyal_markers_hippo_human_legend.png", width = 15, height = 15, units = "cm", res = 300)
DotPlot(hippo, features = c("DNAH12", "FOXJ1", "PIFO", "DYNLRB2", "AQP4", "GFAP"), group.by = "celltype") + RotatedAxis()
dev.off()

png(file = "DotPlot_ependyal_markers2_hippo_human_legend.png", width = 17, height = 15, units = "cm", res = 300)
DotPlot(hippo, features = c("DNAH12", "FOXJ1", "PIFO", "DYNLRB2", "RARRES2", "CCDC153", "ACTA2", "VIM", "AQP4", "GFAP"), group.by = "celltype") + RotatedAxis()
dev.off()

png(file = "DotPlot_ependyal_markers3_hippo_human_legend.png", width = 17, height = 15, units = "cm", res = 300)
DotPlot(hippo, features = c("DNAH12", "FOXJ1", "PIFO", "DYNLRB2", "CFAP126", "ODF3B", "AQP4", "GFAP"), group.by = "celltype") + RotatedAxis()
dev.off()


#DEGs for hippo-----------------------------------------------
library(ggVennDiagram)
library(ggplot2)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/DEG")

for(cell_type in unique(hippo@meta.data$celltype)){ 
	if (cell_type %in% c("Ependymal", "Microglia", "Oligodendrocyte", "Astrocyte")) {
		print(cell_type)
		Idents(hippo) <- "celltype"
		ofinterest_cell = subset(hippo, idents = cell_type)
		Idents(ofinterest_cell) <- "sample_ident"
		markers <- FindMarkers(ofinterest_cell, ident.1 = "AD", ident.2 = "Ctrl", min.pct = 0.1, logfc.threshold = 0)
		write.csv(markers,file = paste("Hippo_human_AD_",cell_type,"_markers.total.csv", sep = ""))
	}
}

#Volcano------------------------------------------
cell_type_list <- c("Ependymal")
library(ggrepel)
for (cell_type in cell_type_list) {
	
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/DEG")
	deg <- read.csv(paste0("Hippo_human_AD_",cell_type,"_markers.total.csv"))
	colnames(deg)[1] <- "gene"
	#Remove (num g) part
	deg$gene <- gsub("\\s*\\(\\d+g\\)\\s*", "", deg$gene)
	#Shorten extended
	deg$gene <- gsub("-extended", "-ext", deg$gene)
	
	#Threshold
	deg_sub <- deg[deg$p_val_adj<0.05 & (abs(deg$avg_log2FC) > 0.1), ]
	
	#Get top to label
	deg_label_desc <- deg_sub %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 15) %>% pull(gene)
	deg_label_asc <- deg_sub %>% arrange(avg_log2FC) %>% slice_head(n = 15) %>% pull(gene)
	deg_label_p_val <- deg_sub %>% arrange(p_val_adj) %>% slice_head(n = 15) %>% pull(gene)
	
	deg_label <- c(deg_label_desc, deg_label_asc, deg_label_p_val)
	
	#Volcano Plot
	deg$label <- "not important"
	deg$label[deg$gene %in% deg_label] <- "important"
	print(cell_type)
	
	
	png(file = paste0(cell_type,"_deg_volcano.png"), width = 18, height = 16, units = "cm", res = 900)
	p1 <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), label=gene)) +
	  scale_color_manual(values = c("<0.05" = "red", "Not Sig" = "black")) +
	  labs(x = "Log2(Fold Change)", y = "-log10(adj P-value)") +
	  theme_bw() +
	  guides(color = guide_legend(title = "Significance")) +
	  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
	  geom_vline(xintercept=c(0.25,-0.25), linetype="dashed", color = "red") +
	  ggtitle(cell_type) +
	  geom_point(aes(color = ifelse(p_val_adj < 0.05, "<0.05", "Not Sig")), alpha = 0.6) +
	  geom_label_repel(data=deg[deg$label == "important",], 
		min.segment.length = 0, box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), max.overlaps = Inf)
	print(p1)
	dev.off()	
}

#DGE for gsea-----------------------------
for(cell_type in unique(hippo@meta.data$celltype)){ 
	if (cell_type %in% c("Ependymal", "Microglia", "Oligodendrocyte", "Astrocyte")) {
		print(cell_type)
		Idents(hippo) <- "celltype"
		ofinterest_cell = subset(hippo, idents = cell_type)
		Idents(ofinterest_cell) <- "sample_ident"
		markers <- FindMarkers(ofinterest_cell, ident.1 = "AD", ident.2 = "Ctrl", min.pct = 0, logfc.threshold = 0)
		write.csv(markers,file = paste("Hippo_human_AD_",cell_type,"_markers.gsea.csv", sep = ""))
	}
}








