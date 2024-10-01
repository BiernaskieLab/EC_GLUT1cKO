##SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=5:00:00
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate seurat4

###Seurat--------------------------------------------------------------------------------------------------------
R
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)

library(purrr)
library(tibble)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")

#Review object
GLUTKO <- readRDS(file = "GLUTKO.rds")
save(GLUTKO,file="GLUTKO.Robj")
load("GLUTKO.Robj", verbose = TRUE)
head(GLUTKO)
head(GLUTKO@meta.data)
GLUTKO$orig.ident <- NULL

#Fix naming
Idents(GLUTKO) <- "sample"
GLUTKO <- RenameIdents(GLUTKO,"6mo female control" = "6 months female control")
GLUTKO <- RenameIdents(GLUTKO,"6mo female Glut1KO" = "6 months female Glut1KO")
GLUTKO <- RenameIdents(GLUTKO,"12mo male control" = "12 months male control")
GLUTKO <- RenameIdents(GLUTKO,"12mo male Glut1KO" = "12 months male Glut1KO")
GLUTKO <- RenameIdents(GLUTKO,"12mo female control" = "12 months female control")
GLUTKO <- RenameIdents(GLUTKO,"12mo female Glut1KO" = "12 months female Glut1KO")
GLUTKO$sample <- Idents(object = GLUTKO)

##Quality control--------------------------------------------------------------------------------------------------------
png(file = paste0("VlnPlot_QC_subset_glut.png"), width = 20, height = 17, units = "cm", res = 500)
p <- VlnPlot(GLUTKO, group.by = "sample_id", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()

png(file = paste0("Scatter_QC_subset_glut.png"), width = 35, height = 17, units = "cm", res = 500)
p1 <- FeatureScatter(GLUTKO, group.by = "sample_id", feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(GLUTKO, group.by = "sample_id", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2)
dev.off()

png(file = paste0("VlnPlot_QC_subset_glut_sample.png"), width = 20, height = 17, units = "cm", res = 500)
p <- VlnPlot(GLUTKO, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()
png(file = paste0("Scatter_QC_subset_glut_sample.png"), width = 35, height = 17, units = "cm", res = 500)
p1 <- FeatureScatter(GLUTKO, group.by = "sample", feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(GLUTKO, group.by = "sample", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2)
dev.off()

mean(GLUTKO$nFeature_RNA) + 3 * sd(GLUTKO$nFeature_RNA)
mean(GLUTKO$nCount_RNA) + 3 * sd(GLUTKO$nCount_RNA)
mean(GLUTKO$percent.mt) + 3 * sd(GLUTKO$percent.mt)

table(GLUTKO$sample)
GLUTKO <- subset(GLUTKO, subset = nFeature_RNA < 5500 & nCount_RNA < 20000 & percent.mt < 25)
table(GLUTKO$sample)

png(file = paste0("VlnPlot_QC_subset_glut_sample_QC.png"), width = 20, height = 17, units = "cm", res = 500)
p <- VlnPlot(GLUTKO, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()
png(file = paste0("Scatter_QC_subset_glut_sample-QC.png"), width = 35, height = 17, units = "cm", res = 500)
p1 <- FeatureScatter(GLUTKO, group.by = "sample", feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(GLUTKO, group.by = "sample", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2)
dev.off()

##Processing--------------------------------------------------------------------------------------------------------
GLUTKO[["RNA"]] <- split(GLUTKO[["RNA"]], f = GLUTKO$sample)
#Normalization 
GLUTKO <- NormalizeData(GLUTKO)
GLUTKO <- FindVariableFeatures(GLUTKO)
GLUTKO <- ScaleData(GLUTKO, verbose = FALSE)
GLUTKO <- RunPCA(GLUTKO, npcs = 30, verbose = FALSE)
#New seurat5 integration method
options(future.globals.maxSize = 8000 * 1024^2)
GLUTKO <- IntegrateLayers(object = GLUTKO, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
GLUTKO[["RNA"]] <- JoinLayers(GLUTKO[["RNA"]])
#Run UMAP with desired dimension
GLUTKO <- FindNeighbors(GLUTKO, reduction = "integrated.cca", dims = 1:30)
GLUTKO <- FindClusters(GLUTKO, resolution = c(0.5))
GLUTKO <- RunUMAP(GLUTKO, reduction = "integrated.cca", dims = 1:30)
#Check clusters
head(GLUTKO@meta.data)
#Check clusters through all resolutions (Change range to match min/max resolution)
for (i in 1:100) {
	num <- i*0.01
	if (paste0("RNA_snn_res.", num) %in% colnames(GLUTKO@meta.data)) {
		Idents(object = GLUTKO) <- paste0("RNA_snn_res.", num)
		png(file = paste0("Cluster_glut_m_res", num,".png"), width = 50, height = 25, units = "cm", res = 500)
		p1 <- DimPlot(GLUTKO, reduction = "umap", group.by = "sample")
		p2 <- DimPlot(GLUTKO, reduction = "umap", label = TRUE, repel = TRUE)
		p1 + p2
		#need to have print to make png work in a loop
		print(p1+p2)
		dev.off()
		png(file = paste0("Cluster_glut_m_res", num,"_split.png"), width = 50, height = 25, units = "cm", res = 500)
		p <- DimPlot(GLUTKO, reduction = "umap", split.by = "sample")
		print(p)
		dev.off()
	}
}

GLUTKO$seurat_clusters <- GLUTKO$RNA_snn_res.0.5
Idents(object = GLUTKO) <- GLUTKO$seurat_clusters


##Annotation--------------------------------------------------------------------------------------------------------

##Descriptions for genes--------------------------------------
#https://hbctraining.github.io/scRNA-seq_online/lessons/fetching_annotations.html
#Connect to AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()
#Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)
#Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
#Download the appropriate Ensembldb database
edb <- ah[[id]]
#Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")
#Select annotations of interest
annotations <- annotations %>%
    dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
write.csv(annotations,file="annotations_Mus_EnsDb.csv")
annotations <- read.csv("annotations_Mus_EnsDb.csv")


annotations <- read.csv("annotations_Mus_EnsDb.csv")
Idents(object = GLUTKO) <- "seurat_clusters"
GLUTKO.markers <- FindAllMarkers(GLUTKO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.csv(GLUTKO.markers,'GLUTKO_findallmarkers.csv')


#Check misc markers
png(file = "Feature_markers_check.png", width = 40, height = 40, units = "cm", res = 300)
FeaturePlot(GLUTKO, features = c("Rgs5", "Cnn1", "Pdgfrb", "Tagln"))
dev.off()
png(file = "Feature_NSC_markers_check.png", width = 60, height = 60, units = "cm", res = 200)
FeaturePlot(GLUTKO, features = c("S100b", "Tspan18", "Hes5", "Id4", "Sox9", "Hopx", "Thbs4"))
dev.off()
png(file = "Feature_NSC2_markers_check.png", width = 40, height = 60, units = "cm", res = 200)
FeaturePlot(GLUTKO, features = c("S100a6", "Gfap", "Slc6a11", "Aqp4", "Clu", "Slc1a3"))
dev.off()
png(file = "Dot_NSC_markers_check.png", width = 30, height = 10, units = "cm", res = 200)
DotPlot(GLUTKO, features = c("S100b", "Tspan18", "Hes5", "Id4", "Sox9", "Hopx", "Thbs4", "S100a6", "Gfap", "Slc6a11", "Aqp4", "Clu", "Slc1a3"), group.by = "cell_ident")
dev.off()
png(file = "Feature_NSC3_markers_check.png", width = 10, height = 10, units = "cm", res = 200)
FeaturePlot(GLUTKO, features = c("Cdk1"))
dev.off()
png(file = "Feature_epen_markers_check.png", width = 20, height = 30, units = "cm", res = 200)
FeaturePlot(GLUTKO, features = c("Ccdc153","Acta2","Rarres2","Odf3b", "Pifo", "Foxj1", "Aif1"))
dev.off()

png(file = "Vln_sox2_split.png", width = 25, height = 15, units = "cm", res = 200)
VlnPlot(GLUTKO, features = c("Sox2"), group.by = "cell_ident", split.by = "genotype")
dev.off()

#Examine Sox2+ cells
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
Idents(GLUTKO) <- "gender"
females <- subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
GLUTKO_6mo_fem <- subset(females, idents = "6 months")

png(file = "Vln_sox2_split_6mo_fem.png", width = 25, height = 15, units = "cm", res = 200)
VlnPlot(GLUTKO_6mo_fem, features = c("Sox2"), group.by = "cell_ident", split.by = "genotype")
dev.off()

Idents(females) <- "age"
GLUTKO_12mo_fem <- subset(females, idents = "12 months")

png(file = "Vln_sox2_split_12mo_fem.png", width = 25, height = 15, units = "cm", res = 200)
VlnPlot(GLUTKO_12mo_fem, features = c("Sox2"), group.by = "cell_ident", split.by = "genotype")
dev.off()

png(file = "Vln_sox2_split_age_fem.png", width = 17, height = 15, units = "cm", res = 200)
VlnPlot(females, features = c("Sox2"), group.by = "age", split.by = "genotype")
dev.off()

#Glial progenitors-1 vs 2 check
Idents(object = GLUTKO) <- GLUTKO$cell_ident

annotations <- read.csv("annotations_Mus_EnsDb.csv")

markers <- FindMarkers(GLUTKO, ident.1 = "Glial progenitors-1", ident.2 = "Glial progenitors-2", min.pct = 0.1, logfc.threshold = 0.25) %>%
		rownames_to_column(var = "gene") %>% 
		left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
		dplyr::arrange(p_val_adj)
  write.csv(markers,file = paste("Glut1_Glial progenitors-1_pct1_vs_gp-2_markers.total.csv", sep = ""))

#Oligodendorcytes check
Idents(object = GLUTKO) <- "cell_ident"
oligclus <- subset(GLUTKO, idents = c("Oligodendrocytes-1","Oligodendrocytes-2","Oligodendrocytes-3","Oligodendrocytes-4"))

annotations <- read.csv("annotations_Mus_EnsDb.csv")
Idents(object = oligclus) <- "cell_ident"
GLUTKO.markers <- FindAllMarkers(oligclus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
write.csv(GLUTKO.markers,'GLUTKO_findallmarkers_oligclus.csv')

#Slca plot
slca <- c("Slc2a1", "Slc2a2", "Slc2a3", "Slc2a4", "Slc2a5", "Slc2a6", "Slc2a7", "Slc2a8", "Slc2a9", "Slc2a10", "Slc2a12", "Slc2a13")
png(file = "Glut1ko_Slc2a_genes_celltype_dot_total.png", width = 18, height = 13, units = "cm", res = 900)
DotPlot(GLUTKO, features = slca, group.by = "cell_ident", dot.scale = 6, dot.min = 0.01, scale=FALSE) + RotatedAxis()
dev.off()

#Rename---------------------------------------------------------
Idents(object = GLUTKO) <- GLUTKO$seurat_clusters
GLUTKO <- RenameIdents(GLUTKO, 
	`0` = "Pericytes", 
	`1` = "Astrocytes", 
	`2` = "Pericytes",
    `3` = "Oligodendrocytes-1", 
	`4` = "Oligodendrocytes-2", 
	`5` = "Smooth muscle cells", 
	`6` = "Endothelial", 
	`7` = "Microglia", 
	`8` = "Neuroblasts-1", 
	`9` = "aNSCs", 
	`10` = "Ependymal", 
	`11` = "Pericytes", 
	`12` = "Pericytes", 
	`13` = "Oligodendrocytes-3", 
	`14` = "qNSCs", 
	`15` = "Choroid plexus", 
	`16` = "TAPs", 
	`17` = "Oligodendrocytes-4", 
	`18` = "Macrophages", 
	`19` = "Pericytes", 
	`20` = "Glial progenitors-1", 
	`21` = "Glial progenitors-2", 
	`22` = "Erythrocytes",
	`23` = "Microglia",
	`24` = "OPCs",
	`25` = "T cells",
	`26` = "Neuroblasts-2",
	`27` = "Unknown-1",
	`28` = "B cells",
	`29` = "Endothelial",
	`30` = "Unknown-2",
	`31` = "Tanycytes",
	`32` = "Endothelial")
GLUTKO$cell_ident <- Idents(object = GLUTKO)

Idents(object = GLUTKO) <- GLUTKO$seurat_clusters
GLUTKO <- RenameIdents(GLUTKO, 
	`0` = "Pericytes", 
	`1` = "Astrocytes", 
	`2` = "Pericytes",
    `3` = "Oligodendrocytes", 
	`4` = "Oligodendrocytes", 
	`5` = "Smooth muscle cells", 
	`6` = "Endothelial", 
	`7` = "Microglia", 
	`8` = "Neuroblasts-1", 
	`9` = "aNSCs", 
	`10` = "Ependymal", 
	`11` = "Pericytes", 
	`12` = "Pericytes", 
	`13` = "Oligodendrocytes", 
	`14` = "qNSCs", 
	`15` = "Choroid plexus", 
	`16` = "TAPs", 
	`17` = "Oligodendrocytes", 
	`18` = "Macrophages", 
	`19` = "Pericytes", 
	`20` = "Glial progenitors-1", 
	`21` = "Glial progenitors-2", 
	`22` = "Erythrocytes",
	`23` = "Microglia",
	`24` = "OPCs",
	`25` = "T cells",
	`26` = "Neuroblasts-2",
	`27` = "Unknown-1",
	`28` = "B cells",
	`29` = "Endothelial",
	`30` = "Unknown-2",
	`31` = "Tanycytes",
	`32` = "Endothelial")
GLUTKO$cell_ident <- Idents(object = GLUTKO)

#Overall marker plot----------------------------------------------------------

mural <- c("Pdgfrb", "Acta2", "Rgs5", "Tagln")

epen <- c("Rarres2", "Ccdc153", "Foxj1")
tan <- c("Slit2", "Col23a1", "Rax")
choroid <- c("Ttr", "Folr1", "Prlr")

nsc <- c("Id4", "Sox9", "Hopx", "Aldoc", "Clu", "Ntsr2", "Ptprz1")
gp <- c("Gria2", "Hes5", "Cspg5")
tap <- c("Top2a", "Mki67", "Birc5")
neuroblast <- c("Tubb3", "Dcx", "Dlx1")

opc <- c("Rtn1", "Nsg2", "Sox10")
odc <- c("Mbp", "Plp1", "Opalin")
astro <- c("Slc6a11", "Aqp4", "Agt", "Aldh1l1")

cd45 <- c("Ptprc")
micro_macro <- c("Csf1r", "Itgam", "Cx3cr1", "Trem2")
micro <- c("Hexb", "Tmem119", "Sall1")
macro <- c("Mrc1")
b <- c("Cd79a", "H2-Aa", "Ly6d")
t <- c("Cd3e", "Ccl5")
endo <- c("Pecam1", "Cd34", "Icam2")
eryth <- c("Hba-a1")

total_markers <- c(mural, epen, tan, choroid, nsc, gp, tap, neuroblast, opc, odc, astro, cd45, micro_macro, micro, macro, b, t, endo, eryth)
	
total_order <- c("Pericytes", "Smooth muscle cells", "Ependymal", "Tanycytes", "Choroid plexus", "qNSCs", "aNSCs", "Glial progenitors-1", "Glial progenitors-2", "TAPs", "Neuroblasts-1", "Neuroblasts-2", "OPCs", "Oligodendrocytes", "Astrocytes", "Microglia", "Macrophages", "B cells", "T cells", "Endothelial", "Erythrocytes", "Unknown-1", "Unknown-2")
#Reorganize levels 
Idents(object = GLUTKO) <- GLUTKO$cell_ident
GLUTKO$cell_ident <- factor(GLUTKO$cell_ident, levels= total_order)

png(file = "Annotations_dot_total.png", width = 43, height = 17, units = "cm", res = 900)
DotPlot(GLUTKO, features = total_markers, group.by = "cell_ident", idents = total_order, dot.scale = 6, dot.min = 0.01, scale=TRUE) + RotatedAxis()
dev.off()

#Colours
library(scales)
print(hue_pal()(23))

# [1] "#F8766D" "#EC823C" "#DD8D00" "#CA9700" "#B3A000" "#97A900" "#71B000"
# [8] "#2FB600" "#00BB4B" "#00BF76" "#00C098" "#00C0B7" "#00BDD1" "#00B7E8"
#[15] "#00AEFA" "#3DA1FF" "#8F91FF" "#BE80FF" "#DE71F9" "#F265E7" "#FE61CF"
#[22] "#FF64B3" "#FF6C92"


#UMAP FInal------------------------------------------------------------------
#Solo cell ident
Idents(object = GLUTKO) <- "cell_ident"
png(file = "Cell_ident_glut_0.5_small_legend.png", width = 34, height = 21, units = "cm", res = 700)
DimPlot(GLUTKO, reduction = "umap", label = FALSE, repel = TRUE, label.size = 4, pt.size = 0.05)
dev.off()
png(file = "Cell_ident_glut_0.5_small_label.png", width = 34, height = 21, units = "cm", res = 700)
DimPlot(GLUTKO, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.05)
dev.off()
png(file = "Cell_ident_glut_0.5_small_nolegend.png", width = 25, height = 25, units = "cm", res = 700)
DimPlot(GLUTKO, reduction = "umap", label = FALSE, repel = TRUE, label.size = 4, pt.size = 0.05) + theme(legend.position = "none")
dev.off()

#Percent plots----------------------------------------------------------

##Calculate cell type % total
#females (6 + 12)
Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_6mo = subset(females, idents = "6 months")
Idents(females_6mo) <- "cell_ident"

cell_counts_type <- table(females_6mo@meta.data$cell_ident, females_6mo@meta.data$genotype, females_6mo@meta.data$age, females_6mo@meta.data$gender)
df <- as.data.frame.table(cell_counts_type)
df$Percentage <- df$Freq / tapply(df$Freq, df$Var1, sum)[df$Var1] * 100
df <- df %>% group_by(Var2) %>% mutate(total_count = sum(Freq),Percentage = Freq/total_count*100)

Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_12mo = subset(females, idents = "12 months")
Idents(females_12mo) <- "cell_ident"

cell_counts_type <- table(females_12mo@meta.data$cell_ident, females_12mo@meta.data$genotype, females_12mo@meta.data$age, females_12mo@meta.data$gender)
df2 <- as.data.frame.table(cell_counts_type)
df2$Percentage <- df2$Freq / tapply(df2$Freq, df2$Var1, sum)[df2$Var1] * 100
df2 <- df2 %>% group_by(Var2) %>% mutate(total_count = sum(Freq),Percentage = Freq/total_count*100)

merged_df <- merge(df, df2, all = TRUE)

write.csv(merged_df , file = "Celltype_percent_total_fem.csv")

#stacked + vert
png(file = "Celltype_percent_total_stacked_fem_shuff.png", width = 30, height = 20, units = "cm", res = 700)
ggplot(merged_df, aes(x = Var2, y = Percentage, fill = Var1)) +
	theme_bw() +
	geom_bar(stat = "identity", width = 0.75) +
	scale_x_discrete(labels=c("Control", "Glut1KO")) +
	labs(y = "Percent Composition", fill = "Category") +
	facet_wrap(~ Var3, ncol = 2) +
	#Add black bar visibility  
	geom_col(position = "stack", color = "black", width = 0.75) +
	theme(legend.title = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size = 20), 
	axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 20), legend.text=element_text(size=15))+
	#Add values to bar
	geom_text(aes(label = ifelse(Percentage > 2.5, round(Percentage), "")), size = 3, position = position_stack(vjust = 0.5)) +
	facet_wrap(~ Var3, ncol = 2)
dev.off()

#females (6 + 12) - remove pericytes and SM
Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_6mo = subset(females, idents = "6 months")
Idents(females_6mo) <- "cell_ident"
females_6mo <- subset(x = females_6mo, idents = c("Smooth muscle cells", "Pericytes"), invert = TRUE)
cell_counts_type <- table(females_6mo@meta.data$cell_ident, females_6mo@meta.data$genotype, females_6mo@meta.data$age, females_6mo@meta.data$gender)
df <- as.data.frame.table(cell_counts_type)
df$Percentage <- df$Freq / tapply(df$Freq, df$Var1, sum)[df$Var1] * 100
df <- df %>% group_by(Var2) %>% mutate(total_count = sum(Freq),Percentage = Freq/total_count*100)

Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_12mo = subset(females, idents = "12 months")
Idents(females_12mo) <- "cell_ident"
females_12mo <- subset(x = females_12mo, idents = c("Smooth muscle cells", "Pericytes"), invert = TRUE)
cell_counts_type <- table(females_12mo@meta.data$cell_ident, females_12mo@meta.data$genotype, females_12mo@meta.data$age, females_12mo@meta.data$gender)
df2 <- as.data.frame.table(cell_counts_type)
df2$Percentage <- df2$Freq / tapply(df2$Freq, df2$Var1, sum)[df2$Var1] * 100
df2 <- df2 %>% group_by(Var2) %>% mutate(total_count = sum(Freq),Percentage = Freq/total_count*100)

merged_df <- merge(df, df2, all = TRUE)
write.csv(merged_df , file = "Celltype_percent_total_fem_noperi&SM.csv")

#stacked + vert
png(file = "Celltype_percent_total_stacked_fem_noperi&SM.png", width = 30, height = 20, units = "cm", res = 700)
ggplot(merged_df, aes(x = Var2, y = Percentage, fill = Var1)) +
	theme_bw() +
	geom_bar(stat = "identity", width = 0.75) +
	scale_x_discrete(labels=c("Control", "Glut1KO")) +
	labs(y = "Percent Composition", fill = "Category") +
	facet_wrap(~ Var3, ncol = 2) +
	#Add black bar visibility  
	geom_col(position = "stack", color = "black", width = 0.75) +
	theme(legend.title = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size = 20), 
	axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 20), legend.text=element_text(size=15))+
	#Add values to bar
	geom_text(aes(label = ifelse(Percentage > 2.5, round(Percentage), "")), size = 3, position = position_stack(vjust = 0.5)) +
	facet_wrap(~ Var3, ncol = 2)
dev.off()

write.csv(merged_df , file = "Celltype_percent_total_fem_noperi&SM_oligclus.csv")
png(file = "Celltype_percent_total_stacked_fem_noperi&SM_oligclus.png", width = 30, height = 20, units = "cm", res = 700)
ggplot(merged_df, aes(x = Var2, y = Percentage, fill = Var1)) +
	theme_bw() +
	geom_bar(stat = "identity", width = 0.75) +
	scale_x_discrete(labels=c("Control", "Glut1KO")) +
	labs(y = "Percent Composition", fill = "Category") +
	facet_wrap(~ Var3, ncol = 2) +
	#Add black bar visibility  
	geom_col(position = "stack", color = "black", width = 0.75) +
	theme(legend.title = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size = 20), 
	axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 20), legend.text=element_text(size=15))+
	#Add values to bar
	geom_text(aes(label = ifelse(Percentage > 2.5, round(Percentage), "")), size = 3, position = position_stack(vjust = 0.5)) +
	facet_wrap(~ Var3, ncol = 2)
dev.off()



#DEG cell_ident------------------------------------------------------------------------------------------------------------

#Total-----------------------
annotations <- read.csv("annotations_Mus_EnsDb.csv")
for(cell_type in unique(GLUTKO$cell_ident)){ 
  Idents(GLUTKO) <- "cell_ident"
  ofinterest = subset(GLUTKO, idents = cell_type)
  Idents(ofinterest) <- "genotype"
  markers <- 
    FindMarkers(ofinterest, ident.1 = "KO", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0.25) %>%
		rownames_to_column(var = "gene") %>% 
		left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
		dplyr::arrange(p_val_adj)
  write.csv(markers,file = paste("GLUT_", cell_type, "_markers.total.csv", sep = ""))
}


##Conditional DGE--------------
library(tidyr)
library(tidyverse)

annotations <- read.csv("annotations_Mus_EnsDb.csv")

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/DEG")

Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")

for(age in unique(females@meta.data$age)) {

	Idents(females) <- "age"
	ofinterest = subset(females, idents = age)

	for(cell_type in unique(ofinterest@meta.data$cell_ident)){ 
	
		Idents(ofinterest) <- "cell_ident"
		ofinterest_cell = subset(ofinterest, idents = cell_type)
		Idents(ofinterest_cell) <- "genotype"
		markers <- FindMarkers(ofinterest_cell, ident.1 = "KO", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0) %>%
			rownames_to_column(var = "gene") %>% 
			left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
			dplyr::arrange(p_val_adj)

		write.csv(markers,file = paste("Glutko_", cell_type,"_",age, "_markers_fem.total.csv", sep = ""))
	  
		markers <- markers[markers$p_val_adj<0.05,]
		markers.list.mat <- markers[,c("gene", "avg_log2FC")]
		up_score <- sum(markers.list.mat[markers.list.mat$avg_log2FC >= 0,]$avg_log2FC)
		down_score <- sum(markers.list.mat[markers.list.mat$avg_log2FC <= 0,]$avg_log2FC)
		
		GLUTKO@meta.data$enrichment_score[GLUTKO@meta.data$cell_ident == cell_type & GLUTKO@meta.data$age == age & GLUTKO@meta.data$gender == "Female"] <- up_score
		GLUTKO@meta.data$depletion_score[GLUTKO@meta.data$cell_ident == cell_type & GLUTKO@meta.data$age == age & GLUTKO@meta.data$gender == "Female"] <- down_score

	}
}

#DGE for GSEA (+ extra line just for olig subclusters)--------------
for(age in unique(females@meta.data$age)) {
	Idents(females) <- "age"
	ofinterest = subset(females, idents = age)
	for(cell_type in unique(ofinterest@meta.data$cell_ident)){ 
		#if (cell_type %in% c("Oligodendrocytes-1", "Oligodendrocytes-2", "Oligodendrocytes-3", "Oligodendrocytes-4")) {
			Idents(ofinterest) <- "cell_ident"
			ofinterest_cell = subset(ofinterest, idents = cell_type)
			Idents(ofinterest_cell) <- "genotype"
			markers <- FindMarkers(ofinterest_cell, ident.1 = "KO", ident.2 = "control", min.pct = 0, logfc.threshold = 0) %>%
				rownames_to_column(var = "gene") %>% 
				left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
				dplyr::arrange(p_val_adj)
			write.csv(markers,file = paste("Glutko_", cell_type,"_",age, "_markers_fem.gsea.csv", sep = ""))
		#}
	}
}

#Density scoring ----------------------------
#Neuroblasts (inhibitory) 12 months less than 3 cells, assign score 0 
GLUTKO@meta.data$enrichment_score[GLUTKO@meta.data$cell_ident == "Neuroblasts (inhibitory)" & GLUTKO@meta.data$age == "12 months"] <- 0
GLUTKO@meta.data$depletion_score[GLUTKO@meta.data$cell_ident == "Neuroblasts (inhibitory)" & GLUTKO@meta.data$age == "12 months"] <- 0

GLUTKO@meta.data$enrichment_score[GLUTKO@meta.data$gender == "Male"] <- NA
GLUTKO@meta.data$depletion_score[GLUTKO@meta.data$gender == "Male"] <- NA


for(age in unique(GLUTKO@meta.data$age)) {
	for(cell_type in unique(GLUTKO@meta.data$cell_ident)){ 
		up_score <- GLUTKO@meta.data$enrichment_score[GLUTKO@meta.data$cell_ident == cell_type & GLUTKO@meta.data$age == age & GLUTKO@meta.data$gender == "Female"]
		down_score <- GLUTKO@meta.data$depletion_score[GLUTKO@meta.data$cell_ident == cell_type & GLUTKO@meta.data$age == age & GLUTKO@meta.data$gender == "Female"]
		total_score <- abs(up_score) + abs(down_score)
		GLUTKO@meta.data$total_magnitude[GLUTKO@meta.data$cell_ident == cell_type & GLUTKO@meta.data$age == age & GLUTKO@meta.data$gender == "Female"] <- total_score
	}
}

#Single DGE 
Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
ofinterest = subset(females, idents = "6 months")
Idents(ofinterest) <- "cell_ident"
ofinterest_cell = subset(ofinterest, idents = "Ependymal")
Idents(ofinterest_cell) <- "genotype"
markers <- FindMarkers(ofinterest_cell, ident.1 = "KO", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0) %>%
	rownames_to_column(var = "gene") %>% 
	left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name")) %>% 
	dplyr::arrange(p_val_adj)
write.csv(markers,file = paste("Glutko_Ependymal_6_months_markers_fem.total.csv", sep = ""))



##Stacked enrichment+depletion scores------------------------------

diff_table <- GLUTKO@meta.data %>%
  select(cell_ident, sample, age, gender, genotype, enrichment_score, depletion_score) %>%
  distinct() %>%
  arrange(cell_ident) %>%
  rownames_to_column(var = "row_id") %>%
  select(-row_id) %>%
  mutate(
    net_score_clus = enrichment_score + depletion_score,
    total_magnitude = abs(enrichment_score) + abs(depletion_score))

write.csv(diff_table,file = "DEG_total_celltype.csv")

diff_table <- read.csv(file="DEG_total_celltype.csv")
diff_table <- diff_table[diff_table$gender != "Male", ]
#Only need to keep one row, duplicated value
diff_table <- diff_table[diff_table$genotype != "control", ]

#6 months---
diff_table_6_months <- diff_table[diff_table$age == "6 months", ]
diff_table_2_6mo <- diff_table_6_months %>% select(cell_ident, enrichment_score, depletion_score)

#Remove 0 scores
diff_table_2_6mo <- diff_table_2_6mo[diff_table_2_6mo$enrichment_score != 0 & diff_table_2_6mo$depletion_score != 0, ]
#Remove Pericyte/Smooth muscle
diff_table_2_6mo <- diff_table_2_6mo[diff_table_2_6mo$cell_ident != "Smooth muscle cells", ]
diff_table_2_6mo <- diff_table_2_6mo[diff_table_2_6mo$cell_ident != "Pericytes", ] 

melted_diff_table_6 <- reshape2::melt(diff_table_2_6mo, id.vars = "cell_ident")
diff_table_levels_6 <- diff_table_6_months %>% group_by(cell_ident) %>% mutate(total_magnitude = abs(enrichment_score) + abs(depletion_score)) %>% arrange(desc(total_magnitude))

#Create a stacked bar graph
png(file = "DEG_score_cell_sorted_6_mo_fem_thin.png", width = 20, height = 14, units = "cm", res = 700)
ggplot(melted_diff_table_6 , aes(x = factor(cell_ident, levels = diff_table_levels_6$cell_ident), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  ylim(-110,50) +
  scale_fill_manual(values = c("enrichment_score" = "blue", "depletion_score" = "red"),
	name = "Score Type",  # Set the legend title
    labels = c("Enrichment", "Depletion")) +
  labs(title = "6 months Enrichment and Depletion Scores",
       x = "Cell Type",
       y = "Perturbation score") +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size = 14, angle = 45, vjust = 1, hjust = 1), , axis.title.x = element_text(face="bold", size = 14), axis.text.y = element_text(face="bold", size = 14), axis.title.y = element_text(face="bold", size = 14), plot.title = element_text(face = "bold", size=14),
  plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
dev.off()

#create the final plot
Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_6mo = subset(females, idents = "6 months")
#Remove peri & smooth
Idents(females_6mo) <- "cell_ident"
females_6mo <- subset(x = females_6mo, idents = c("Smooth muscle cells", "Pericytes"), invert = TRUE)
library(Nebulosa)
png(file = "Nebulosa_Glut_fem_total_DEG_6mo.jpeg", width = 22, height = 20, units = "cm", res = 700)
#need limits argument or legend scale and colour fill was off by a bit for both
plot_density(females_6mo, "total_magnitude", pal = "inferno", limits = c(0, 0.05))
dev.off()

#12 months---
diff_table_12_months <- diff_table[diff_table$age == "12 months", ]
diff_table_2_12mo <- diff_table_12_months %>% select(cell_ident, enrichment_score, depletion_score)

#Remove 0 scores
diff_table_2_12mo <- diff_table_2_12mo[diff_table_2_12mo$enrichment_score != 0 & diff_table_2_12mo$depletion_score != 0, ]
#Remove Pericyte/Smooth muscle
diff_table_2_12mo <- diff_table_2_12mo[diff_table_2_12mo$cell_ident != "Smooth muscle cells", ]
diff_table_2_12mo <- diff_table_2_12mo[diff_table_2_12mo$cell_ident != "Pericytes", ] 

melted_diff_table_12 <- reshape2::melt(diff_table_2_12mo, id.vars = "cell_ident")
diff_table_levels_12 <- diff_table_12_months %>% group_by(cell_ident) %>% mutate(total_magnitude = abs(enrichment_score) + abs(depletion_score)) %>% arrange(desc(total_magnitude))

#Create a stacked bar graph
png(file = "DEG_score_cell_sorted_12_mo_fem_thin.png", width = 20, height = 14, units = "cm", res = 700)
ggplot(melted_diff_table_12 , aes(x = factor(cell_ident, levels = diff_table_levels_12$cell_ident), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  ylim(-110,50) +
  scale_fill_manual(values = c("enrichment_score" = "blue", "depletion_score" = "red"),
	name = "Score Type",  # Set the legend title
    labels = c("Enrichment", "Depletion")) +
  labs(title = "12 months Enrichment and Depletion Scores",
       x = "Cell Type",
       y = "Perturbation score") +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size = 14, angle = 45, vjust = 1, hjust = 1), , axis.title.x = element_text(face="bold", size = 14), 
	axis.text.y = element_text(face="bold", size = 14), axis.title.y = element_text(face="bold", size = 14), plot.title = element_text(face = "bold", size=14),
	plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
dev.off()

#create the final plot
Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
females_12mo = subset(females, idents = "12 months")
#Remove peri & smooth
Idents(females_12mo) <- "cell_ident"
females_12mo <- subset(x = females_12mo, idents = c("Smooth muscle cells", "Pericytes"), invert = TRUE)
library(Nebulosa)
png(file = "Nebulosa_Glut_fem_total_DEG_12mo.jpeg", width = 22, height = 20, units = "cm", res = 700)
#need limits argument or legend scale and colour fill was off by a bit for both
plot_density(females_12mo, "total_magnitude", pal = "inferno", limits = c(0, 0.05))
dev.off()


#Ependymal DEG dotplot 6mo--------------------------------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/DEG/cell")

Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")
Idents(females) <- "age"
ofinterest = subset(females, idents = "6 months")
Idents(ofinterest) <- "cell_ident"
ofinterest_cell = subset(ofinterest, idents = "Ependymal")

cilia_sig <- c("Daw1","Spag1","Myl9","Bicc1","Myo1e","Crocc","Stpg1","Frmd6","Synpo2")

#to reorder ggplot ctrl top
ofinterest_cell$genotype <- factor(ofinterest_cell$genotype , levels = c("control", "KO"))

png(file = "Ependymal_cilia_dot_vert.png", width = 10.25, height = 11, units = "cm", res = 900)
DotPlot(ofinterest_cell, features = cilia_sig, group.by = "genotype", dot.scale = 6, dot.min = 0.01, scale=FALSE) + 
	RotatedAxis() +
	#Care order
	scale_y_discrete(labels = c("Control", "Glut1KO"))  +
	coord_flip()
dev.off()

#Volcano Plots DEG---------------------------------------------------
library(ggrepel)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/DEG/cell/fem_deg")
cell_type_list <- c("Glial progenitors-2_6 months", "Glial progenitors-2_12 months", "OPCs_6 months", "OPCs_12 months")

cell_type_list <- c("Ependymal_6 months", "Ependymal_12 months")
names(cell_type_list) <- c("Ependymal 6 months", "Ependymal 12 months")

cell_type_list <- c("Astrocytes_6 months", "Astrocytes_12 months")
names(cell_type_list) <- c("Astrocytes 6 months", "Astrocytes 12 months")

cell_type_list <- c("Microglia_6 months", "Microglia_12 months")
names(cell_type_list) <- c("Microglia 6 months", "Microglia 12 months")

for (cell_type_name in names(cell_type_list)) {
	cell_type <- cell_type_list[cell_type_name]
	deg <- read.csv(paste0("Glutko_",cell_type,"_markers_fem.total.csv"))

	#Threshold
	deg_sub <- deg[deg$p_val_adj<0.05, ]
	
	#Get top to label
	deg_label_desc <- deg_sub %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% pull(gene)
	deg_label_asc <- deg_sub %>% arrange(avg_log2FC) %>% slice_head(n = 10) %>% pull(gene)
	deg_label_p_val <- deg_sub %>% arrange(p_val_adj) %>% slice_head(n = 10) %>% pull(gene)
	
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
	  ggtitle(cell_type_name) +
	  geom_point(aes(color = ifelse(p_val_adj < 0.05, "<0.05", "Not Sig")), alpha = 0.6) +
	  geom_label_repel(data=deg[deg$label == "important",], 
		min.segment.length = 0, box.padding = unit(0.2, "lines"), point.padding = unit(0.3, "lines"), max.overlaps = Inf, force = 6, max.iter = 500000, max.time = 5) +
		xlim(-8,6)
	print(p1)
	dev.off()	
}

