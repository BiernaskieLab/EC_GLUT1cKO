#SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"

source activate scenic

#SCENIC-------------------------------------------------------------------------------------------------------------------------------

#Prep----------------------------------------------------------------------------
R
library(Seurat)
library(SeuratObject)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)

runSCENICprep <- function (object, cellInfoList, orgin) {
	cellInfo = object@meta.data[cellInfoList]
	dir.create("int")
	saveRDS(cellInfo, file= "int/cellInfo.Rds")
	#Colors
	colVars <- list(orig.ident=c("Glut1ko"="blue", "Control"="red"))
	colVars$orig.ident <- colVars$orig.ident[intersect(names(colVars$orig.ident), cellInfo$orig.ident)]
	saveRDS(colVars, file="int/colVars.Rds")
	### Initialize settings
	myDatasetTitle <- "GLUTKO"
	dbDir <- "/work/biernaskie_lab/apun/scenic_database/cisTarget_databases"
	dbs <- defaultDbNames[[origin]]
	scenicOptions <- initializeScenic(org= origin, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=16)
	scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
	scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
	saveRDS(scenicOptions, file="int/scenicOptions.Rds")
}

#setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
#load("GLUTKO.Robj", verbose = TRUE)
load("females_glutko.Robj", verbose = TRUE)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic")

object <- females
cellInfoList <- c("age", "gender", "genotype", "cell_ident", "seurat_clusters") #info to put in cellInfo
origin <- "mgi" #for mouse 

runSCENICprep(object, cellInfoList, origin)

#Bugfix run this to fix vairable naming
motifAnnotations_mgi <- motifAnnotations
#Bugfix save preload motifs for later
save(motifAnnotations, file="motifAnnotations.Robj")
load("motifAnnotations.Robj", verbose=TRUE)
motifAnnotations_mgi <- motifAnnotations

#Rerun after bug halt
runSCENICprep(object, cellInfoList, origin)
scenicOptions <- readRDS("int/scenicOptions.Rds")
#ExprMat
exprMat <- as.matrix(object[["RNA"]]$counts)
dim(exprMat)
#saves as 1.1
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered_log <- log2(exprMat_filtered+1)
#saves as 1.2
runCorrelation(exprMat_filtered, scenicOptions)

##Export for arboreto/grnboost2
#ExportforArboreto function incorrect missing row names do it yourself: https://github.com/aertslab/pySCENIC/issues/67
allTFs <- getDbTfs(scenicOptions)
inputTFs <- allTFs[allTFs %in% rownames(exprMat_filtered_log)]
write.table(inputTFs, file = "1.1_inputTFs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
exprMat_filtered_t <- t(exprMat_filtered_log)
write.table(exprMat_filtered_t, file = "1.1_exprMatrix_filtered_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

q()

#Grnboost2----------------------------------------- 	
#Python3 script---------------------------------------------------------------------------- 	
####################################
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
import os

DATA_FOLDER="/work/biernaskie_lab/apun/nilesh_glut1/scenic"
SC_EXP_FNAME = os.path.join(DATA_FOLDER, "1.1_exprMatrix_filtered_t.tsv")
TFS_FNAME = os.path.join(DATA_FOLDER, '1.1_inputTFs.txt')
ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "grn_output.tsv")

if __name__ == '__main__':
	ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0)
	tf_names = load_tf_names(TFS_FNAME)
	adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)

	adjancencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
####################################

#SLURM job Grnboost----------------------------------------------------------------------------
####################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --partition=cpu2023
#SBATCH --mem=120G
#SBATCH --cpus-per-task=18

#SBATCH --job-name=grnboost_glut
#SBATCH --output=grnboost_glut.out
#SBATCH --error=grnboost_glut.err

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
source activate pyscenic
cd /work/biernaskie_lab/apun/nilesh_glut1/scenic
python3 grnboost_glut.py 
####################################


#Reconvert grnboost pyscenic output back to scenic------------------------------------------------------------
conda activate scenic
R
library(SCENIC)
library(RcisTarget)
library(AUCell)
library(SCopeLoomR)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic")
# Importing GRNBoost output https://github.com/aertslab/SCENIC/discussions/262

grn_output <- read.table(file="grn_output.tsv", sep = '\t', header = TRUE)
colnames(grn_output) <- c("TF", "Target", "Importance")
write.table(grn_output, file="grn_output.tsv", sep = '\t', row.names = FALSE, quote = FALSE)

grnboost_linklist <- importArboreto("grn_output.tsv")

#Check this part make sure accurate
grnboost_linklist <- grnboost_linklist[c(1,2,4)]
names(grnboost_linklist) <- c("TF", "Target", "weight")
#Confirm correct name
scenicOptions <- readRDS("int/scenicOptions.Rds")
getIntName(scenicOptions,"genie3ll")

#move to int/
saveRDS(grnboost_linklist, "1.4_GENIE3_linkList.Rds")


#Scenic processing steps script---------------------------------------------------------------------------------
####################################
library(Seurat)
library(SeuratObject)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic")

load("motifAnnotations.Robj", verbose=TRUE)
motifAnnotations_mgi <- motifAnnotations
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo <- readRDS(file= "int/cellInfo.Rds")

load("females_glutko.Robj", verbose = TRUE)
object <- females

genesKept <- readRDS("int/1.1_genesKept.Rds")
exprMat <- as.matrix(object[["RNA"]]$counts)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered_log <- log2(exprMat_filtered+1) 
scenicOptions@settings$nCores <- 24
#runGenie3(exprMat_filtered_log, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
##Reduce parallelization for SCENIC_3 https://github.com/aertslab/SCENIC/issues/50 (do so 2nd attempt)
scenicOptions@settings$nCores <- 4
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
#scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions@fileNames$output["loomFile",] <- "output/fem_glutko_SCENIC.loom"

saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)
####################################

#SLURM job scenic------------------------------------------------- 
####################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --partition=cpu2023
#SBATCH --mem=185G
#SBATCH --cpus-per-task=24

#SBATCH --job-name=scenicG123_glutko
#SBATCH --output=scenicG123_glutko.out
#SBATCH --error=scenicG123_glutko.err

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
source activate scenic
R CMD BATCH scenicG123_glutko.R 
####################################



#Scenic analysis ----------------------------------------------------------------------------------------
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
salloc --mem=40G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=8:00:00

source activate scenic
R
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(SCopeLoomR)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic/")

#Open scenic output
scenicLoomPath <- "/work/biernaskie_lab/apun/nilesh_glut1/scenic/output/fem_glutko_SCENIC.loom"
#Care need surprising amount of mem will mem out
loom <- open_loom(scenicLoomPath, mode = "r+")
    # Read information from loom file:
    exprMat <- get_dgem(loom)
        exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
    regulons_incidMat <- get_regulons(loom, column.attr.name="MotifRegulons")
        regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom)
    regulonAucThresholds <- get_regulon_thresholds(loom)
    embeddings <- get_embeddings(loom)
    #cellClusters <- get_clusterings(loom)
	cellInfo <- get_cell_annotation(loom)
close_loom(loom)

#Adjust aertslab function to return a subsetable matrix instead of nonadjustable AUCellResults (SummarizedExperiments)
get_regulons_AUC1 <- function(
  loom,
  column.attr.name="MotifRegulonsAUC",
  rows="regulons",
  columns="cells"
) {
  if(!column.attr.name %in% names(loom[["col_attrs"]]))
  {
    msg <- paste("The attribute '", column.attr.name, "' is not available in this loom file.", sep='')
    possible_values <- grep("egulon", names(x = loom[["col_attrs"]]), value=T)
    if(length(x = possible_values)>0) {
      msg <- c(
        msg, 
        " Possible values include: ",
        paste(possible_values, collapse=", "),
        paste(". Try setting the 'column.attr.name' argument to one of these values (i.e., get_regulons_AUC(loom, column.attr.name='", possible_values[1], "')).",sep="")
      )
    }
    if(length(x = possible_values) == 0) {
      msg <- c(
        msg, 
        " The loom doesn't contain regulon information."
      )
    }
    stop(msg)
  }
  
  mtx <- loom[["col_attrs"]][[column.attr.name]][]
  rownames(x = mtx) <- get_cell_ids(loom = loom)
  mtx <- t(x = mtx)
  names(x = dimnames(x = mtx)) <- c(rows, columns)
  return(mtx)
}

regulonAUC_matrix <- get_regulons_AUC1(loom)
write.csv(regulonAUC_matrix, file="glutko_regulonAUC_matrix.csv")

###Seurat ---------------------------------------------------
source activate seurat4
R
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)

#Convert matrix-----------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic/")
regulonAUC_matrix <- read.csv("glutko_regulonAUC_matrix.csv", row.names = 1)
#Weird bug? - got replaced by . for some reason 
colnames(regulonAUC_matrix) <- sub("\\.", "-", colnames(regulonAUC_matrix))
regulonAUC_df <- as.data.frame(regulonAUC_matrix)
AUC_assay <- CreateAssayObject(counts = regulonAUC_df)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic")
load("females_glutko.Robj", verbose = TRUE)

females[["AUC"]] <- AUC_assay
DefaultAssay(females) <- "AUC"
setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic")
save(females,file="fem_glutko_scenic_seurat.Robj")

load("fem_glutko_scenic_seurat.Robj", verbose=TRUE)


#Analysis-------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic/deg")

#FindAllMarkers-------------
DefaultAssay(females) <- "AUC"
Idents(object = females) <- "cell_ident"
rtx_m.markers <- FindAllMarkers(females, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(rtx_m.markers,'fem_glutko_findallmarkers.csv')

#DEG all-------------
DefaultAssay(females) <- "AUC"
for (cell_type in unique(females$cell_ident)) {
	for(age in unique(females$age)) {

		Idents(object = females) <- "age"
		females_age <- subset(females, idents = age)

		Idents(object = females_age) <- "cell_ident"
		females_age_cell <- subset(females_age, idents = cell_type)
		DefaultAssay(females_age_cell) <- "AUC"

		Idents(females_age_cell) <- "genotype"
		#Do all for volcano
		markers <- FindMarkers(females_age_cell, ident.1 = "KO", ident.2 = "control", logfc.threshold = 0, min.pct = 0.1)
		write.csv(markers, paste0("fem_", age, "_",cell_type, "_AUC_geno_pct1KO.csv"))
	}
}

#Individual DEG
age = "6 months"
Idents(object = females) <- "age"
females_age <- subset(females, idents = age)

cell_type = "Neuroblasts-2"
Idents(object = females_age) <- "cell_ident"
females_age_cell <- subset(females_age, idents = cell_type)
DefaultAssay(females_age_cell) <- "AUC"

Idents(females_age_cell) <- "genotype"
#Do all for volcano
markers <- FindMarkers(females_age_cell, ident.1 = "KO", ident.2 = "control", logfc.threshold = 0, min.pct = 0.1)
write.csv(markers, paste0("fem_", age, "_",cell_type, "_AUC_geno_pct1KO.csv"))

#Volcano plot-------------
cell_type_list <- c("fem_6 monthsqNSCs","fem_6 monthsaNSCs","fem_6 monthsGlial progenitors-1","fem_6 monthsGlial progenitors-2","fem_6 monthsTAPs")
names(cell_type_list) <- c("qNSCs", "aNSCs", "Glial progenitors-1", "Glial progenitors-2","TAPs")

cell_type_list <- c("fem_6 monthsOPCs", "fem_12 monthsOPCs")
names(cell_type_list) <- c("6 months OPCs", "12 months OPCs")

cell_type_list <- c("fem_6 monthsMicroglia", "fem_12 monthsMicroglia")
names(cell_type_list) <- c("6 months Microglia", "12 months Microglia")

cell_type_list <- c("fem_6 monthsNeuroblasts-1", "fem_12 monthsNeuroblasts-1", "fem_6 monthsAstrocytes", "fem_12 monthsAstrocytes",  "fem_6 monthsOligodendrocytes", "fem_12 monthsOligodendrocytes")
names(cell_type_list) <- c("6 months Neuroblasts-1", "12 months Neuroblasts-1", "6 months Astrocytes", "12 months Astrocytes", "6 months Oligodendrocytes", "12 months Oligodendrocytes")

cell_type_list <- c("fem_6 months_Oligodendrocytes-1", "fem_12 months_Oligodendrocytes-1", "fem_6 months_Oligodendrocytes-2", "fem_12 months_Oligodendrocytes-2", "fem_6 months_Oligodendrocytes-3", "fem_12 months_Oligodendrocytes-3", "fem_6 months_Oligodendrocytes-4", "fem_12 months_Oligodendrocytes-4")
names(cell_type_list) <- c("6 months Oligodendrocytes-1", "12 months Oligodendrocytes-1", "6 months Oligodendrocytes-2", "12 months Oligodendrocytes-2", "6 months Oligodendrocytes-3", "12 months Oligodendrocytes-3", "6 months Oligodendrocytes-4", "12 months Oligodendrocytes-4")

cell_type_list <- c("fem_6 monthsEpendymal", "fem_12 monthsEpendymal")
names(cell_type_list) <- c("6 months Ependymal", "12 months Ependymal")

library(ggrepel)
for (cell_type_name in names(cell_type_list)) {
	
	cell_type <- cell_type_list[cell_type_name]
	
	deg <- read.csv(paste0(cell_type,"_AUC_geno_pct1KO.csv"))
	colnames(deg)[1] <- "gene"
	#Remove (num g) part
	deg$gene <- gsub("\\s*\\(\\d+g\\)\\s*", "", deg$gene)
	#Shorten extended
	deg$gene <- gsub("-extended", "-ext", deg$gene)
	
	#Threshold
	#& (abs(deg$avg_log2FC) > 0.1)
	deg_sub <- deg[deg$p_val_adj<0.05, ]
	
	#Get top to label, choose 
	deg_label_desc <- deg_sub %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 7) %>% pull(gene)
	deg_label_asc <- deg_sub %>% arrange(avg_log2FC) %>% slice_head(n = 7) %>% pull(gene)
	deg_label_p_val <- deg_sub %>% arrange(p_val_adj) %>% slice_head(n = 5) %>% pull(gene)
	
	deg_label <- c(deg_label_desc, deg_label_asc, deg_label_p_val)
	
	#Volcano Plot
	deg$label <- "not important"
	deg$label[deg$gene %in% deg_label] <- "important"
	print(cell_type)
	print(cell_type_name)
	
	png(file = paste0(cell_type_name,"_auc_volcano.png"), width = 18, height = 16, units = "cm", res = 900)
	p1 <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), label=gene)) +
	  scale_color_manual(values = c("<0.05" = "red", "Not Sig" = "black")) +
	  labs(x = "Log2(Fold Change)", y = "-log10(adj P-value)") +
	  theme_bw() +
	  guides(color = guide_legend(title = "Significance")) +
	  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
	  geom_vline(xintercept=c(0.25,-0.25), linetype="dashed", color = "red") +
	  ggtitle(paste0(cell_type_name)) +
	  geom_point(aes(color = ifelse(p_val_adj < 0.05, "<0.05", "Not Sig")), alpha = 0.6) +
	  geom_label_repel(data=deg[deg$label == "important",], 
		min.segment.length = 0, box.padding = unit(0.5, "lines"), point.padding = unit(0.3, "lines"), max.overlaps = Inf)
	print(p1)
	dev.off()	
}

#NSC UMAP
DefaultAssay(females) <- "AUC"
females[["AUC"]] <- split(females[["AUC"]], f = females$sample_id)

#Normalization (Consider sctransform alternative:https://satijalab.org/seurat/articles/sctransform_vignette.html and bottom of https://satijalab.org/seurat/articles/integration_introduction.html and https://github.com/satijalab/seurat/issues/672)
females <- NormalizeData(females)
females <- FindVariableFeatures(females)
females <- ScaleData(females, verbose = FALSE)
females <- RunPCA(females, npcs = 30, verbose = FALSE)

#Approx number of PCs to select
png(file = "Elbow_ndim_glutko_fem_scenic.png", width = 35, height = 25, units = "cm", res = 500)
ElbowPlot(females, ndims = 30, reduction = "pca")
dev.off()

#New seurat5 integration method
options(future.globals.maxSize = 8000 * 1024^2)
females <- IntegrateLayers(object = females, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# re-join layers after integration
females[["AUC"]] <- JoinLayers(females[["AUC"]])
#Run UMAP with desired dimension
females <- RunUMAP(females, reduction = "integrated.cca", dims = 1:30)

Idents(object = females) <- "cell_ident"
females_NSC <- subset(females, idents = c("TAPs", "Neuroblasts-1", "Neuroblasts-2","aNSCs", "qNSCs"))

females_NSC <- RunUMAP(females_NSC, reduction = "integrated.cca", dims = 1:30)

Idents(object = females_NSC) <- "cell_ident"
png(file = "Scenic_umap_glut_fem_NSC_cell_ident_small.png", width = 12, height = 9, units = "cm", res = 700)
DimPlot(females_NSC, reduction = "umap", label = FALSE, pt.size = 0.15, order = c("qNSCs", "aNSCs", "Neuroblasts-2", "Neuroblasts-1", "TAPs"), cols = c("#00C0B7", "#71B000", "#DE71F9", "#2FB600", "#00BF76"))
dev.off()

Idents(object = females_NSC) <- "genotype"
png(file = "Scenic_umap_glut_fem_NSC_genotype_small.png", width = 10, height = 9, units = "cm", res = 700)
DimPlot(females_NSC, reduction = "umap", label = FALSE, pt.size = 0.15)
dev.off()

Idents(object = females_NSC) <- "cell_ident"
png(file = "Scenic_umap_glut_fem_NSC_cell_ident_small_noleg.png", width = 9, height = 9, units = "cm", res = 700)
DimPlot(females_NSC, reduction = "umap", label = FALSE, pt.size = 0.15, order = c("qNSCs", "aNSCs", "Neuroblasts-2", "Neuroblasts-1", "TAPs"), cols = c("#00C0B7", "#71B000", "#DE71F9", "#2FB600", "#00BF76")) + theme(legend.position = "none")
dev.off()

Idents(object = females_NSC) <- "genotype"
png(file = "Scenic_umap_glut_fem_NSC_genotype_small_noleg.png", width = 9, height = 9, units = "cm", res = 700)
DimPlot(females_NSC, reduction = "umap", label = FALSE, pt.size = 0.15) + theme(legend.position = "none")
dev.off()

setwd("/work/biernaskie_lab/apun/nilesh_glut1/scenic/")
save(females_NSC, file = "fem_glutko_NSC_scenic_seurat_newumap.Robj")
load(file = "fem_glutko_NSC_scenic_seurat_newumap.Robj")