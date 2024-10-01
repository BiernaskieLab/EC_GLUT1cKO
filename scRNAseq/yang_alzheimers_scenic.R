#Yang alzheimers scenic---------------------------------------------------------------------------------------------

#SLURM---------
salloc --mem=160G --cpus-per-task=32 --nodes=1 --ntasks=1 --time=07:30:00
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"

#Scenic --------------------------------------

source activate seurat4
R
library(Seurat)
library(SeuratObject)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human")
hippo <- readRDS("x_hippocampus.rds")
DefaultAssay(hippo) <- "RNA"

#Downsample sample too large-------------------
Idents(hippo) <- "celltype"
hippo_1 <- subset(x = hippo, idents = c("BEC_Capillary", "BEC_Venous", "BEC_Arterial", "Pericyte", "SMC"), downsample = 3000)
hippo_2 <- subset(x = hippo, idents = c("Oligodendrocyte", "Astrocyte"), downsample = 10000)
hippo_3 <- subset(x = hippo, idents = c("BEC_Capillary", "BEC_Venous", "BEC_Arterial", "Pericyte", "SMC", "Oligodendrocyte", "Astrocyte"), invert=TRUE)

hippo4 <- merge(x = hippo_3, y = list(hippo_1, hippo_2))

saveRDS(hippo4, "x_hippocampus_downsample.rds")

q()

#Scenic prep ----------------------------------------------------------------------------
source activate scenic
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
	colVars <- list(orig.ident=c("AD"="blue", "Ctrl"="red"))
	colVars$orig.ident <- colVars$orig.ident[intersect(names(colVars$orig.ident), cellInfo$orig.ident)]
	saveRDS(colVars, file="int/colVars.Rds")
	### Initialize settings
	myDatasetTitle <- "AD_hippo"
	dbDir <- "/work/biernaskie_lab/apun/scenic_database/cisTarget_databases"
	dbs <- defaultDbNames[[origin]]
	scenicOptions <- initializeScenic(org= origin, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=16)
	scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
	scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
	saveRDS(scenicOptions, file="int/scenicOptions.Rds")
}

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human")
hippo <- readRDS("x_hippocampus_downsample.rds")
DefaultAssay(hippo) <- "RNA"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")

object <- hippo
cellInfoList <- c("sample_ident", "celltype", "seurat_clusters") #info to put in cellInfo
origin <- "hgnc" #for human 

runSCENICprep(object, cellInfoList, origin)

motifAnnotations_hgnc <- motifAnnotations
#Bugfix run this to preload motifs later
save(motifAnnotations, file="motifAnnotations.Robj")
load("motifAnnotations.Robj", verbose=TRUE)
motifAnnotations_hgnc <- motifAnnotations

#Rerun
runSCENICprep(object, cellInfoList, origin)
scenicOptions <- readRDS("int/scenicOptions.Rds")
#ExprMat
exprMat <- as.matrix(object[["RNA"]]$counts)
dim(exprMat)
#saves as 1.1
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
genesKept <- readRDS("int/1.1_genesKept.Rds")
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered_log <- log2(exprMat_filtered+1)
#saves as 1.2
runCorrelation(exprMat_filtered, scenicOptions)

##Export for arboreto/grnboost2----------------------------
#ExportforArboreto function incorrect missing row names do it yourself: https://github.com/aertslab/pySCENIC/issues/67
allTFs <- getDbTfs(scenicOptions)
inputTFs <- allTFs[allTFs %in% rownames(exprMat_filtered_log)]
write.table(inputTFs, file = "1.1_inputTFs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
exprMat_filtered_t <- t(exprMat_filtered_log)
write.table(exprMat_filtered_t, file = "1.1_exprMatrix_filtered_t.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

q()

#Grnboost2----------------------------------------------------------------------------
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
source activate pyscenic
cd /work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic

python
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
import os

DATA_FOLDER="/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic"
SC_EXP_FNAME = os.path.join(DATA_FOLDER, "1.1_exprMatrix_filtered_t.tsv")
TFS_FNAME = os.path.join(DATA_FOLDER, '1.1_inputTFs.txt')
ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "grn_output.tsv")

if __name__ == '__main__':
	ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0)
	tf_names = load_tf_names(TFS_FNAME)
	adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
	adjancencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
quit()

#Convert grnboost output to R scenic version ----------------------------
conda activate scenic
R
library(SCENIC)
library(RcisTarget)
library(AUCell)
library(SCopeLoomR)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")
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
saveRDS(grnboost_linklist, "int/1.4_GENIE3_linkList.Rds")


#Rest of scenic processing R script--------------------------------------
############################################
library(Seurat)
library(SeuratObject)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")

load("motifAnnotations.Robj", verbose=TRUE)
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo <- readRDS(file= "int/cellInfo.Rds")

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human")
hippo <- readRDS("x_hippocampus_downsample.rds")
DefaultAssay(hippo) <- "RNA"

object <- hippo

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")

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

scenicOptions@fileNames$output["loomFile",] <- "output/hippo_SCENIC.loom"
#did not test temp line
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)
############################################

#SLURM job-name
#Scenic processing 
############################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --partition=cpu2023
#SBATCH --mem=185G
#SBATCH --cpus-per-task=24

#SBATCH --job-name=scenicG123_hippo
#SBATCH --output=scenicG123_hippo.out
#SBATCH --error=scenicG123_hippo.err

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
source activate scenic
R CMD BATCH scenic_hippo.R 
############################################

#Scenic analysis ----------------------------------------------------------------------------------------
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
salloc --mem=120G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=8:00:00

source activate scenic
R
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(SCopeLoomR)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")

scenicLoomPath <- "/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic/output/hippo_SCENIC.loom"
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
write.csv(regulonAUC_matrix, file="hippo_regulonAUC_matrix.csv")


###Seurat analysis---------------------------------------------------
source activate seurat4
R
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)

#Convert and add matrix---------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")
regulonAUC_matrix <- read.csv("hippo_regulonAUC_matrix.csv", row.names = 1)
#Weird bug? - got replaced by . for some reason 
colnames(regulonAUC_matrix) <- sub("\\.", "-", colnames(regulonAUC_matrix))
regulonAUC_df <- as.data.frame(regulonAUC_matrix)
AUC_assay <- CreateAssayObject(counts = regulonAUC_df)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/")
hippo <- readRDS("x_hippocampus_downsample.rds")

hippo[["AUC"]] <- AUC_assay
DefaultAssay(hippo) <- "AUC"
Idents(hippo) <- "sample_ident"
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")
saveRDS(hippo,file="x_hippocampus_downsample_AUC.rds")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic")
hippo <- readRDS(hippo,file="x_hippocampus_downsample_AUC.rds")
#DEG
DefaultAssay(hippo) <- "AUC"
Idents(object = hippo) <- "celltype"
hippo <- subset(hippo, idents = "Astrocyte/microglia", invert = TRUE)

#DEG---------------------
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic/AUC_deg")
for (cell_type in unique(hippo$celltype)) {
	Idents(object = hippo) <- "celltype"
	sub <- subset(hippo, idents = cell_type)
	DefaultAssay(sub) <- "AUC"

	Idents(sub) <- "sample_ident"
	#Do all for volcano
	markers <- FindMarkers(sub, ident.1 = "AD", ident.2 = "Ctrl", logfc.threshold = 0, min.pct = 0.1)
	write.csv(markers, paste0(cell_type, "_AD_AUC_orig_pct1AD.csv"))
}

#Volcano plot------------
cell_type_list <- c("Ependymal", "Oligodendrocyte", "Microglia", "Astrocyte")
library(ggrepel)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/scenic/AUC_deg")
for (cell_type in cell_type_list) {
	deg <- read.csv(paste0(cell_type,"_AD_AUC_orig_pct1AD.csv"))
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
	deg_label_p_val <- deg_sub %>% arrange(p_val_adj) %>% slice_head(n = 10) %>% pull(gene)
	
	deg_label <- c(deg_label_desc, deg_label_asc, deg_label_p_val)
	
	#Volcano Plot
	deg$label <- "not important"
	deg$label[deg$gene %in% deg_label] <- "important"
	print(cell_type)
	
	
	png(file = paste0(cell_type,"_AD_auc_volcano.png"), width = 18, height = 16, units = "cm", res = 900)
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
		size = 3, min.segment.length = 0, box.padding = unit(0.7, "lines"), point.padding = unit(0.2, "lines"), max.overlaps = Inf, max.time = 20, max.iter = 100000)
	print(p1)
	dev.off()	
}
