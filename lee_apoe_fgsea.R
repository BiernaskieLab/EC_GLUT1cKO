#Lee Apoe----------------------------------------------------------------

#SLURM--------

salloc --mem=50G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate fgsea

#fGSEA-------------------------------
R
library(fgsea)
library(tidyverse)
library(ggplot2)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/")

GO_BP <- gmtPathways("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/m5.go.bp.v2023.2.Mm.symbols.gmt")

cells.list <- c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")

for (i in 1:length(cells.list)) {
		
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/deg")
	filename <- paste0(cells.list[i],"_apoe_markers_old.gsea.csv")			
	allgene_list <- read.csv(file= filename, header=TRUE)
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/gsea")

	#grab only gene id and log2FC
	allgene_list <- allgene_list[,c(2,4)]

	allgene_list <- allgene_list[order(-allgene_list$avg_log2FC),]
	gene_list <- allgene_list$avg_log2FC
	names(gene_list) <- allgene_list$gene
	print(head(gene_list))
	print(cells.list[i])

	fgsea <- fgseaMultilevel(GO_BP, gene_list, minSize = 5, maxSize = 500, gseaParam = 0)
	saveRDS(fgsea, paste0(cells.list[i],"_apoe_fgsea.rds"))
	
	fgsea <- apply(fgsea,2,as.character)
	fgsea_df <- as.data.frame(fgsea)
	
	fgsea_df <- fgsea_df[order(fgsea_df$NES, decreasing = T),]
	write.csv(fgsea_df, paste0(cells.list[i],"_apoe_fgsea.csv"))
}
