#Yang human AD---------------------------------------------------------------------------------------------------------
#SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=09:30:00
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate fgsea


#fGSEA----------------------------------------------------------------------------------------
R
library(fgsea)
library(tidyverse)
library(ggplot2)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/DEG")

GO_BP <- gmtPathways("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/c5.go.bp.v2023.2.Hs.symbols.gmt")

cells.list <- c("Ependymal", "Microglia", "Oligodendrocyte", "Astrocyte")

for (i in 1:length(cells.list)) {
		
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/DEG")
	filename <- paste0("Hippo_human_AD_",cells.list[i],"_markers.gsea.csv")		
	allgene_list <- read.csv(file= filename, header=TRUE)
	
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/gsea")

	#grab only gene id and log2FC
	allgene_list <- allgene_list[,c(1,3)]
	colnames(allgene_list)[1] <- "gene"

	allgene_list <- allgene_list[order(-allgene_list$avg_log2FC),]
	gene_list <- allgene_list$avg_log2FC
	names(gene_list) <- allgene_list$gene
	print(head(gene_list))
	print(cells.list[i])

	fgsea <- fgseaMultilevel(GO_BP, gene_list, minSize = 5, maxSize = 500, gseaParam = 0)
	saveRDS(fgsea, paste0("Hippo_human_AD_",cells.list[i],"_fgsea.rds"))
	
	fgsea <- apply(fgsea,2,as.character)
	fgsea_df <- as.data.frame(fgsea)
	
	fgsea_df <- fgsea_df[order(fgsea_df$NES, decreasing = T),]
	write.csv(fgsea_df, paste0("Hippo_human_AD_",cells.list[i],"_fgsea.csv"))
}


