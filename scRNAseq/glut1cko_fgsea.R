#SLURM
salloc --mem=50G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate fgsea


#fGSEA-------------------------------------------------------------------------------------------------------------------
R

library(fgsea)
library(tidyverse)
library(ggplot2)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
annotations <- read.csv("annotations_Mus_EnsDb.csv")
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO")

GO_BP <- gmtPathways("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/m5.go.bp.v2023.2.Mm.symbols.gmt")

#"Microglia", "aNSCs", "qNSCs",
cells.list <- c("Ependymal",  "Oligodendrocytes")
cells.list <- c("Microglia", "Astrocytes")
cells.list <- c("Oligodendrocytes-1","Oligodendrocytes-2","Oligodendrocytes-3","Oligodendrocytes-4")
cells.list <- c("Glial progenitors-1", "Glial progenitors-2", "TAPs", "Neuroblasts-1", "Neuroblasts-2")
cells.list <- c("aNSCs", "qNSCs", "OPCs")

times.list <- c("6 months", "12 months")


for (i in 1:length(cells.list)) {
	for (j in 1:length(times.list)) {
		
		setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/DEG/cell/fem_gsea_in")
		filename <- paste0("Glutko_",cells.list[i], "_", times.list[j],"_markers_fem.gsea.csv")		
		
		#setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/DEG/cell/ODC v2")
		#filename <- paste0("Glutko_",cells.list[i],"_", times.list[j],"_markers_fem.gsea.csv")	
		
		allgene_list <- read.csv(file= filename, header=TRUE)
		setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")

		#grab only gene id and log2FC
		allgene_list <- allgene_list[,c(2,4)]

		allgene_list <- allgene_list[order(-allgene_list$avg_log2FC),]
		gene_list <- allgene_list$avg_log2FC
		names(gene_list) <- allgene_list$gene
		print(head(gene_list))
		print(cells.list[i])
		print(times.list[j])

		fgsea <- fgseaMultilevel(GO_BP, gene_list, minSize = 5, maxSize = 500, gseaParam = 0)
		saveRDS(fgsea, paste0(cells.list[i],"_",times.list[j],"_glut_fem_fgsea.rds"))
		
		fgsea <- apply(fgsea,2,as.character)
		fgsea_df <- as.data.frame(fgsea)
		
		fgsea_df <- fgsea_df[order(fgsea_df$NES, decreasing = T),]
		write.csv(fgsea_df, paste0(cells.list[i],"_",times.list[j],"_glut_fem_fgsea.csv"))
	}
}	