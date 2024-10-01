#Glut1KO Cellchat----------------------------------------------------------------------------------------------------
source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
source activate seurat4
salloc --mem=24G --cpus-per-task=12 --nodes=1 --ntasks=1 --time=08:00:00

#Seurat
#Split and analyze each condition separately then merge https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.Rmd
R
library(Seurat)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
load("GLUTKO.Robj", verbose = TRUE)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/cellchat")

Idents(GLUTKO) <- "gender"
females = subset(GLUTKO, idents = "Female")

#Remove peri and SM
Idents(females) <- "cell_ident"
females <- subset(x = females, idents = c("Smooth muscle cells", "Pericytes"), invert = TRUE)
females@meta.data$cell_ident <- droplevels(females@meta.data$cell_ident)
levels(females@meta.data$cell_ident)

#6 months
Idents(females) <- "age"
females_6mo = subset(females, idents = "6 months")

datasets_6mo <- SplitObject(females_6mo, split.by = "genotype")

#Confirm which dataset is which
head(datasets_6mo[[1]]@meta.data)
ctrl_seu <- datasets_6mo[[1]]
KO_seu  <- datasets_6mo[[2]]

# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
ctrl.input <- ctrl_seu[["RNA"]]$data
Idents(object = ctrl_seu) <- ctrl_seu$cell_ident
labels <- Idents(ctrl_seu)
ctrl_meta <- data.frame(labels = labels, row.names = names(labels))
save(ctrl.input,file="fem_6mo_ctrl_input.Robj")
save(ctrl_meta,file="fem_6mo_ctrl_meta.Robj")

KO_seu  <- datasets_6mo[[2]]
KO.input <- KO_seu[["RNA"]]$data
Idents(object = KO_seu) <- KO_seu$cell_ident
labels <- Idents(KO_seu)
KO_meta <- data.frame(labels = labels, row.names = names(labels))
save(KO.input,file="fem_6mo_KO_input.Robj")
save(KO_meta,file="fem_6mo_KO_meta.Robj")

#6 months analysis ------------------------------------------

#Cellchat 6 months
source activate cellchat
R
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/cellchat")

load(file="fem_6mo_KO_input.Robj")
load(file="fem_6mo_KO_meta.Robj")

load(file="fem_6mo_ctrl_input.Robj")
load(file="fem_6mo_ctrl_meta.Robj")

setwd("/work/biernaskie_lab/apun/nilesh_glut1/cellchat/6mo_fem")

cellchat.KO <- createCellChat(object = KO.input, meta = KO_meta, group.by = "labels")
cellchat.ctrl <- createCellChat(object = ctrl.input, meta = ctrl_meta, group.by = "labels")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
#Image required
showDatabaseCategory(CellChatDB)
#"Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use secreted CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

#KO Analysis
cellchat <- cellchat.KO

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#Threshold need pc/fc? eg thresh.pc = 0.1, thresh.fc = 0.1?
cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 50)
	
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = "KO_cellchat_all_df.net.csv")
df.net <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, file = "KO_cellchat_all_df.netP.csv")
	
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "cellchat_KO_6mo.rds")

#Ctrl Analysis
cellchat <- cellchat.ctrl

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#Threshold
cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.05)
cellchat <- identifyOverExpressedInteractions(cellchat)
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 50)
	
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = "ctrl_cellchat_all_df.net.csv")
df.net <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, file = "ctrl_cellchat_all_df.netP.csv")
	
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "cellchat_ctrl_6mo.rds")

#Merged analysis

cellchat.ctrl <- readRDS("cellchat_ctrl_6mo.rds")
cellchat.KO <- readRDS("cellchat_KO_6mo.rds")

#Care order (red = second dataset increased)
object.list <- list(Control = cellchat.ctrl, Glut1KO = cellchat.KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat, file = "cellchat_glut1ko_6mo_combined.rds")
cellchat <- readRDS("cellchat_glut1ko_6mo_combined.rds")

#NSC neuronal lineage
nsc_neu <- c("TAPs", "Neuroblasts-1", "Neuroblasts-2", "aNSCs", "qNSCs")
nsc_neu_epen <- c("Ependymal", "TAPs", "Neuroblasts-1", "Neuroblasts-2", "aNSCs", "qNSCs")

cellchat_nsc_epen <- subsetCellChat(cellchat, idents.use = nsc_neu_epen)

#Interactions
png(file = "Glut_fem_6mo_epen_nsc_neuronal_interaction_num.png", width = 7, height = 7, units = "cm", res = 900)
compareInteractions(cellchat_nsc_epen, show.legend = F, group = c(1,2))
dev.off()
png(file = "Glut_fem_6mo_epen_nsc_neuronal_interaction_weight.png", width = 7, height = 7, units = "cm", res = 900)
compareInteractions(cellchat_nsc_epen, show.legend = F, group = c(1,2), measure = "weight")
dev.off()

#NetVisual_diff
png(file = "Glut_fem_6mo_epen_nsc_neuronal_netDiff.png", width = 12, height = 12, units = "cm", res = 900)
par(xpd=TRUE)
netVisual_diffInteraction(cellchat_nsc_epen, weight.scale = T, label.edge = T)
dev.off()
png(file = "Glut_fem_6mo_epen_nsc_neuronal_netDiff_weight.png", width = 12, height = 12, units = "cm", res = 900)
par(xpd=TRUE)
netVisual_diffInteraction(cellchat_nsc_epen, weight.scale = T, measure = "weight", label.edge = T)
dev.off()

#Ranknet 
#no default reticulate env
#Compare the overall information flow of each signaling pathway or ligand-receptor pair
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#Can't do parallel
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
rankSimilarity(cellchat, type = "functional")

png(file = "Glut_fem_6mo_epen_nsc_neuronal_ranknet.png", width = 17, height = 8, units = "cm", res = 900)
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Ependymal", targets.use = nsc_neu, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Ependymal", targets.use = nsc_neu, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

png(file = "Glut_fem_6mo_nsc_neuronal_to_epen_ranknet.png", width = 17, height = 8, units = "cm", res = 900)
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = nsc_neu, targets.use = "Ependymal", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = nsc_neu, targets.use = "Ependymal", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
