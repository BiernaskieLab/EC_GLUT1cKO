#GSEA comparison all Glut, 5xfad, AD, apoe-------------------------------------------------------------------------------------------------------

#Overlap and venn---------------------------------------------------------------------------------

#Print out intersection
#https://www.biostars.org/p/9554414/
######################################
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
      }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
  }

overlapGroups <- function (listInput, sort = TRUE) {
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  for (i in 1:nrow(listInputunique)) {
	currentRow <- listInputunique[i,]
	myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
	attr(myelements, "groups") <- currentRow
	grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
	myelements
  }
  if (sort) {
	grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
}


cell_type <- c("Ependymal", "Astrocyte", "Oligodendrocyte", "Microglia")
for (i in 1:length(cell_type)) {

	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/gsea")
	AD <- read.csv(paste0("Hippo_human_AD_",cell_type[[i]],"_fgsea.csv"))

	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad/gsea")
	fad <- read.csv(paste0(cell_type[[i]],"_5xfad_fgsea.csv"))
	
	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/gsea")
	apoe <- read.csv(paste0(cell_type[[i]],"_apoe_fgsea.csv"))

	setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
	if(cell_type[[i]] %in% c("Astrocyte", "Oligodendrocyte")) {
		glut_6 <- read.csv(paste0(cell_type[[i]],"s_6 months_glut_fem_fgsea.csv"))
		glut_12 <- read.csv(paste0(cell_type[[i]],"s_12 months_glut_fem_fgsea.csv"))
	} else {
		glut_6 <- read.csv(paste0(cell_type[[i]],"_6 months_glut_fem_fgsea.csv"))
		glut_12 <- read.csv(paste0(cell_type[[i]],"_12 months_glut_fem_fgsea.csv"))
	}

	setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/venn")
	
	#Upregulated filter
	AD_filtered_up <- AD[AD$padj < 0.05 & AD$NES > 0, ]
	fad_filtered_up <- fad[fad$padj < 0.05 & fad$NES > 0, ]
	apoe_filtered_up <- apoe[apoe$padj < 0.05 & apoe$NES > 0, ]	
	glut_6_filtered_up <- glut_6[glut_6$padj < 0.05 & glut_6$NES > 0, ]	
	glut_12_filtered_up <- glut_12[glut_12$padj < 0.05 & glut_12$NES > 0, ]	
	
	
	#Downregulated filter
	AD_filtered_down <- AD[AD$padj < 0.05 & AD$NES < 0, ]
	fad_filtered_down <- fad[fad$padj < 0.05 & fad$NES < 0, ]
	apoe_filtered_down <- apoe[apoe$padj < 0.05 & apoe$NES < 0, ]
	glut_6_filtered_down <- glut_6[glut_6$padj < 0.05 & glut_6$NES < 0, ]	
	glut_12_filtered_down <- glut_12[glut_12$padj < 0.05 & glut_12$NES < 0, ]	

	#For writting to csv
	intersection_6_up <- list("AD" = AD_filtered_up$pathway, "5xFAD" = fad_filtered_up$pathway, "ApoE4" = apoe_filtered_up$pathway, "Glut1KO 6mo" = glut_6_filtered_up$pathway)
	intersection_12_up <- list("AD" = AD_filtered_up$pathway, "5xFAD" = fad_filtered_up$pathway, "ApoE4" = apoe_filtered_up$pathway,"Glut1KO 12mo" = glut_12_filtered_up$pathway)

	intersection_6_down <- list("AD" = AD_filtered_down$pathway, "5xFAD" = fad_filtered_down$pathway, "ApoE4" = apoe_filtered_down$pathway, "Glut1KO 6mo" = glut_6_filtered_down$pathway)
	intersection_12_down <- list("AD" = AD_filtered_down$pathway, "5xFAD" = fad_filtered_down$pathway, "ApoE4" = apoe_filtered_down$pathway, "Glut1KO 12mo" = glut_12_filtered_down$pathway)

	intersections <- list(intersection_6_up = intersection_6_up, intersection_12_up = intersection_12_up, intersection_6_down = intersection_6_down, intersection_12_down = intersection_12_down)

	for (u in 1:length(intersections)) {
		overlap_intersections <- overlapGroups(intersections[[u]])
		result_list <- purrr::map(overlap_intersections, ~ attr(overlap_intersections, "elements")[.x])
		name_list <- names(intersections[u])
		capture.output(result_list, file = paste0(cell_type[[i]],"_",name_list,".txt")) 
	}
	#Upreg Venn 6 mo
	png(file = paste0(cell_type[[i]],"_6_coupregulated_venn.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_6_up, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-upregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none",) + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()
	
	#12 mo
	png(file = paste0(cell_type[[i]],"_12_coupregulated_venn.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_12_up, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-upregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none") + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()
	
	
	#Downreg 6 mo
	png(file = paste0(cell_type[[i]],"_6_codownregulated_venn.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_6_down, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-downregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none") + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()
	
	#12 mo
	png(file = paste0(cell_type[[i]],"_12_codownregulated_venn.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_12_down, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-downregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none") + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()	
	
	png(file = paste0(cell_type[[i]],"_down_clean.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_12_down, show_percentage = FALSE, text_size = 0, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-downregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none") + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()

	png(file = paste0(cell_type[[i]],"_up_clean.jpeg"), width = 19, height = 13, units = "cm", res = 400)
	venn = ggvenn(intersection_12_up, show_percentage = FALSE, text_size = 0, fill_alpha = 0.15, set_name_size = 5) + labs(title = paste0(cell_type[[i]]," co-upregulated")) +
	theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)), legend.position = "none") + scale_x_continuous(expand = expansion(mult = .2))
	print(venn)
	dev.off()	
}	
######################################

#Overlap----------------------------------------
conda activate seurat4
library(ggvenn)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/")
source("overlap.R")
	


#DotPlot gsea terms------------------------------------------------------------------ 
library(stringr)
library(forcats)
library(cowplot)

#Ependymal------------------------------------------------------------------ ------------------------------------------------------------------ 
setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/gsea")
AD_epen <- read.csv(paste0("Hippo_human_AD_Ependymal_fgsea.csv"))
AD_epen$origin <- "Human AD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad/gsea")
fad_epen <- read.csv(paste0("Ependymal_5xfad_fgsea.csv"))
fad_epen$origin <- "5xFAD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/gsea")
apoe_epen <- read.csv(paste0("Ependymal_apoe_fgsea.csv"))
apoe_epen$origin <- "ApoE4"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_epen <- read.csv(paste0("Ependymal_6 months_glut_fem_fgsea.csv"))
glut_6_epen$origin <- "Glut1KO 6mo"

glut_12_epen <- read.csv(paste0("Ependymal_12 months_glut_fem_fgsea.csv"))
glut_12_epen$origin <- "Glut1KO 12mo"

merged_epen_12 <- rbind(AD_epen, fad_epen, apoe_epen, glut_12_epen)
merged_epen_6 <- rbind(AD_epen, fad_epen, apoe_epen, glut_6_epen)


setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/")

#Get count of leading edge genes	
merged_epen_12$count <- 0
for(i in 1:nrow(merged_epen_12)){
  edge <- strsplit(merged_epen_12$leadingEdge[i], ",")
  merged_epen_12$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_epen_12$pathway <- gsub("GOBP_", "", merged_epen_12$pathway)
#Only take sig terms
merged_epen_12 <- merged_epen_12[merged_epen_12$padj<0.05, ]

#Get count of leading edge genes	
merged_epen_6$count <- 0
for(i in 1:nrow(merged_epen_6)){
  edge <- strsplit(merged_epen_6$leadingEdge[i], ",")
  merged_epen_6$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_epen_6$pathway <- gsub("GOBP_", "", merged_epen_6$pathway)
#Only take sig terms
merged_epen_6 <- merged_epen_6[merged_epen_6$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

#Test decreasing 
decreasing_12_months_epen <- c("MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
"NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY",
"CILIUM_ORGANIZATION",
"AEROBIC_RESPIRATION",
"REGULATION_OF_AUTOPHAGY",
"MICROTUBULE_BASED_TRANSPORT",
"MOTILE_CILIUM_ASSEMBLY",
"DNA_INTEGRITY_CHECKPOINT_SIGNALING",
"REGULATION_OF_CYTOSKELETON_ORGANIZATION",
"RESPONSE_TO_OXIDATIVE_STRESS",
"CELL_MATRIX_ADHESION",
"REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
"REGULATION_OF_NEURON_APOPTOTIC_PROCESS",
"GLIOGENESIS",
"GLUCOSE_METABOLIC_PROCESS",
"POSITIVE_REGULATION_OF_NEUROGENESIS")

merged_epen_sub_down_12 <- subset(merged_epen_12, pathway %in% c(decreasing_12_months_epen))


#png(file = "Ependymal_12_gsea_glut_downreg_sig_venncol_test.png", width = 20, height = 16, units = "cm", res = 500)
p1 <- merged_epen_sub_down_12 %>%
  mutate(pathway=factor(pathway, levels=decreasing_12_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(-7, 5) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
#print(p1)
#dev.off()


decreasing_6_months_epen <- c("REGULATION_OF_CYTOSKELETON_ORGANIZATION",
"TRANSMEMBRANE_RECEPTOR_PROTEIN_SERINE_THREONINE_KINASE_SIGNALING_PATHWAY",
"REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_DISASSEMBLY",
"REGULATION_OF_NEUROGENESIS",
"EPITHELIAL_CELL_DEVELOPMENT",
"MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
"SECOND_MESSENGER_MEDIATED_SIGNALING",
"IMPORT_ACROSS_PLASMA_MEMBRANE",
"HIPPO_SIGNALING",
"CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS")  

merged_epen_sub_down_6 <- subset(merged_epen_6, pathway %in% c(decreasing_6_months_epen))


#png(file = "Ependymal_6_gsea_glut_downreg_sig_venncol.png", width = 26, height = 11, units = "cm", res = 500)
p2 <- merged_epen_sub_down_6 %>%
  mutate(pathway=factor(pathway, levels=decreasing_6_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("green", "red", "blue")) +
	ylim(-7, 5) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
#dev.off()


increasing_12_months_epen <- c("TYPE_II_INTERFERON_PRODUCTION",
"ADAPTIVE_IMMUNE_RESPONSE",
"B_CELL_PROLIFERATION",
"REGULATION_OF_LEUKOCYTE_PROLIFERATION",
"DEFENSE_RESPONSE_TO_GRAM_POSITIVE_BACTERIUM")

merged_epen_sub_up_12 <- subset(merged_epen_12, pathway %in% c(increasing_12_months_epen))


#png(file = "Ependymal_12_gsea_glut_upreg_sig_venncol.png", width = 19, height = 7, units = "cm", res = 500)
p3 <- merged_epen_sub_up_12 %>%
  mutate(pathway=factor(pathway, levels=increasing_12_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
#dev.off()


increasing_6_months_epen <- c("DEFENSE_RESPONSE_TO_GRAM_POSITIVE_BACTERIUM",
"RIBOSOMAL_SMALL_SUBUNIT_ASSEMBLY",
"TYPE_II_INTERFERON_PRODUCTION",
"POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE",
"ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY")

merged_epen_sub_up_6 <- subset(merged_epen_6, pathway %in% c(increasing_6_months_epen))

#png(file = "Ependymal_6_gsea_glut_upreg_sig_venncol.png", width = 26, height = 7, units = "cm", res = 500)
p4 <- merged_epen_sub_up_6 %>%
  mutate(pathway=factor(pathway, levels=increasing_6_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("green", "red", "blue")) + 
	ylim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
#dev.off()


#Standardize size

p_all <- align_plots(p1, p2, p3, p4, align="v")

ggsave2("Ependymal_12_gsea_glut_downreg_sig_venncol.png", p_all[[1]], dpi = 500, width = 7, height = 7)
ggsave2("Ependymal_6_gsea_glut_downreg_sig_venncol.png", p_all[[2]], dpi = 500, width = 7, height = 4.5)
ggsave2("Ependymal_12_gsea_glut_upreg_sig_venncol.png", p_all[[3]], dpi = 500, width = 7, height = 3)
ggsave2("Ependymal_6_gsea_glut_upreg_sig_venncol_new.png", p_all[[4]], dpi = 500, width = 7, height = 3)

#Legend bottom origin
png(file = "Legend bottom_origin.png", width = 7, height = 7, units = "in", res = 500)
p1 <- merged_epen_sub_down_12 %>%
  mutate(pathway=factor(pathway, levels=decreasing_12_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	scale_color_manual(values = c("yellow", "green", "red", "blue")) + 
	ylim(-7, 5) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
	theme(legend.position = "bottom", legend.box.just = "left") + 
	guides(
		color = guide_legend(order = 1, title = "Origin"),
		size = guide_legend(order = 2, title = "Count")) +
    coord_flip()
print(p1)
dev.off()
#Legend bottom count
png(file = "Legend bottom_count.png", width = 7, height = 7, units = "in", res = 500)
p1 <- merged_epen_sub_down_12 %>%
  mutate(pathway=factor(pathway, levels=decreasing_12_months_epen)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line 
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	ylim(-7, 5) +
	scale_color_manual(values = c("yellow", "green", "red", "blue")) + 
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
	theme(legend.position = "bottom", legend.box.just = "left") + 
	guides(
		color = guide_legend(order = 2, title = "Origin"),
		size = guide_legend(order = 1, title = "Count")) +
    coord_flip()
print(p1)
dev.off()


#Original padj colour matches
png(file = "test.png", width = 20, height = 18, units = "cm", res = 500)
ggplot(test, aes_string(x="NES", y="pathway", size="count", color="padj")) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE)) +
	ylab(NULL) + ggtitle("Test") + scale_size(range=c(3, 8))
dev.off()



#Astrocytes--------------------------------------------------------------------------------------

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/gsea")
AD_astro <- read.csv(paste0("Hippo_human_AD_Astrocyte_fgsea.csv"))
AD_astro$origin <- "Human AD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad/gsea")
fad_astro <- read.csv(paste0("Astrocyte_5xfad_fgsea.csv"))
fad_astro$origin <- "5xFAD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/gsea")
apoe_astro <- read.csv(paste0("Astrocyte_apoe_fgsea.csv"))
apoe_astro$origin <- "ApoE4"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_astro <- read.csv(paste0("Astrocytes_6 months_glut_fem_fgsea.csv"))
glut_6_astro$origin <- "Glut1KO 6mo"

glut_12_astro <- read.csv(paste0("Astrocytes_12 months_glut_fem_fgsea.csv"))
glut_12_astro$origin <- "Glut1KO 12mo"

merged_astro_12 <- rbind(AD_astro, fad_astro, apoe_astro, glut_12_astro)
merged_astro_6 <- rbind(AD_astro, fad_astro, apoe_astro, glut_6_astro)


setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/dot/astrocyte")

#Get count of leading edge genes	
merged_astro_12$count <- 0
for(i in 1:nrow(merged_astro_12)){
  edge <- strsplit(merged_astro_12$leadingEdge[i], ",")
  merged_astro_12$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_astro_12$pathway <- gsub("GOBP_", "", merged_astro_12$pathway)
#Only take sig terms
merged_astro_12 <- merged_astro_12[merged_astro_12$padj<0.05, ]

#Get count of leading edge genes	
merged_astro_6$count <- 0
for(i in 1:nrow(merged_astro_6)){
  edge <- strsplit(merged_astro_6$leadingEdge[i], ",")
  merged_astro_6$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_astro_6$pathway <- gsub("GOBP_", "", merged_astro_6$pathway)
#Only take sig terms
merged_astro_6 <- merged_astro_6[merged_astro_6$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

astro_6mo_up <- c("HUMORAL_IMMUNE_RESPONSE",
"NEUTROPHIL_CHEMOTAXIS",
"POSITIVE_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE",
"ACUTE_PHASE_RESPONSE",
"REGULATION_OF_AMIDE_METABOLIC_PROCESS",
"TYPE_II_INTERFERON_PRODUCTION",
"T_CELL_MEDIATED_IMMUNITY",
"GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION",
"ALPHA_BETA_T_CELL_ACTIVATION",
"LEUKOCYTE_MEDIATED_CYTOTOXICITY",
"LEUKOCYTE_MEDIATED_IMMUNITY",
"ADAPTIVE_IMMUNE_RESPONSE",
"LYMPHOCYTE_MEDIATED_IMMUNITY",
"INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS")

merged_astro_sub_up_6 <- subset(merged_astro_6, pathway %in% c(astro_6mo_up))

p1 <- merged_astro_sub_up_6 %>%
  mutate(pathway=factor(pathway, levels=astro_6mo_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

astro_6mo_down <- c("CYTOPLASMIC_TRANSLATION",
"RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
"CELLULAR_RESPIRATION",
"CELL_JUNCTION_ASSEMBLY",
"CELL_SUBSTRATE_ADHESION",
"REGULATION_OF_GTPASE_ACTIVITY",
"TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
"TRANSMEMBRANE_RECEPTOR_PROTEIN_SERINE_THREONINE_KINASE_SIGNALING_PATHWAY",
"RAS_PROTEIN_SIGNAL_TRANSDUCTION",
"CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS")

merged_astro_sub_down_6 <- subset(merged_astro_6, pathway %in% c(astro_6mo_down))

p2 <- merged_astro_sub_down_6 %>%
  mutate(pathway=factor(pathway, levels=astro_6mo_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

astro_12mo_up <- c("ACUTE_PHASE_RESPONSE",
"LEUKOCYTE_CELL_CELL_ADHESION",
"NEUTROPHIL_CHEMOTAXIS",
"HUMORAL_IMMUNE_RESPONSE",
"ADAPTIVE_IMMUNE_RESPONSE",
"LYMPHOCYTE_MEDIATED_IMMUNITY",
"LEUKOCYTE_MEDIATED_CYTOTOXICITY",
"MYELOID_LEUKOCYTE_ACTIVATION",
"TYPE_II_INTERFERON_PRODUCTION",
"T_CELL_MEDIATED_IMMUNITY")

merged_astro_sub_up_12 <- subset(merged_astro_12, pathway %in% c(astro_12mo_up))

p3 <- merged_astro_sub_up_12 %>%
  mutate(pathway=factor(pathway, levels=astro_12mo_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

astro_12mo_down <- c("TRANSLATION_AT_SYNAPSE",
"CYTOPLASMIC_TRANSLATION",
"RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
"NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
"CELLULAR_RESPIRATION",
"VESICLE_ORGANIZATION",
"GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",
"CELL_SUBSTRATE_ADHESION",
"CELL_JUNCTION_ASSEMBLY",
"TISSUE_MIGRATION",
"ACTOMYOSIN_STRUCTURE_ORGANIZATION",
"REGULATION_OF_GTPASE_ACTIVITY",
"TRANSMEMBRANE_RECEPTOR_PROTEIN_SERINE_THREONINE_KINASE_SIGNALING_PATHWAY",
"TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
"RAS_PROTEIN_SIGNAL_TRANSDUCTION",
"CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS",
"ORGANIC_ACID_TRANSMEMBRANE_TRANSPORT",
"NOTCH_SIGNALING_PATHWAY")

merged_astro_sub_down_12 <- subset(merged_astro_12, pathway %in% c(astro_12mo_down))

p4 <- merged_astro_sub_down_12 %>%
  mutate(pathway=factor(pathway, levels=astro_12mo_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	#Colorblind safe RBYG
	scale_color_manual(values = c("yellow", "green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

#Standardize size

p_all <- align_plots(p1, p2, p3, p4, align="v")
ggsave2("Astrocyte_6_gsea_glut_upreg_sig_venncol.png", p_all[[1]], dpi = 500, width = 7, height = 6.5)
ggsave2("Astrocyte_6_gsea_glut_downreg_sig_venncol.png", p_all[[2]], dpi = 500, width = 7, height = 4.25)
ggsave2("Astrocyte_12_gsea_glut_upreg_sig_venncol.png", p_all[[3]], dpi = 500, width = 7, height = 4.25)
ggsave2("Astrocyte_12_gsea_glut_downreg_sig_venncol_new.png", p_all[[4]], dpi = 500, width = 7, height = 7.25)



#Oligodendrocytes------------------------------------------------------------------ 
oligo_6mo_up <- c(
"MACROAUTOPHAGY",
"PROTEIN_POLYUBIQUITINATION",
"RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
"AUTOPHAGOSOME_ORGANIZATION",
"CELLULAR_RESPIRATION",
"AUTOPHAGY_OF_MITOCHONDRION",
"INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
"SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
"MEMBRANE_FISSION",
"MITOCHONDRIAL_FISSION")

oligo_6mo_down <- c(
"TISSUE_MIGRATION",
"CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
"FOCAL_ADHESION_ASSEMBLY",
"REGULATION_OF_NEUROGENESIS",
"OLIGODENDROCYTE_DIFFERENTIATION",
"CELL_FATE_COMMITMENT",
"ACTIN_FILAMENT_ORGANIZATION",
"RESPONSE_TO_CELL_CYCLE_CHECKPOINT_SIGNALING",
"MICROTUBULE_BASED_MOVEMENT",
"AXON_GUIDANCE")


oligo_12mo_up <- c(
"ADAPTIVE_IMMUNE_RESPONSE",
"NEUTROPHIL_MIGRATION",
"HUMORAL_IMMUNE_RESPONSE",
"LEUKOCYTE_MEDIATED_IMMUNITY",
"LYMPHOCYTE_MEDIATED_IMMUNITY",
"COMPLEMENT_ACTIVATION_ALTERNATIVE_PATHWAY",
"INTERLEUKIN_1_ALPHA_PRODUCTION",
"TYPE_II_INTERFERON_PRODUCTION",
"MESENCHYMAL_CELL_APOPTOTIC_PROCESS",
"PYROPTOSIS")

oligo_12mo_down <- c(
"TISSUE_MIGRATION",
"CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
"FOCAL_ADHESION_ASSEMBLY",
"PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWA",
"INTERMEDIATE_FILAMENT_BASED_PROCESS",
"REGULATION_OF_PLATELET_ACTIVATION",
"AEROBIC_RESPIRATION",
"REGULATION_OF_NEUROGENESIS",
"LYSOSOMAL_TRANSPORT",
"GLIOGENESIS")


setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/human/gsea")
AD_oligo <- read.csv(paste0("Hippo_human_AD_Oligodendrocyte_fgsea.csv"))
AD_oligo$origin <- "Human AD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/5xfad/gsea")
fad_oligo <- read.csv(paste0("Oligodendrocyte_5xfad_fgsea.csv"))
fad_oligo$origin <- "5xFAD"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/apoe/gsea")
apoe_oligo <- read.csv(paste0("Oligodendrocyte_apoe_fgsea.csv"))
apoe_oligo$origin <- "ApoE4"

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_oligo <- read.csv(paste0("Oligodendrocytes_6 months_glut_fem_fgsea.csv"))
glut_6_oligo$origin <- "Glut1KO 6mo"

glut_12_oligo <- read.csv(paste0("Oligodendrocytes_12 months_glut_fem_fgsea.csv"))
glut_12_oligo$origin <- "Glut1KO 12mo"

merged_oligo_12 <- rbind(AD_oligo, fad_oligo, apoe_oligo, glut_12_oligo)
merged_oligo_6 <- rbind(AD_oligo, fad_oligo, apoe_oligo, glut_6_oligo)


setwd("/work/biernaskie_lab/apun/nilesh_glut1/human_AD/dot/oligodendrocyte")

#Get count of leading edge genes	
merged_oligo_12$count <- 0
for(i in 1:nrow(merged_oligo_12)){
  edge <- strsplit(merged_oligo_12$leadingEdge[i], ",")
  merged_oligo_12$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_oligo_12$pathway <- gsub("GOBP_", "", merged_oligo_12$pathway)
#Only take sig terms
merged_oligo_12 <- merged_oligo_12[merged_oligo_12$padj<0.05, ]

#Get count of leading edge genes	
merged_oligo_6$count <- 0
for(i in 1:nrow(merged_oligo_6)){
  edge <- strsplit(merged_oligo_6$leadingEdge[i], ",")
  merged_oligo_6$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_oligo_6$pathway <- gsub("GOBP_", "", merged_oligo_6$pathway)
#Only take sig terms
merged_oligo_6 <- merged_oligo_6[merged_oligo_6$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

merged_oligo_sub_up_6 <- subset(merged_oligo_6, pathway %in% c(oligo_6mo_up))

p1 <- merged_oligo_sub_up_6 %>%
  mutate(pathway=factor(pathway, levels=oligo_6mo_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	scale_color_manual(values = c("green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

merged_oligo_sub_up_12 <- subset(merged_oligo_12, pathway %in% c(oligo_12mo_up))
	
p2 <- merged_oligo_sub_up_12 %>%
  mutate(pathway=factor(pathway, levels=oligo_12mo_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-upregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	scale_color_manual(values = c("green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
	

merged_oligo_sub_down_6 <- subset(merged_oligo_6, pathway %in% c(oligo_6mo_down))
	
p3 <- merged_oligo_sub_down_6 %>%
  mutate(pathway=factor(pathway, levels=oligo_6mo_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	scale_color_manual(values = c("green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()
	

merged_oligo_sub_down_12 <- subset(merged_oligo_12, pathway %in% c(oligo_12mo_down))
	
p4 <- merged_oligo_sub_down_12 %>%
  mutate(pathway=factor(pathway, levels=oligo_12mo_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=origin)) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO Co-downregulated GO Biological Process Terms") +
    theme(plot.title = element_text(hjust=1), axis.text.x = element_text(size = 14)) +
	scale_color_manual(values = c("green", "red", "blue")) +
	ylim(-10, 10) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

#Standardize size
p_all <- align_plots(p1, p2, p3, p4, align="v")
ggsave2("Oligodendrocyte_6_gsea_glut_upreg_sig_venncol.png", p_all[[1]], dpi = 500, width = 7, height = 4.25)
ggsave2("Oligodendrocyte_12_gsea_glut_upreg_sig_venncol.png", p_all[[2]], dpi = 500, width = 7, height = 4.25)
ggsave2("Oligodendrocyte_6_gsea_glut_downreg_sig_venncol.png", p_all[[3]], dpi = 500, width = 7, height = 4.25)
ggsave2("Oligodendrocyte_12_gsea_glut_downreg_sig_venncol_new.png", p_all[[4]], dpi = 500, width = 7, height = 4.25)





#Astrocytes just glut------------------------------------------------------------------ 

astro_6_up <- c("HUMORAL_IMMUNE_RESPONSE",
"NEUTROPHIL_CHEMOTAXIS",
"POSITIVE_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE",
"ACUTE_PHASE_RESPONSE",
"TYPE_II_INTERFERON_PRODUCTION",
"T_CELL_MEDIATED_IMMUNITY",
"GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION",
"ALPHA_BETA_T_CELL_ACTIVATION",
"LEUKOCYTE_MEDIATED_CYTOTOXICITY",
"LEUKOCYTE_MEDIATED_IMMUNITY",
"ADAPTIVE_IMMUNE_RESPONSE",
"LYMPHOCYTE_MEDIATED_IMMUNITY",
"INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS")

astro_6_down <- c("CYTOPLASMIC_TRANSLATION",
"RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
"CELLULAR_RESPIRATION",
"CELL_JUNCTION_ASSEMBLY",
"CELL_SUBSTRATE_ADHESION",
"REGULATION_OF_GTPASE_ACTIVITY",
"TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
"TRANSMEMBRANE_RECEPTOR_PROTEIN_SERINE_THREONINE_KINASE_SIGNALING_PATHWAY",
"RAS_PROTEIN_SIGNAL_TRANSDUCTION",
"CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS")

astro_12_up <- c("ACUTE_PHASE_RESPONSE",
"LEUKOCYTE_CELL_CELL_ADHESION",
"NEUTROPHIL_CHEMOTAXIS",
"HUMORAL_IMMUNE_RESPONSE",
"ADAPTIVE_IMMUNE_RESPONSE",
"LYMPHOCYTE_MEDIATED_IMMUNITY",
"LEUKOCYTE_MEDIATED_CYTOTOXICITY",
"MYELOID_LEUKOCYTE_ACTIVATION",
"TYPE_II_INTERFERON_PRODUCTION",
"T_CELL_MEDIATED_IMMUNITY")

astro_12_down <- c("TRANSLATION_AT_SYNAPSE",
"CYTOPLASMIC_TRANSLATION",
"RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
"NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
"CELLULAR_RESPIRATION",
"VESICLE_ORGANIZATION",
"GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",
"CELL_SUBSTRATE_ADHESION",
"CELL_JUNCTION_ASSEMBLY",
"TISSUE_MIGRATION",
"ACTOMYOSIN_STRUCTURE_ORGANIZATION",
"REGULATION_OF_GTPASE_ACTIVITY",
"TRANSMEMBRANE_RECEPTOR_PROTEIN_SERINE_THREONINE_KINASE_SIGNALING_PATHWAY",
"TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
"RAS_PROTEIN_SIGNAL_TRANSDUCTION",
"CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS",
"ORGANIC_ACID_TRANSMEMBRANE_TRANSPORT",
"NOTCH_SIGNALING_PATHWAY")

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_astro <- read.csv(paste0("Astrocytes_6 months_glut_fem_fgsea.csv"))
glut_6_astro$origin <- "Glut1KO 6mo"

glut_12_astro <- read.csv(paste0("Astrocytes_12 months_glut_fem_fgsea.csv"))
glut_12_astro$origin <- "Glut1KO 12mo"

#Get count of leading edge genes	
glut_12_astro$count <- 0
for(i in 1:nrow(glut_12_astro)){
  edge <- strsplit(glut_12_astro$leadingEdge[i], ",")
  glut_12_astro$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_12_astro$pathway <- gsub("GOBP_", "", glut_12_astro$pathway)
#Only take sig terms
glut_12_astro <- glut_12_astro[glut_12_astro$padj<0.05, ]

#Get count of leading edge genes	
glut_6_astro$count <- 0
for(i in 1:nrow(glut_6_astro)){
  edge <- strsplit(glut_6_astro$leadingEdge[i], ",")
  glut_6_astro$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_6_astro$pathway <- gsub("GOBP_", "", glut_6_astro$pathway)
#Only take sig terms
glut_6_astro <- glut_6_astro[glut_6_astro$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

glut_6_astro_up_sub <- subset(glut_6_astro, pathway %in% c(astro_6_up))
levels_astro = glut_6_astro_up_sub$pathway[order(glut_6_astro_up_sub$NES)]
p1 <- glut_6_astro_up_sub %>%
  mutate(pathway=factor(pathway, levels=levels_astro)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Astrocyte Upregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

glut_6_astro_down_sub <- subset(glut_6_astro, pathway %in% c(astro_6_down))
levels_astro = glut_6_astro_down_sub$pathway[order(glut_6_astro_down_sub$NES)]
p2 <- glut_6_astro_down_sub %>%
  mutate(pathway=factor(pathway, levels=levels_astro)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Astrocyte Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")
	
glut_12_astro_up_sub <- subset(glut_12_astro, pathway %in% c(astro_12_up))
levels_astro = glut_12_astro_up_sub$pathway[order(glut_12_astro_up_sub$NES)]
p3 <- glut_12_astro_up_sub %>%
  mutate(pathway=factor(pathway, levels=levels_astro)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Astrocyte Upregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")
	
glut_12_astro_down_sub <- subset(glut_12_astro, pathway %in% c(astro_12_down))
levels_astro = glut_12_astro_down_sub$pathway[order(glut_12_astro_down_sub$NES)]
p4 <- glut_12_astro_down_sub %>%
  mutate(pathway=factor(pathway, levels=levels_astro)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Astrocyte Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

p_all <- align_plots(p1, p2, p3, p4, align="v")
ggsave2("Astrocyte_6_gsea_glut_upreg_sig_venncol.png", p_all[[1]], dpi = 500, width = 7, height = 6)
ggsave2("Astrocyte_6_gsea_glut_downreg_sig_venncol.png", p_all[[2]], dpi = 500, width = 7, height = 4.25)
ggsave2("Astrocyte_12_gsea_glut_upreg_sig_venncol.png", p_all[[3]], dpi = 500, width = 7, height = 4.25)
ggsave2("Astrocyte_12_gsea_glut_downreg_sig_venncol.png", p_all[[4]], dpi = 500, width = 7, height = 7.25)

#Microglia just glut--------------------------------------------------------------------------------------

micro_6_up <- c("ANTIMICROBIAL_HUMORAL_RESPONSE",
"NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION",
"HUMORAL_IMMUNE_RESPONSE",
"ORGAN_OR_TISSUE_SPECIFIC_IMMUNE_RESPONSE",
"PHOSPHATIDYLSERINE_EXPOSURE_ON_APOPTOTIC_CELL_SURFACE",
"B_CELL_ADHESION",
"INNATE_IMMUNE_RESPONSE_IN_MUCOSA",
"RESPONSE_TO_PROSTAGLANDIN_D",
"DISACCHARIDE_METABOLIC_PROCESS",
"POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
"POSITIVE_REGULATION_OF_CORTICOSTEROID_HORMONE_SECRETION",
"KILLING_OF_CELLS_OF_ANOTHER_ORGANISM")

micro_6_down <- c(
"CELLULAR_RESPIRATION",
"CELLULAR_RESPONSE_TO_LIPID",
"LYSOSOMAL_TRANSPORT",
"REGULATION_OF_NEUROGENESIS",
"GLIOGENESIS",
"GLYCEROPHOSPHOLIPID_METABOLIC_PROCESS",
"CELLULAR_RESPONSE_TO_ORGANIC_CYCLIC_COMPOUND",
"LIPOPROTEIN_METABOLIC_PROCESS",
"CELLULAR_LIPID_CATABOLIC_PROCESS",
"FATTY_ACID_CATABOLIC_PROCESS",
"AMYLOID_BETA_METABOLIC_PROCESS",
"SPHINGOMYELIN_METABOLIC_PROCESS",
"LIPID_STORAGE")

micro_12_down <- c("CELLULAR_RESPIRATION",
"GLIOGENESIS",
"MICROTUBULE_BASED_TRANSPORT",
"LYSOSOMAL_TRANSPORT",
"ACTIN_FILAMENT_BASED_MOVEMENT",
"ACTIN_POLYMERIZATION_OR_DEPOLYMERIZATION",
"INTRACELLULAR_LIPID_TRANSPORT",
"PHOSPHOLIPID_METABOLIC_PROCESS",
"INTRACELLULAR_STEROL_TRANSPORT",
"PHOSPHOLIPID_BIOSYNTHETIC_PROCESS",
"LATE_ENDOSOME_TO_LYSOSOME_TRANSPORT",
"REGULATION_OF_LYSOSOME_SIZE",
"SPHINGOMYELIN_METABOLIC_PROCESS",
"FATTY_ACID_BETA_OXIDATION")

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_micro <- read.csv(paste0("Microglia_6 months_glut_fem_fgsea.csv"))
glut_6_micro$origin <- "Glut1KO 6mo"

glut_12_micro <- read.csv(paste0("Microglia_12 months_glut_fem_fgsea.csv"))
glut_12_micro$origin <- "Glut1KO 12mo"

#Get count of leading edge genes	
glut_12_micro$count <- 0
for(i in 1:nrow(glut_12_micro)){
  edge <- strsplit(glut_12_micro$leadingEdge[i], ",")
  glut_12_micro$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_12_micro$pathway <- gsub("GOBP_", "", glut_12_micro$pathway)
#Only take sig terms
glut_12_micro <- glut_12_micro[glut_12_micro$padj<0.05, ]

#Get count of leading edge genes	
glut_6_micro$count <- 0
for(i in 1:nrow(glut_6_micro)){
  edge <- strsplit(glut_6_micro$leadingEdge[i], ",")
  glut_6_micro$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_6_micro$pathway <- gsub("GOBP_", "", glut_6_micro$pathway)
#Only take sig terms
glut_6_micro <- glut_6_micro[glut_6_micro$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

glut_6_micro_up_sub <- subset(glut_6_micro, pathway %in% c(micro_6_up))

p1 <- glut_6_micro_up_sub %>%
  mutate(pathway=factor(pathway, levels=micro_6_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Microglia Upregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(0, 3) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

glut_6_micro_down_sub <- subset(glut_6_micro, pathway %in% c(micro_6_down))

p2 <- glut_6_micro_down_sub %>%
  mutate(pathway=factor(pathway, levels=micro_6_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Microglia Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")
	
glut_12_micro_down_sub <- subset(glut_12_micro, pathway %in% c(micro_12_down))

p3 <- glut_12_micro_down_sub %>%
  mutate(pathway=factor(pathway, levels=micro_12_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Microglia Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")
	
p_all <- align_plots(p1, p2, p3, align="v")

ggsave2("Microglia_6_gsea_glut_upreg_sig.png", p_all[[1]], dpi = 500, width = 7, height = 4.75)
ggsave2("Microglia_6_gsea_glut_downreg_sig.png", p_all[[2]], dpi = 500, width = 7, height = 4.75)
ggsave2("Microglia_12_gsea_glut_downreg_sig.png", p_all[[3]], dpi = 500, width = 7, height = 5)


#Oligodendrocytes clus just glut
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	

glut_6_oligo1 <- read.csv(paste0("Oligodendrocytes-1_6 months_glut_fem_fgsea.csv"))
glut_6_oligo1$cluster <- "Glut1KO 6mo Oligo-1"

glut_6_oligo2 <- read.csv(paste0("Oligodendrocytes-2_6 months_glut_fem_fgsea.csv"))
glut_6_oligo2$cluster <- "Glut1KO 6mo Oligo-2"

glut_6_oligo3 <- read.csv(paste0("Oligodendrocytes-3_6 months_glut_fem_fgsea.csv"))
glut_6_oligo3$cluster <- "Glut1KO 6mo Oligo-3"

glut_6_oligo4 <- read.csv(paste0("Oligodendrocytes-4_6 months_glut_fem_fgsea.csv"))
glut_6_oligo4$cluster <- "Glut1KO 6mo Oligo-4"

glut_12_oligo1 <- read.csv(paste0("Oligodendrocytes-1_12 months_glut_fem_fgsea.csv"))
glut_12_oligo1$cluster <- "Glut1KO 12mo Oligo-1"

glut_12_oligo2 <- read.csv(paste0("Oligodendrocytes-2_12 months_glut_fem_fgsea.csv"))
glut_12_oligo2$cluster <- "Glut1KO 12mo Oligo-2"

glut_12_oligo3 <- read.csv(paste0("Oligodendrocytes-3_12 months_glut_fem_fgsea.csv"))
glut_12_oligo3$cluster <- "Glut1KO 12mo Oligo-3"

glut_12_oligo4 <- read.csv(paste0("Oligodendrocytes-4_12 months_glut_fem_fgsea.csv"))
glut_12_oligo4$cluster <- "Glut1KO 12mo Oligo-4"

merged_oligo_12 <- rbind(glut_12_oligo1, glut_12_oligo2, glut_12_oligo3, glut_12_oligo4)
merged_oligo_6 <- rbind(glut_6_oligo1, glut_6_oligo2, glut_6_oligo3, glut_6_oligo4)


#Get count of leading edge genes	
merged_oligo_12$count <- 0
for(i in 1:nrow(merged_oligo_12)){
  edge <- strsplit(merged_oligo_12$leadingEdge[i], ",")
  merged_oligo_12$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_oligo_12$pathway <- gsub("GOBP_", "", merged_oligo_12$pathway)
#Only take sig terms
merged_oligo_12 <- merged_oligo_12[merged_oligo_12$padj<0.05, ]

#Get count of leading edge genes	
merged_oligo_6$count <- 0
for(i in 1:nrow(merged_oligo_6)){
  edge <- strsplit(merged_oligo_6$leadingEdge[i], ",")
  merged_oligo_6$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
merged_oligo_6$pathway <- gsub("GOBP_", "", merged_oligo_6$pathway)
#Only take sig terms
merged_oligo_6 <- merged_oligo_6[merged_oligo_6$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}
  
lipid_meta <- c(
"STEROL_METABOLIC_PROCESS",
"SPHINGOLIPID_METABOLIC_PROCESS",
"GLYCEROPHOSPHOLIPID_METABOLIC_PROCESS",
"GLYCEROLIPID_METABOLIC_PROCESS",
"FATTY_ACID_METABOLIC_PROCESS")

merged_oligo_6 <- subset(merged_oligo_6, pathway %in% c(lipid_meta))

p1 <- merged_oligo_6 %>%
  mutate(pathway=factor(pathway, levels=lipid_meta)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=cluster)) +
	scale_x_discrete(drop = FALSE) +
    scale_colour_discrete(drop = FALSE, limits = c("Glut1KO 6mo Oligo-1", "Glut1KO 6mo Oligo-2", "Glut1KO 6mo Oligo-3", "Glut1KO 6mo Oligo-4")) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Oligodendorcyte Downregulated Lipid GO Biological Process Terms") + scale_size(range=c(3, 8)) +
    theme(plot.title = element_text(hjust=1, size=8), axis.text.x = element_text(size = 14)) +
	ylim(-5, 5) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

merged_oligo_12 <- subset(merged_oligo_12, pathway %in% c(lipid_meta))

p2 <- merged_oligo_12 %>%
  mutate(pathway=factor(pathway, levels=lipid_meta)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=pathway, y=NES, size=count, color=cluster)) +
	scale_x_discrete(drop = FALSE) +
	    scale_colour_discrete(drop = FALSE) +
	geom_point(alpha = 0.3) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Oligodendorcyte Downregulated Lipid GO Biological Process Terms") + scale_size(range=c(3, 8)) +
    theme(plot.title = element_text(hjust=1, size=8), axis.text.x = element_text(size = 14)) +
	ylim(-5, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_hline(yintercept=c(0), linetype="dashed", color = "red") +
    coord_flip()

p_all <- align_plots(p1, p2, align="v")

ggsave2("Oligodendrocyte_6_gsea_glut_lipid_sig.png", p_all[[1]], dpi = 500, width = 7, height = 3.5)
ggsave2("Oligodendrocyte_12_gsea_glut_lipid_sig.png", p_all[[2]], dpi = 500, width = 7, height = 3.5)

#leading edge dot plot
fatty_acid_meta <- c("Ppargc1a", "Tnfrsf1a", "Alox5ap", "Ltc4s", "Gstp1")
glycerolip_glycerophos_meta <- c("Inpp4b", "Plcg2", "Pla2g7", "Plscr1", "Pla2g3")
sphingolipid_meta <- c("Cerkl", "Hacd4", "B4galt4", "Acer2", "St6galnac6")
sterol_meta <- c("Ephx2", "Cyp46a1", "Soat2", "Lrp5", "Tm7sf2")
#"Ephx2" top sterol but duplicated in FA so replaced
lipid_meta <- c(fatty_acid_meta, glycerolip_glycerophos_meta, sphingolipid_meta, sterol_meta)

#Oligodendrocyte glut------------------------------------------------------------------ 
library(readr)
setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/")

#Subset
load("GLUTKO.Robj")
Idents(GLUTKO) <- "gender"
females <- subset(GLUTKO, idents = "Female")
Idents(females) <- "cell_ident"
females_oligo <- subset(females, idents = c("Oligodendrocytes"))
#Idents(females_oligo) <- "age"
#females_oligo_6mo <- subset(females_oligo, idents = c("6 months"))
#females_oligo_12mo <- subset(females_oligo, idents = c("12 months"))

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")
#Combined 
png(file = "fem_oligo_lipid_meta_6_12_fil_top_down_dot_pval.png", width = 13, height = 22, units = "cm", res = 700)
DotPlot(females_oligo, features = lipid_meta, group.by = "sample", scale=FALSE, dot.min = 0, dot.scale = 8) + coord_flip() + RotatedAxis() + theme(axis.text.x = element_text(size = 8))
dev.off()

png(file = "fem_oligo_lipid_meta_6_12_fil_top_down_dot_pval_scaled.png", width = 13, height = 22, units = "cm", res = 700)
DotPlot(females_oligo, features = lipid_meta, group.by = "sample", scale=TRUE, dot.min = 0, dot.scale = 8) + coord_flip() + RotatedAxis() + theme(axis.text.x = element_text(size = 8))
dev.off()



#12 mo downreg glut overall
oligo_12_down <- c("CELLULAR_RESPIRATION",
"PHOSPHOLIPID_METABOLIC_PROCESS",
"GLYCEROLIPID_METABOLIC_PROCESS",
"GLYCEROPHOSPHOLIPID_METABOLIC_PROCESS",
"POST_GOLGI_VESICLE_MEDIATED_TRANSPORT",
"GLYCEROPHOSPHOLIPID_BIOSYNTHETIC_PROCESS",
"GOLGI_TO_PLASMA_MEMBRANE_TRANSPORT",
"LIPOPROTEIN_METABOLIC_PROCESS",
"FATTY_ACID_METABOLIC_PROCESS",
"MEMBRANE_LIPID_METABOLIC_PROCESS",
"PHOSPHATIDYLINOSITOL_BIOSYNTHETIC_PROCESS",
"REGULATION_OF_LIPID_METABOLIC_PROCESS",
"PHOSPHATIDYLCHOLINE_ACYL_CHAIN_REMODELING",
"STEROL_METABOLIC_PROCESS",
"GLYCOLIPID_BIOSYNTHETIC_PROCESS",
"REGULATION_OF_FATTY_ACID_BIOSYNTHETIC_PROCESS",
"LIPID_DROPLET_ORGANIZATION",
"REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS",
"SPHINGOLIPID_METABOLIC_PROCESS",
"NEUTRAL_LIPID_METABOLIC_PROCESS")

glut_12_oligo <- read.csv(paste0("Oligodendrocytes_12 months_glut_fem_fgsea.csv"))

#Get count of leading edge genes	
glut_12_oligo$count <- 0
for(i in 1:nrow(glut_12_oligo)){
  edge <- strsplit(glut_12_oligo$leadingEdge[i], ",")
  glut_12_oligo$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_12_oligo$pathway <- gsub("GOBP_", "", glut_12_oligo$pathway)
#Only take sig terms
glut_12_oligo <- glut_12_oligo[glut_12_oligo$padj<0.05, ]

str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

glut_12_oligo_down_sub <- subset(glut_12_oligo, pathway %in% c(oligo_12_down))
levels_oligo = glut_12_oligo_down_sub$pathway[order(glut_12_oligo_down_sub$NES)]
p0 <- glut_12_oligo_down_sub %>%
  mutate(pathway=factor(pathway, levels=levels_oligo)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.05)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Oligodendrocyte Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

ggsave2("Oligodendrocyte_12_gsea_glut_downreg_sig.png", p0, dpi = 500, width = 7, height = 8)

#Neuroblasts just glut1ko 6--------------------------

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_neurob <- read.csv(paste0("Neuroblasts-1_6 months_glut_fem_fgsea.csv"))
glut_6_neurob$origin <- "Glut1KO 6mo"

#Get count of leading edge genes	
glut_6_neurob$count <- 0
for(i in 1:nrow(glut_6_neurob)){
  edge <- strsplit(glut_6_neurob$leadingEdge[i], ",")
  glut_6_neurob$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_6_neurob$pathway <- gsub("GOBP_", "", glut_6_neurob$pathway)
#Only take sig terms
glut_6_neurob <- glut_6_neurob[glut_6_neurob$padj<0.05, ]


down_6_months_neurob <- c(
"REGULATION_OF_MITOTIC_CELL_CYCLE",
"DENDRITE_DEVELOPMENT",
"REGULATION_OF_CYTOSKELETON_ORGANIZATION",
"REGULATION_OF_CELL_GROWTH",
"REGULATION_OF_NEUROGENESIS",
"REGULATION_OF_AXONOGENESIS",
"SPINDLE_ORGANIZATION",
"CELL_CYCLE_G1_S_PHASE_TRANSITION",
"REGULATION_OF_SYNAPTIC_PLASTICITY",
"NEURON_MIGRATION")

glut_6_neurobsub_down_6 <- subset(glut_6_neurob, pathway %in% c(down_6_months_neurob))

png(file = "Neuroblast1_6_gsea_glut_downreg_sig_count.png", width = 5.5, height = 5, units = "in", res = 500)
glut_6_neurobsub_down_6 %>%
  mutate(pathway=factor(pathway, levels=down_6_months_neurob)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " "))) %>%
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Neuroblast-1 Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red") +
    theme(legend.position = "bottom", legend.box.just = "left")
dev.off()

png(file = "Neuroblast1_6_gsea_glut_downreg_sig_padj.png", width = 5.5, height = 5, units = "in", res = 500)
glut_6_neurobsub_down_6 %>%
  mutate(pathway=factor(pathway, levels=down_6_months_neurob)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " "))) %>%
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Neuroblast-1 Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-6, 0) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red") +
    theme(legend.position = "bottom", legend.box.just = "left") +
	guides(colour = guide_colourbar(order = 1, barwidth = 10), size = guide_legend(order = 2))
dev.off()

#Ependymal 6 and 12 months glut---------------------------------------
library(dplyr)
library(stringr)
library(forcats)


epen_6months_up <- c("MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
    "NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY",
    "VERY_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
    "PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
    "FATTY_ACID_BIOSYNTHETIC_PROCESS",
    "POSITIVE_REGULATION_OF_LEUKOCYTE_MEDIATED_IMMUNITY",
    "ACTIVATION_OF_INNATE_IMMUNE_RESPONSE",
    "FATTY_ACID_METABOLIC_PROCESS")

epen_6months_down <- c("REGULATION_OF_NEUROGENESIS",
	"REGULATION_OF_CYTOSKELETON_ORGANIZATION",
	"CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",
	"IMPORT_ACROSS_PLASMA_MEMBRANE",
	"MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
	"CILIUM_ORGANIZATION")	

epen_12months_up <- c("HUMORAL_IMMUNE_RESPONSE",
    "ADAPTIVE_IMMUNE_RESPONSE",
    "LEUKOCYTE_MEDIATED_IMMUNITY",
    "LYMPHOCYTE_MEDIATED_IMMUNITY",
    "REGULATION_OF_IMMUNE_EFFECTOR_PROCESS",
    "LEUKOCYTE_MEDIATED_CYTOTOXICITY",
    "TYPE_II_INTERFERON_PRODUCTION",
    "REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
    "B_CELL_ADHESION",
    "T_HELPER_1_TYPE_IMMUNE_RESPONSE",
    "B_CELL_PROLIFERATION",
    "T_HELPER_1_CELL_CYTOKINE_PRODUCTION",
    "REGULATION_OF_MYELOID_DENDRITIC_CELL_ACTIVATION")

epen_12months_down <- c("CILIUM_ORGANIZATION",
	"ACTIN_FILAMENT_ORGANIZATION",
    "REGULATION_OF_AUTOPHAGY",
    "REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
    "PHOSPHOLIPID_METABOLIC_PROCESS",
    "REGULATION_OF_NEUROGENESIS",
    "CELL_CELL_SIGNALING_BY_WNT",
    "MOTILE_CILIUM_ASSEMBLY",
    "GLIOGENESIS",
    "ESTABLISHMENT_OF_CELL_POLARITY",
    "MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
	"CILIUM_MOVEMENT",
    "CELLULAR_SENESCENCE",
    "CELL_MATRIX_ADHESION",
	"CEREBROSPINAL_FLUID_CIRCULATION",
	"CELL_CELL_JUNCTION_ORGANIZATION",
	"REGULATION_OF_FATTY_ACID_OXIDATION",
	"GLUCOSE_IMPORT",
	"OLIGODENDROCYTE_DIFFERENTIATION",
	"REGULATION_OF_CELL_JUNCTION_ASSEMBLY")
	
	



setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat/GO/fgsea_alex_new")	
glut_6_epen <- read.csv(paste0("Ependymal_6 months_glut_fem_fgsea.csv"))
glut_6_epen$origin <- "Glut1KO 6mo"

glut_12_epen <- read.csv(paste0("Ependymal_12 months_glut_fem_fgsea.csv"))
glut_12_epen$origin <- "Glut1KO 12mo"

#Get count of leading edge genes	
glut_12_epen$count <- 0
for(i in 1:nrow(glut_12_epen)){
  edge <- strsplit(glut_12_epen$leadingEdge[i], ",")
  glut_12_epen$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_12_epen$pathway <- gsub("GOBP_", "", glut_12_epen$pathway)
#Only take sig terms
#glut_12_epen <- glut_12_epen[glut_12_epen$padj<0.05, ]

#Get count of leading edge genes	
glut_6_epen$count <- 0
for(i in 1:nrow(glut_6_epen)){
  edge <- strsplit(glut_6_epen$leadingEdge[i], ",")
  glut_6_epen$count[i] = length(edge[[1]])
}
#Strip GOBP_ from all terms in pathway
glut_6_epen$pathway <- gsub("GOBP_", "", glut_6_epen$pathway)
#Only take sig terms
#glut_6_epen <- glut_6_epen[glut_6_epen$padj<0.05, ]


str_wrap_factor <- function(x, ...) {
  levels(x) <- str_wrap(levels(x), ...)
  x
}

glut_6_epen_up_sub <- subset(glut_6_epen, pathway %in% c(epen_6months_up))

p1 <- glut_6_epen_up_sub %>%
  mutate(pathway=factor(pathway, levels=epen_6months_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Ependymal Upregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

glut_6_epen_down_sub <- subset(glut_6_epen, pathway %in% c(epen_6months_down))

p2 <- glut_6_epen_down_sub %>%
  mutate(pathway=factor(pathway, levels=epen_6months_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 6 mo Ependymal Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-3, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

glut_12_epen_up_sub <- subset(glut_12_epen, pathway %in% c(epen_12months_up))

p3 <- glut_12_epen_up_sub %>%
  mutate(pathway=factor(pathway, levels=epen_12months_up)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Ependymal Upregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(0, 4) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")

glut_12_epen_down_sub <- subset(glut_12_epen, pathway %in% c(epen_12months_down))

p4 <- glut_12_epen_down_sub %>%
  mutate(pathway=factor(pathway, levels=epen_12months_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Ependymal Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10)) +
	xlim(-7, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red")


#Standardize size
p_all <- align_plots(p1, p2, p3, p4, align="v")

ggsave2("Ependymal_6_gsea_glut_upreg_sig_venncol_0.07.png", p_all[[1]], dpi = 500, width = 7, height = 3.5)
ggsave2("Ependymal_6_gsea_glut_downreg_sig_venncol_0.07.png", p_all[[2]], dpi = 500, width = 7, height = 2.3)
ggsave2("Ependymal_12_gsea_glut_upreg_sig_venncol_0.07.png", p_all[[3]], dpi = 500, width = 7, height = 5)
ggsave2("Ependymal_12_gsea_glut_downreg_sig_venncol_0.07.png", p_all[[4]], dpi = 500, width = 7, height = 6)

#Just for bottom legend
png(file = "Ependymal_12_gsea_glut_downreg_sig_venncol_legend_bot_padj_0.07.png", width = 7, height = 6, units = "in", res = 500)
p5 <- glut_12_epen_down_sub %>%
  mutate(pathway=factor(pathway, levels=epen_12months_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Ependymal Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10), legend.position = "bottom", legend.box.just = "left") +
	xlim(-5, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red") +
	guides(colour = guide_colourbar(order = 1, barwidth = 10), size = guide_legend(order = 2))
print(p5)
dev.off()

png(file = "Ependymal_12_gsea_glut_downreg_sig_venncol_legend_bot_count_0.07.png", width = 7, height = 6, units = "in", res = 500)
p5 <- glut_12_epen_down_sub %>%
  mutate(pathway=factor(pathway, levels=epen_12months_down)) %>%
  mutate(pathway=forcats::fct_relabel(pathway, function(x) stringr::str_replace_all(x, "_", " ")), pathway = str_wrap_factor(pathway, 35)) %>% # Only x characters per line   
ggplot(aes(x=NES, y=pathway, size=count, color=padj)) +
	geom_point() +
	scale_color_continuous(low="red", high="blue", name = "padj", guide=guide_colorbar(reverse=TRUE, show.limits = TRUE), limits=c(0,0.07)) +
	ylab(NULL) + ggtitle("Glut1KO 12 mo Ependymal Downregulated GO Biological Process Terms") + scale_size(range=c(3, 8)) +
	theme(plot.title = element_text(hjust=1, size = 10), legend.position = "bottom", legend.box.just = "left") +
	xlim(-5, 0) +
	scale_size_continuous(limits=c(1,400)) +
	geom_vline(xintercept=c(0), linetype="dashed", color = "red") +
	guides(colour = guide_colourbar(order = 2, barwidth = 10), size = guide_legend(order = 1))
print(p5)
dev.off()
