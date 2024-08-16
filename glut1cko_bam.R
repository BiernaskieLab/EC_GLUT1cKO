#SLURM
salloc --mem=32G --cpus-per-task=16 --nodes=1 --ntasks=1 --time=05:00:00

source ~/software/init-conda.sh
export PATH="/home/alexander.pun/software/miniconda3/bin:$PATH"
conda activate seurat4


#Check Glut1 cKO Bam-----------------------------------------------------------------------------------------------------
#6 months----------------------------
R
library(Seurat)

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
load("GLUTKO.Robj", verbose = TRUE)

Idents(GLUTKO) <- "cell_ident"
GLUTKO_epen = subset(GLUTKO, idents = "Ependymal")
Idents(GLUTKO_epen) <- "sample_id"
N_6K_epen = subset(GLUTKO_epen, idents = "N_6K")

N_6K_epen_umi_id <- rownames(N_6K_epen@meta.data)
N_6K_epen_umi_id_1 <- gsub("-10$", "-1", N_6K_epen_umi_id)
N_6K_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_6K_epen_umi_id_1, Group = "Ependymal_6K")
write.table(N_6K_epen_umi_id_1_df, file = "N_6K_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

N_6W_epen = subset(GLUTKO_epen, idents = "N_6W")
N_6W_epen_umi_id <- rownames(N_6W_epen@meta.data)
N_6W_epen_umi_id_1 <- gsub("-9$", "-1", N_6W_epen_umi_id)
N_6W_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_6W_epen_umi_id_1, Group = "Ependymal_6W")
write.table(N_6W_epen_umi_id_1_df, file = "N_6W_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

q()

#Filter
source activate sinto
cd /work/biernaskie_lab/apun/nilesh_glut1/seurat

sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6K_count_chem/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_6K_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6K_count_chem/
sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6W_count_chem/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_6W_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6W_count_chem/

#Index
source activate velocyto
cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun

cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6K_count_chem/outs
samtools index Ependymal_6K.bam

cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_6W_count_chem/outs
samtools index Ependymal_6W.bam


#12 months ----------------------------
source activate seurat4

R
library(Seurat)
#Control N3, N7

setwd("/work/biernaskie_lab/apun/nilesh_glut1/seurat")
load("GLUTKO.Robj", verbose = TRUE)
Idents(GLUTKO) <- "cell_ident"
GLUTKO_epen = subset(GLUTKO, idents = "Ependymal")

Idents(GLUTKO_epen) <- "sample_id"

N_3_epen = subset(GLUTKO_epen, idents = "N_3")
N_3_epen_umi_id <- rownames(N_3_epen@meta.data)
N_3_epen_umi_id_1 <- gsub("-3$", "-1", N_3_epen_umi_id)
N_3_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_3_epen_umi_id_1, Group = "Ependymal_N_3")
write.table(N_3_epen_umi_id_1_df, file = "N_3_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

N_7_epen = subset(GLUTKO_epen, idents = "N_7")
N_7_epen_umi_id <- rownames(N_7_epen@meta.data)
N_7_epen_umi_id_1 <- gsub("-7$", "-1", N_7_epen_umi_id)
N_7_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_7_epen_umi_id_1, Group = "Ependymal_N_7")
write.table(N_7_epen_umi_id_1_df, file = "N_7_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#KO N4, N8
N_4_epen = subset(GLUTKO_epen, idents = "N_4")
N_4_epen_umi_id <- rownames(N_4_epen@meta.data)
N_4_epen_umi_id_1 <- gsub("-4$", "-1", N_4_epen_umi_id)
N_4_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_4_epen_umi_id_1, Group = "Ependymal_N_4")
write.table(N_4_epen_umi_id_1_df, file = "N_4_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

N_8_epen = subset(GLUTKO_epen, idents = "N_8")
N_8_epen_umi_id <- rownames(N_8_epen@meta.data)
N_8_epen_umi_id_1 <- gsub("-8$", "-1", N_8_epen_umi_id)
N_8_epen_umi_id_1_df <- data.frame(Cell_Barcode = N_8_epen_umi_id_1, Group = "Ependymal_N_8")
write.table(N_8_epen_umi_id_1_df, file = "N_8_epen_umi_id_1_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

q()

#Filter
source activate sinto
cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam

sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_7_chem_deep/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_7_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam
sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_3_chem/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_3_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam

sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_8_chem_deep/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_8_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam
sinto filterbarcodes -p 30 -b /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/N_4_chem/outs/possorted_genome_bam.bam -c /work/biernaskie_lab/apun/nilesh_glut1/seurat/N_4_epen_umi_id_1_df.txt --outdir /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam


#Merge
source activate velocyto
cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam
samtools merge N3_7_ependymal_possorted_genome_bam.bam Ependymal_N_3.bam Ependymal_N_7.bam
samtools merge N4_8_ependymal_possorted_genome_bam.bam Ependymal_N_4.bam Ependymal_N_8.bam

#Index
source activate velocyto
cd /work/biernaskie_lab/sarthak_sinha/Sequencing_Run_Feb_2024/counts/GLUT1/Deeprun/12mo_bam
samtools index N3_7_ependymal_possorted_genome_bam.bam
samtools index N4_8_ependymal_possorted_genome_bam.bam



#samtools view -b possorted_genome_bam.bam "chr4:119,131,976-119,133,887" > possorted_genome_bam_glut1_6K.bam