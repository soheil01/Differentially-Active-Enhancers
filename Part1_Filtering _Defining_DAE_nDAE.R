
###
### loading packages
###
suppressPackageStartupMessages ({library(Biostrings)
  library(stringr)
  library(pvclust)
  library(data.table)
  library(plyr)
  library(ggplot2)
  library(scales)
  library(EDASeq)
  library(reshape2)
  library(ggrepel)
  library(gplots)
  library(edgeR)
  library(RColorBrewer)
  library(ggpubr)
  library(ggbiplot)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  require(gridExtra)
  library(GenomicRanges)
  library(devtools)
  library(ggsignif)
  library(forestplot)
  library(phastCons100way.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19)})


######
###### pCR table is created using "pCR.txt". Counting reads for each sample are applied using `multiBamCov` commands.
######

# Read pCR table

## Data including All pCR regions as rows and epigenome data as columns
pCR <- read.delim("pCR_Raw.Counts.Regions.bed")
## Filter1: Remove enhancers with zero counts across all samples.
pCR_nonzero <- pCR[!rowSums(pCR[,-c(1:3)]) == 0,]
## Filter2: Select enhancers with length between 50-1000 bps.
pCR_nonzero <- pCR_nonzero %>% 
  mutate (Length = end - start)
pCR_NonZero.length.filtered <- pCR_nonzero[(pCR_nonzero$Length >= 50 & pCR_nonzero$Length<=1000), ]
dim(pCR_NonZero.length.filtered)

ggplot(pCR_NonZero.length.filtered, aes(x=Length)) + 
  geom_histogram(stat='count', binwidth=100, color="blue", fill="blue") +
  theme_classic()+
  scale_x_continuous("Length(bp)", expand = c(0, 0)) +
  scale_y_continuous("Frequency", expand = c(0, 0)) + 
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'))

## Filter3: Remove enhancers that are located in 2kb UP- and 1kb Down-stream of TSS (Downloaded from Ensembl). 
### Run it in bash using bedtools intersect command and "-v" flag.
#### write.table(pCR_NonZero.length.filtered,'pCA_Fetal_Raw_filtered.bed', quote=F, sep="\t", row.names= F, col.names=T)

######
###### Defining DAE vs nDAE
######

# Read raw count data and make a coordinate column
## Data including All filtered pCR regions as rows and epigenome data as columns
pCR_raw <- read.delim("pCR_Fetal_Raw_filtered.bed")
rownames(pCR_raw) <- paste0(pCR_raw$chr, ":", pCR_raw$start, "-", pCR_raw$end)

# Both DAE and nDAE are defined for each epigenome mark across time-points and brain parts.

##ATAC
### ATAC_Time-points
ATAC_design <- read.csv("Design_Matrix/ATAC_design.csv", row.names = 1)# Read design table
ATAC_data <- pCR_raw[,colnames(pCR_raw) %in% rownames(ATAC_design)]# select ATAC-seq samples
ATAC_data <- ATAC_data[rowSums(ATAC_data >=10) >=2,]# Filter low count data (< 10 reads) in at least two samples
all(rownames(ATAC_design) %in% colnames(ATAC_data)); all(rownames(ATAC_design) == colnames(ATAC_data))# check samples and their order between count and design tables.
ATAC_design_times <- model.matrix(~0 + Group_Time, ATAC_design)# making model matrix

ATAC_list_timepoints <- DGEList(counts = ATAC_data, group=ATAC_design$Group_Time)# making DGE list 
ATAC_norm_timepoints <- calcNormFactors(ATAC_list_timepoints)# Normalizing data
ATAC_dispersion_timepoints <- estimateDisp(ATAC_norm_timepoints, ATAC_design_times, robust=TRUE)# Estimating distribution
ATAC_fited_timepoints <- glmQLFit(ATAC_dispersion_timepoints, ATAC_design_times, robust=TRUE)# Fitting model

ATAC_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB,
  T1vsT3 = Group_TimeA - Group_TimeC,
  T2vsT3 = Group_TimeB - Group_TimeC, levels=ATAC_design_times)# making contrast table

qlf_ATAC_timepoints <- glmQLFTest(ATAC_fited_timepoints, contrast=ATAC_Con_times)# Calculating Differentially active enhancers based on fitted model
qlf_ATAC_timepoints <- data.frame(topTags(qlf_ATAC_timepoints, 3000000))
DAE_ATAC_timepoints <- qlf_ATAC_timepoints[qlf_ATAC_timepoints$FDR < 0.05,]## significant regions
dim(DAE_ATAC_timepoints)

### ATAC_Brain parts
ATAC_design_Types <- model.matrix(~0 + Group_Type, ATAC_design)
ATAC_list_Types <- DGEList(counts = ATAC_data, group=ATAC_design$Group_Type)
ATAC_norm_Types <- calcNormFactors(ATAC_list_Types)
ATAC_dispersion_Types <- estimateDisp(ATAC_norm_Types, ATAC_design_Types, robust=TRUE)
ATAC_fited_Types <- glmQLFit(ATAC_dispersion_Types, ATAC_design_Types, robust=TRUE)
ATAC_Con_Types <- makeContrasts(T1vsT2 = Group_TypeA - Group_TypeB, levels=ATAC_design_Types)
qlf_ATAC_Types <- glmQLFTest(ATAC_fited_Types, contrast= ATAC_Con_Types)
qlf_ATAC_Types <- data.frame(topTags(qlf_ATAC_Types, 3000000))
DAE_ATAC_Types <- qlf_ATAC_Types[qlf_ATAC_Types$FDR < 0.05, ]
dim(DAE_ATAC_Types)

DAE_ATAC_timepoints <- data.frame(row.names(DAE_ATAC_timepoints)); colnames(DAE_ATAC_timepoints) <- "enhancer"
DAE_ATAC_Types <- data.frame(row.names(DAE_ATAC_Types)); colnames(DAE_ATAC_Types) <- "enhancer"  
DAE_ATAC <- unique(rbind(DAE_ATAC_timepoints[1], DAE_ATAC_Types[1]))
dim(DAE_ATAC)

## DNase
### DNase_Time-points
DNase_design <- read.csv("Design_Matrix/DNase_design.csv", row.names = 1)
DNase_data <- pCR_raw[,colnames(pCR_raw) %in% rownames(DNase_design)]
DNase_data <- DNase_data[rowSums(DNase_data >=10) >=2,]
all(rownames(DNase_design) %in% colnames(DNase_data)); all(rownames(DNase_design) == colnames(DNase_data))
DNase_design_times <- model.matrix(~0 + Group, DNase_design)

DNase_list_timepoints <- DGEList(counts = DNase_data, group=DNase_design$Group)
DNase_norm_timepoints <- calcNormFactors(DNase_list_timepoints)
DNase_dispersion_timepoints <- estimateDisp(DNase_norm_timepoints, DNase_design_times, robust=TRUE)
DNase_fited_timepoints <- glmQLFit(DNase_dispersion_timepoints, DNase_design_times, robust=TRUE)

DNase_Con_times <- makeContrasts(
  T1vsT2 = GroupA - GroupB, T1vsT3 = GroupA - GroupC, T1vsT4 = GroupA - GroupD, T1vsT5 = GroupA - GroupG,
  T1vsT6 = GroupA - GroupH, T1vsT7 = GroupA - GroupK, T2vsT3 = GroupB - GroupC, T2vsT4 = GroupB - GroupD,
  T2vsT5 = GroupB - GroupG, T2vsT6 = GroupB - GroupH, T2vsT7 = GroupB - GroupK, T3vsT4 = GroupC - GroupD,
  T3vsT5 = GroupC - GroupG, T3vsT6 = GroupC - GroupH, T3vsT7 = GroupC - GroupK, T4vsT5 = GroupD - GroupG,
  T4vsT6 = GroupD - GroupH, T4vsT7 = GroupD - GroupK, T5vsT6 = GroupG - GroupH, T5vsT7 = GroupG - GroupK,
  T6vsT7 = GroupH - GroupK, levels=DNase_design_times)

qlf_DNase_timepoints <- glmQLFTest(DNase_fited_timepoints, contrast=DNase_Con_times)
qlf_DNase_timepoints <- data.frame(topTags(qlf_DNase_timepoints, 3000000))
DAE_DNase <- qlf_DNase_timepoints[which(qlf_DNase_timepoints$FDR < 0.05), ]
DAE_DNase <- data.frame(row.names(DAE_DNase)); colnames(DAE_DNase) <- "enhancer"  
dim(DAE_DNase)

### DNase_Brain parts
#### There is no brain parts for DNase.

## H3K27a
### H3K27ac_Time-points
K27ac_design <- read.csv("Design_Matrix/H3K27ac_design.csv", row.names = 1)
K27ac_data <- pCR_raw[,colnames(pCR_raw) %in% rownames(K27ac_design)]
K27ac_data <- K27ac_data[rowSums(K27ac_data >=10) >=2,]
all(rownames(K27ac_design) %in% colnames(K27ac_data)); all(rownames(K27ac_design) == colnames(K27ac_data))
K27ac_design_times <- model.matrix (~0 + Group_Time, K27ac_design)

K27ac_list_timepoints <- DGEList(counts = K27ac_data, group=K27ac_design$Group_Time)
K27ac_norm_timepoints <- calcNormFactors(K27ac_list_timepoints)
K27ac_dispersion_timepoints <- estimateDisp(K27ac_norm_timepoints, K27ac_design_times, robust=TRUE)
K27ac_fited_timepoints <- glmQLFit(K27ac_dispersion_timepoints, K27ac_design_times, robust=TRUE)

K27ac_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB, 
  T1vsT3 = Group_TimeA - Group_TimeC,
  T1vsT4 = Group_TimeA - Group_TimeD,
  T1vsT5 = Group_TimeA - Group_TimeE,
  T1vsT6 = Group_TimeA - Group_TimeF,
  T1vsT7 = Group_TimeA - Group_TimeG,
  T1vsT8 = Group_TimeA - Group_TimeH,
  T1vsT9 = Group_TimeA - Group_TimeI,
  T1vsT10 = Group_TimeA - Group_TimeJ,
  T1vsT11 = Group_TimeA - Group_TimeK,
  T2vsT3 = Group_TimeB - Group_TimeC,
  T2vsT4 = Group_TimeB - Group_TimeD,
  T2vsT5 = Group_TimeB - Group_TimeE,
  T2vsT6 = Group_TimeB - Group_TimeF,
  T2vsT7 = Group_TimeB - Group_TimeG,
  T2vsT8 = Group_TimeB - Group_TimeH,
  T2vsT9 = Group_TimeB - Group_TimeI,
  T2vsT10 = Group_TimeB - Group_TimeJ,
  T2vsT11 = Group_TimeB - Group_TimeK,
  T3vsT4 = Group_TimeC - Group_TimeD,
  T3vsT5 = Group_TimeC - Group_TimeE,
  T3vsT6 = Group_TimeC - Group_TimeF,
  T3vsT7 = Group_TimeC - Group_TimeG,
  T3vsT8 = Group_TimeC - Group_TimeH,
  T3vsT9 = Group_TimeC - Group_TimeI,
  T3vsT10 = Group_TimeC - Group_TimeJ,
  T3vsT11 = Group_TimeC - Group_TimeK,
  T4vsT5 = Group_TimeD - Group_TimeE,
  T4vsT6 = Group_TimeD - Group_TimeF,
  T4vsT7 = Group_TimeD - Group_TimeG,
  T4vsT8 = Group_TimeD - Group_TimeH,
  T4vsT9 = Group_TimeD - Group_TimeI,
  T4vsT10 = Group_TimeD - Group_TimeJ,
  T4vsT11 = Group_TimeD - Group_TimeK,
  T5vsT6 = Group_TimeE - Group_TimeF,
  T5vsT7 = Group_TimeE - Group_TimeG,
  T5vsT8 = Group_TimeE - Group_TimeH,
  T5vsT9 = Group_TimeE - Group_TimeI,
  T5vsT10 = Group_TimeE - Group_TimeJ,
  T5vsT11 = Group_TimeE - Group_TimeK,
  T6vsT7 = Group_TimeF - Group_TimeG,
  T6vsT8 = Group_TimeF - Group_TimeH,
  T6vsT9 = Group_TimeF - Group_TimeI,
  T6vsT10 = Group_TimeF - Group_TimeJ,
  T6vsT11 = Group_TimeF - Group_TimeK,
  T7vsT8 = Group_TimeG - Group_TimeH,
  T7vsT9 = Group_TimeG - Group_TimeI,
  T7vsT10 = Group_TimeG - Group_TimeJ,
  T7vsT11 = Group_TimeG - Group_TimeK,
  T8vsT9 = Group_TimeH - Group_TimeI,
  T8vsT10 = Group_TimeH - Group_TimeJ,
  T8vsT11 = Group_TimeH - Group_TimeK,
  T9vsT10 = Group_TimeI - Group_TimeJ,
  T9vsT11 = Group_TimeI - Group_TimeK,
  T10vsT11 = Group_TimeJ - Group_TimeK,
  levels=K27ac_design_times)

qlf_K27ac_timepoints <- glmQLFTest(K27ac_fited_timepoints, contrast=K27ac_Con_times)
qlf_K27ac_timepoints <- data.frame(topTags(qlf_K27ac_timepoints, 3000000))
DAE_K27ac_timepoints <- qlf_K27ac_timepoints[which(qlf_K27ac_timepoints$FDR < 0.05), ]
dim(DAE_K27ac_timepoints)

### H3K27ac_Brain parts
K27ac_design_Type <- model.matrix(~0 + Group_Type, K27ac_design)
K27ac_list_Type <- DGEList(counts = K27ac_data, group=K27ac_design$Group_Type)
K27ac_norm_Type <- calcNormFactors(K27ac_list_Type)
K27ac_dispersion_Type <- estimateDisp(K27ac_norm_Type, K27ac_design_Type, robust=TRUE)
K27ac_fited_Type <- glmQLFit(K27ac_dispersion_Type, K27ac_design_Type, robust=TRUE)

K27ac_Con_Type <- makeContrasts(
  T1vsT2 = Group_TypeA - Group_TypeB,
  T1vsT3 = Group_TypeA - Group_TypeC,
  T1vsT4 = Group_TypeA - Group_TypeD,
  T1vsT5 = Group_TypeA - Group_TypeE,
  T1vsT6 = Group_TypeA - Group_TypeF,
  T1vsT7 = Group_TypeA - Group_TypeG,
  T1vsT8 = Group_TypeA - Group_TypeH,
  T1vsT9 = Group_TypeA - Group_TypeI,
  T2vsT3 = Group_TypeB - Group_TypeC,
  T2vsT4 = Group_TypeB - Group_TypeD,
  T2vsT5 = Group_TypeB - Group_TypeE,
  T2vsT6 = Group_TypeB - Group_TypeF,
  T2vsT7 = Group_TypeB - Group_TypeG,
  T2vsT8 = Group_TypeB - Group_TypeH,
  T2vsT9 = Group_TypeB - Group_TypeI,
  T3vsT4 = Group_TypeC - Group_TypeD,
  T3vsT5 = Group_TypeC - Group_TypeE,
  T3vsT6 = Group_TypeC - Group_TypeF,
  T3vsT7 = Group_TypeC - Group_TypeG,
  T3vsT8 = Group_TypeC - Group_TypeH,
  T3vsT9 = Group_TypeC - Group_TypeI,
  T4vsT5 = Group_TypeD - Group_TypeE,
  T4vsT6 = Group_TypeD - Group_TypeF,
  T4vsT7 = Group_TypeD - Group_TypeG,
  T4vsT8 = Group_TypeD - Group_TypeH,
  T4vsT9 = Group_TypeD - Group_TypeI,
  T5vsT6 = Group_TypeE - Group_TypeF,
  T5vsT7 = Group_TypeE - Group_TypeG,
  T5vsT8 = Group_TypeE - Group_TypeH,
  T5vsT9 = Group_TypeE - Group_TypeI,
  T6vsT7 = Group_TypeF - Group_TypeG,
  T6vsT8 = Group_TypeF - Group_TypeH,
  T6vsT9 = Group_TypeF - Group_TypeI,
  T7vsT8 = Group_TypeG - Group_TypeH,
  T7vsT9 = Group_TypeG - Group_TypeI,
  T8vsT9 = Group_TypeH - Group_TypeI,
  levels=K27ac_design_Type)

qlf_K27ac_Type <- glmQLFTest(K27ac_fited_Type, contrast=K27ac_Con_Type)
qlf_K27ac_Type <- data.frame(topTags(qlf_K27ac_Type, 3000000))
DAE_K27ac_Type <- qlf_K27ac_Type[which(qlf_K27ac_Type$FDR < 0.05), ]
dim(DAE_K27ac_Type)

K27ac_Con_Type.Reilly <- makeContrasts(T1vsT2 = Group_TypeF - Group_TypeG,levels=K27ac_design_Type)
qlf_K27ac_Type.Reilly <- glmQLFTest(K27ac_fited_Type, contrast=K27ac_Con_Type.Reilly)
qlf_K27ac_Type.Reilly <- data.frame(topTags(qlf_K27ac_Type.Reilly, 3000000))
DAE_K27ac_Type.Reilly <- qlf_K27ac_Type.Reilly[which(qlf_K27ac_Type.Reilly$FDR < 0.05), ]
dim(DAE_K27ac_Type.Reilly)

K27ac_Con_Type.PSychENCODE <- makeContrasts(T1vsT2 = Group_TypeA - Group_TypeB, levels=K27ac_design_Type)
qlf_K27ac_Type.PSychENCODE <- glmQLFTest(K27ac_fited_Type, contrast=K27ac_Con_Type.PSychENCODE)
qlf_K27ac_Type.PSychENCODE <- data.frame(topTags(qlf_K27ac_Type.PSychENCODE, 3000000))
DAE_K27ac_Type.PSychENCODE <- qlf_K27ac_Type.PSychENCODE[which(qlf_K27ac_Type.PSychENCODE$FDR < 0.05), ]
dim(DAE_K27ac_Type.PSychENCODE)

K27ac_Con_Type.Li <- makeContrasts(T1vsT2 = Group_TypeC - Group_TypeD,levels=K27ac_design_Type)
qlf_K27ac_Type.Li <- glmQLFTest(K27ac_fited_Type, contrast=K27ac_Con_Type.Li)
qlf_K27ac_Type.Li <- data.frame(topTags(qlf_K27ac_Type.Li, 3000000))
DAE_K27ac_Type.Li <- qlf_K27ac_Type.Li[which(qlf_K27ac_Type.Li$FDR < 0.05), ]
dim(DAE_K27ac_Type.Li)

table(rownames(DAE_K27ac_Type.Li) %in% rownames(DAE_K27ac_Type))
table(rownames(DAE_K27ac_Type) %in% rownames(DAE_K27ac_timepoints))

DAE_K27ac_timepoints <- data.frame(row.names(DAE_K27ac_timepoints)); colnames(DAE_K27ac_timepoints) <- "enhancer"
DAE_K27ac_Type <- data.frame(row.names(DAE_K27ac_Type)); colnames(DAE_K27ac_Type) <- "enhancer"  
DAE_K27ac_Type.Li <- data.frame(row.names(DAE_K27ac_Type.Li)); colnames(DAE_K27ac_Type.Li) <- "enhancer"  
DAE_K27ac <- unique(rbind(DAE_K27ac_timepoints[1], DAE_K27ac_Type[1], DAE_K27ac_Type.Li[1]))
dim(DAE_K27ac)

## H3K4me1
### H3K4me1_Time-points
K4me1_design <- read.csv("Design_Matrix/H3K4me1_design.csv", row.names = 1)
K4me1_data <- pCR_raw[, colnames(pCR_raw) %in% rownames(K4me1_design)]
K4me1_data <- K4me1_data[rowSums(K4me1_data >=10) >=2, ]
all(rownames(K4me1_design) %in% colnames(K4me1_data)); all(rownames(K4me1_design) == colnames(K4me1_data))
K4me1_design_times <- model.matrix(~0 + Group_Time, K4me1_design)

K4me1_list_timepoints <- DGEList(counts = K4me1_data, group=K4me1_design$Group_Time)
K4me1_norm_timepoints <- calcNormFactors(K4me1_list_timepoints)
K4me1_dispersion_timepoints <- estimateDisp(K4me1_norm_timepoints, K4me1_design_times, robust=TRUE)
K4me1_fited_timepoints <- glmQLFit(K4me1_dispersion_timepoints, K4me1_design_times, robust=TRUE)

K4me1_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB,
  T1vsT3 = Group_TimeA - Group_TimeC,
  T1vsT4 = Group_TimeA - Group_TimeD,
  T2vsT3 = Group_TimeB - Group_TimeC,
  T2vsT4 = Group_TimeB - Group_TimeD,
  T3vsT4 = Group_TimeC - Group_TimeD,
  levels=K4me1_design_times)

qlf_K4me1_timepoints <- glmQLFTest(K4me1_fited_timepoints, contrast=K4me1_Con_times)
qlf_K4me1_timepoints <- data.frame(topTags(qlf_K4me1_timepoints, 3000000))
DAE_K4me1_timepoints <- qlf_K4me1_timepoints[which(qlf_K4me1_timepoints$FDR < 0.05), ]
dim(DAE_K4me1_timepoints)

### H3K4me1_Brain parts
K4me1_design_Type <- model.matrix(~0 + Group_Type, K4me1_design)
K4me1_list_Type <- DGEList(counts = K4me1_data, group=K4me1_design$Group_Type)
K4me1_norm_Type <- calcNormFactors(K4me1_list_Type)
K4me1_dispersion_Type <- estimateDisp(K4me1_norm_Type, K4me1_design_Type, robust=TRUE)
K4me1_fited_Type <- glmQLFit(K4me1_dispersion_Type, K4me1_design_Type, robust=TRUE)

K4me1_Con_Type <- makeContrasts(
  T1vsT2 = Group_TypeA - Group_TypeB,
  T1vsT3 = Group_TypeA - Group_TypeC,
  T1vsT4 = Group_TypeA - Group_TypeD,
  T2vsT3 = Group_TypeB - Group_TypeC,
  T2vsT4 = Group_TypeB - Group_TypeD,
  T3vsT4 = Group_TypeC - Group_TypeD,
  levels=K4me1_design_Type)

qlf_K4me1_Type <- glmQLFTest(K4me1_fited_Type, contrast=K4me1_Con_Type)
qlf_K4me1_Type <- data.frame(topTags(qlf_K4me1_Type, 3000000))
DAE_K4me1_Type <- qlf_K4me1_Type[which(qlf_K4me1_Type$FDR < 0.05),]
dim(unique(DAE_K4me1_Type[1]))

K4me1_Con_Type.CNSPC.GENSPC <- makeContrasts(T1vsT2 = Group_TypeB - Group_TypeC,levels=K4me1_design_Type)
qlf_K4me1_Type.CNSPC.GENSPC <- glmQLFTest(K4me1_fited_Type, contrast=K4me1_Con_Type.CNSPC.GENSPC)
qlf_K4me1_Type.CNSPC.GENSPC <- data.frame(topTags(qlf_K4me1_Type.CNSPC.GENSPC, 3000000))
DAE_K4me1_Type.CNSPC.GENSPC <- qlf_K4me1_Type.CNSPC.GENSPC[which(qlf_K4me1_Type.CNSPC.GENSPC$FDR < 0.05),]
dim(DAE_K4me1_Type.CNSPC.GENSPC)

table(rownames(DAE_K4me1_Type) %in% rownames(DAE_K4me1_timepoints))

DAE_K4me1_timepoints <- data.frame(row.names(DAE_K4me1_timepoints)); colnames(DAE_K4me1_timepoints) <- "enhancer"  
DAE_K4me1_Type <- data.frame(row.names(DAE_K4me1_Type)); colnames(DAE_K4me1_Type) <- "enhancer"  
DAE_K4me1 <- unique(rbind(DAE_K4me1_timepoints[1], DAE_K4me1_Type[1]))
dim(DAE_K4me1)

## H3K4me2
### H3K4me2_Time points
K4me2_design <- read.csv("Design_Matrix/H3K4me2_design.csv", row.names = 1)
K4me2_data <- pCR_raw[, colnames(pCR_raw) %in% rownames(K4me2_design)]
K4me2_data <- K4me2_data[rowSums(K4me2_data >=10) >=2, ]
all(rownames(K4me2_design) %in% colnames(K4me2_data)); all(rownames(K4me2_design) == colnames(K4me2_data))
K4me2_design_times <- model.matrix(~0 + Group_Time, K4me2_design)

K4me2_list_timepoints <- DGEList(counts = K4me2_data, group=K4me2_design$Group_Time)
K4me2_norm_timepoints <- calcNormFactors(K4me2_list_timepoints)
K4me2_dispersion_timepoints <- estimateDisp(K4me2_norm_timepoints, K4me2_design_times, robust=TRUE)
K4me2_fited_timepoints <- glmQLFit(K4me2_dispersion_timepoints, K4me2_design_times, robust=TRUE)

K4me2_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB,
  T1vsT3 = Group_TimeA - Group_TimeC,
  T2vsT3 = Group_TimeB - Group_TimeC, levels=K4me2_design_times)

qlf_K4me2_timepoints <- glmQLFTest(K4me2_fited_timepoints, contrast=K4me2_Con_times)
qlf_K4me2_timepoints <- data.frame(topTags(qlf_K4me2_timepoints, 3000000))
DAE_K4me2_timepoints <- qlf_K4me2_timepoints[which(qlf_K4me2_timepoints$FDR < 0.05), ]
dim(DAE_K4me2_timepoints)

### H3K4me2_Brain parts
K4me2_design_Type <- model.matrix(~0 + Group_Type, K4me2_design)
K4me2_list_Type <- DGEList(counts = K4me2_data, group=K4me2_design$Group_Type)
K4me2_norm_Type <- calcNormFactors(K4me2_list_Type)
K4me2_dispersion_Type <- estimateDisp(K4me2_norm_Type, K4me2_design_Type, robust=TRUE)
K4me2_fited_Type <- glmQLFit(K4me2_dispersion_Type, K4me2_design_Type, robust=TRUE)

K4me2_Con_Type <- makeContrasts(
  T1vsT2 = Group_TypeA - Group_TypeB,
  T1vsT3 = Group_TypeA - Group_TypeC,
  T2vsT3 = Group_TypeB - Group_TypeC, levels=K4me2_design_Type)

qlf_K4me2_Type <- glmQLFTest(K4me2_fited_Type, contrast=K4me2_Con_Type)
qlf_K4me2_Type <- data.frame(topTags(qlf_K4me2_Type, 3000000))
DAE_K4me2_Type <- qlf_K4me2_Type[which(qlf_K4me2_Type$FDR < 0.05), ]
dim(DAE_K4me2_Type)

K4me2_Con_Type.O.F <- makeContrasts(T1vsT2 = Group_TypeB - Group_TypeC, levels=K4me2_design_Type)
qlf_K4me2_Type.O.F <- glmQLFTest(K4me2_fited_Type, contrast=K4me2_Con_Type.O.F)
qlf_K4me2_Type.O.F <- data.frame(topTags(qlf_K4me2_Type.O.F, 3000000))
DAE_K4me2_Type.O.F <- qlf_K4me2_Type.O.F[which(qlf_K4me2_Type.O.F$FDR < 0.05), ]
dim(DAE_K4me2_Type.O.F)

table(rownames(DAE_K4me2_Type) %in% rownames(DAE_K4me2_timepoints))

DAE_K4me2_timepoints <- data.frame(row.names(DAE_K4me2_timepoints)); colnames(DAE_K4me2_timepoints) <- "enhancer"  
DAE_K4me2_Type <- data.frame(row.names(DAE_K4me2_Type)); colnames(DAE_K4me2_Type) <- "enhancer"
DAE_K4me2 <- unique(rbind(DAE_K4me2_timepoints[1], DAE_K4me2_Type[1]))
dim(DAE_K4me2)

## H3K4me3
### H3K4me3_Time points
K4me3_design <- read.csv("Design_Matrix/H3K4me3_design.csv", row.names = 1)
K4me3_data <- pCR_raw[,colnames(pCR_raw) %in% rownames(K4me3_design)]
K4me3_data <- K4me3_data[rowSums(K4me3_data >=10) >=2, ]
all(rownames(K4me3_design) %in% colnames(K4me3_data)); all(rownames(K4me3_design) == colnames(K4me3_data))
K4me3_design_times <- model.matrix(~0 + Group_Time, K4me3_design)

K4me3_list_timepoints <- DGEList(counts = K4me3_data, group=K4me3_design$Group_Time)
K4me3_norm_timepoints <- calcNormFactors(K4me3_list_timepoints)
K4me3_dispersion_timepoints <- estimateDisp(K4me3_norm_timepoints, K4me3_design_times, robust=TRUE)
K4me3_fited_timepoints <- glmQLFit(K4me3_dispersion_timepoints, K4me3_design_times, robust=TRUE)

K4me3_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB,
  T1vsT3 = Group_TimeA - Group_TimeC,
  T1vsT4 = Group_TimeA - Group_TimeD,
  T1vsT5 = Group_TimeA - Group_TimeE,
  T1vsT6 = Group_TimeA - Group_TimeF,
  T1vsT7 = Group_TimeA - Group_TimeG,
  T1vsT8 = Group_TimeA - Group_TimeH,
  T2vsT3 = Group_TimeB - Group_TimeC,
  T2vsT4 = Group_TimeB - Group_TimeD,
  T2vsT5 = Group_TimeB - Group_TimeE,
  T2vsT6 = Group_TimeB - Group_TimeF,
  T2vsT7 = Group_TimeB - Group_TimeG,
  T2vsT8 = Group_TimeB - Group_TimeH,
  T3vsT4 = Group_TimeC - Group_TimeD,
  T3vsT5 = Group_TimeC - Group_TimeE,
  T3vsT6 = Group_TimeC - Group_TimeF,
  T3vsT7 = Group_TimeC - Group_TimeG,
  T3vsT8 = Group_TimeC - Group_TimeH,
  T4vsT5 = Group_TimeD - Group_TimeE,
  T4vsT6 = Group_TimeD - Group_TimeF,
  T4vsT7 = Group_TimeD - Group_TimeG,
  T4vsT8 = Group_TimeD - Group_TimeH,
  T5vsT6 = Group_TimeE - Group_TimeF,
  T5vsT7 = Group_TimeE - Group_TimeG,
  T5vsT8 = Group_TimeE - Group_TimeH,
  T6vsT7 = Group_TimeF - Group_TimeG,
  T6vsT8 = Group_TimeF - Group_TimeH,
  T7vsT8 = Group_TimeG - Group_TimeH,
  levels=K4me3_design_times)

qlf_K4me3_timepoints <- glmQLFTest(K4me3_fited_timepoints, contrast=K4me3_Con_times)
qlf_K4me3_timepoints <- data.frame(topTags(qlf_K4me3_timepoints, 3000000))
DAE_K4me3_timepoints <- qlf_K4me3_timepoints[which(qlf_K4me3_timepoints$FDR < 0.05), ]
dim(DAE_K4me3_timepoints)

### H3K4me3_Brain parts
K4me3_design_type <- model.matrix(~0 + Group_Type, K4me3_design)
K4me3_list_type <- DGEList(counts = K4me3_data, group=K4me3_design$Group_Type)
K4me3_norm_type <- calcNormFactors(K4me3_list_type)
K4me3_dispersion_type <- estimateDisp(K4me3_norm_type, K4me3_design_type, robust=TRUE)
K4me3_fited_type <- glmQLFit(K4me3_dispersion_type, K4me3_design_type, robust=TRUE)

K4me3_Con_type <- makeContrasts(
  T1vsT2 = Group_TypeA - Group_TypeB,
  T1vsT3 = Group_TypeA - Group_TypeC,
  T1vsT4 = Group_TypeA - Group_TypeD,
  T1vsT5 = Group_TypeA - Group_TypeE,
  T1vsT6 = Group_TypeA - Group_TypeF,
  T1vsT7 = Group_TypeA - Group_TypeG,
  T1vsT8 = Group_TypeA - Group_TypeH,
  T1vsT9 = Group_TypeA - Group_TypeI,
  T2vsT3 = Group_TypeB - Group_TypeC,
  T2vsT4 = Group_TypeB - Group_TypeD,
  T2vsT5 = Group_TypeB - Group_TypeE,
  T2vsT6 = Group_TypeB - Group_TypeF,
  T2vsT7 = Group_TypeB - Group_TypeG,
  T2vsT8 = Group_TypeB - Group_TypeH,
  T2vsT9 = Group_TypeB - Group_TypeI,
  T3vsT4 = Group_TypeC - Group_TypeD,
  T3vsT5 = Group_TypeC - Group_TypeE,
  T3vsT6 = Group_TypeC - Group_TypeF,
  T3vsT7 = Group_TypeC - Group_TypeG,
  T3vsT8 = Group_TypeC - Group_TypeH,
  T3vsT9 = Group_TypeC - Group_TypeI,
  T4vsT5 = Group_TypeD - Group_TypeE,
  T4vsT6 = Group_TypeD - Group_TypeF,
  T4vsT7 = Group_TypeD - Group_TypeG,
  T4vsT8 = Group_TypeD - Group_TypeH,
  T4vsT9 = Group_TypeD - Group_TypeI,
  T5vsT6 = Group_TypeE - Group_TypeF,
  T5vsT7 = Group_TypeE - Group_TypeG,
  T5vsT8 = Group_TypeE - Group_TypeH,
  T5vsT9 = Group_TypeE - Group_TypeI,
  T6vsT7 = Group_TypeF - Group_TypeG,
  T6vsT8 = Group_TypeF - Group_TypeH,
  T6vsT9 = Group_TypeF - Group_TypeI,
  T7vsT8 = Group_TypeG - Group_TypeH,
  T7vsT9 = Group_TypeG - Group_TypeI,
  T8vsT9 = Group_TypeH - Group_TypeI,
  levels=K4me3_design_type)

qlf_K4me3_type <- glmQLFTest(K4me3_fited_type, contrast=K4me3_Con_type)
qlf_K4me3_type <- data.frame(topTags(qlf_K4me3_type, 3000000))
DAE_K4me3_type <- qlf_K4me3_type[which(qlf_K4me3_type$FDR < 0.05), ]
dim(DAE_K4me3_type)

K4me3_Con_Type.PSychENCODE <- makeContrasts(T1vsT2 = Group_TypeA - Group_TypeB, levels=K4me3_design_type)
qlf_K4me3_Type.PSychENCODE <- glmQLFTest(K4me3_fited_type, contrast=K4me3_Con_Type.PSychENCODE)
qlf_K4me3_Type.PSychENCODE <- data.frame(topTags(qlf_K4me3_Type.PSychENCODE, 3000000))
DAE_K4me3_Type.PSychENCODE <- qlf_K4me3_Type.PSychENCODE[which(qlf_K4me3_Type.PSychENCODE$FDR < 0.05), ]
dim(DAE_K4me3_Type.PSychENCODE)

table(rownames(DAE_K4me3_type) %in% rownames(DAE_K4me3_timepoints))

DAE_K4me3_timepoints <- data.frame(row.names(DAE_K4me3_timepoints)); colnames(DAE_K4me3_timepoints) <- "enhancer"  
DAE_K4me3_type <- data.frame(row.names(DAE_K4me3_type)); colnames(DAE_K4me3_type) <- "enhancer"
DAE_K4me3 <- unique(rbind(DAE_K4me3_timepoints[1], DAE_K4me3_type[1]))
dim(DAE_K4me3)

## H3K27me3
### H3K27me3_Time-points
K27me3_design <- read.csv("Design_Matrix/H3K27me3_design.csv", row.names = 1)
K27me3_data <- pCR_raw[, colnames(pCR_raw) %in% rownames(K27me3_design)]
K27me3_data <- K27me3_data[rowSums(K27me3_data >=10) >=2, ]
all(rownames(K27me3_design) %in% colnames(K27me3_data)); all(rownames(K27me3_design) == colnames(K27me3_data))
K27me3_design_times <- model.matrix(~0 + Group_Time, K27me3_design)

K27me3_list_timepoints <- DGEList(counts = K27me3_data, group=K27me3_design$Group_Time)
K27me3_norm_timepoints <- calcNormFactors(K27me3_list_timepoints)
K27me3_dispersion_timepoints <- estimateDisp(K27me3_norm_timepoints, K27me3_design_times, robust=TRUE)
K27me3_fited_timepoints <- glmQLFit(K27me3_dispersion_timepoints, K27me3_design_times, robust=TRUE)

K27me3_Con_times <- makeContrasts(
  T1vsT2 = Group_TimeA - Group_TimeB,
  T1vsT3 = Group_TimeA - Group_TimeC,
  T1vsT4 = Group_TimeA - Group_TimeD,
  T1vsT5 = Group_TimeA - Group_TimeE,
  T1vsT6 = Group_TimeA - Group_TimeF,
  T1vsT7 = Group_TimeA - Group_TimeG,
  T1vsT8 = Group_TimeA - Group_TimeH,
  T2vsT3 = Group_TimeB - Group_TimeC,
  T2vsT4 = Group_TimeB - Group_TimeD,
  T2vsT5 = Group_TimeB - Group_TimeE,
  T2vsT6 = Group_TimeB - Group_TimeF,
  T2vsT7 = Group_TimeB - Group_TimeG,
  T2vsT8 = Group_TimeB - Group_TimeH,
  T3vsT4 = Group_TimeC - Group_TimeD,
  T3vsT5 = Group_TimeC - Group_TimeE,
  T3vsT6 = Group_TimeC - Group_TimeF,
  T3vsT7 = Group_TimeC - Group_TimeG,
  T3vsT8 = Group_TimeC - Group_TimeH,
  T4vsT5 = Group_TimeD - Group_TimeE,
  T4vsT6 = Group_TimeD - Group_TimeF,
  T4vsT7 = Group_TimeD - Group_TimeG,
  T4vsT8 = Group_TimeD - Group_TimeH,
  T5vsT6 = Group_TimeE - Group_TimeF,
  T5vsT7 = Group_TimeE - Group_TimeG,
  T5vsT8 = Group_TimeE - Group_TimeH,
  T6vsT7 = Group_TimeF - Group_TimeG,
  T6vsT8 = Group_TimeF - Group_TimeH,
  T7vsT8 = Group_TimeG - Group_TimeH,
  levels=K27me3_design_times)

qlf_K27me3_timepoints <- glmQLFTest(K27me3_fited_timepoints, contrast=K27me3_Con_times)
qlf_K27me3_timepoints <- data.frame(topTags(qlf_K27me3_timepoints, 3000000))
DAE_K27me3_timepoints <- qlf_K27me3_timepoints[which(qlf_K27me3_timepoints$FDR < 0.05), ]
dim(DAE_K27me3_timepoints)

### H3K27me3_Brain parts
K27me3_design_type <- model.matrix(~0 + Group_Type, K27me3_design)
K27me3_list_type <- DGEList(counts = K27me3_data, group=K27me3_design$Group_Type)
K27me3_norm_type <- calcNormFactors(K27me3_list_type)
K27me3_dispersion_type <- estimateDisp(K27me3_norm_type, K27me3_design_type, robust=TRUE)
K27me3_fited_type <- glmQLFit(K27me3_dispersion_type, K27me3_design_type, robust=TRUE)

K27me3_Con_type <- makeContrasts(
  T1vsT2 = Group_TypeA - Group_TypeB,
  T1vsT3 = Group_TypeA - Group_TypeC,
  T1vsT4 = Group_TypeA - Group_TypeD,
  T1vsT5 = Group_TypeA - Group_TypeE,
  T1vsT6 = Group_TypeA - Group_TypeF,
  T1vsT7 = Group_TypeA - Group_TypeG,
  T2vsT3 = Group_TypeB - Group_TypeC,
  T2vsT4 = Group_TypeB - Group_TypeD,
  T2vsT5 = Group_TypeB - Group_TypeE,
  T2vsT6 = Group_TypeB - Group_TypeF,
  T2vsT7 = Group_TypeB - Group_TypeG,
  T3vsT4 = Group_TypeC - Group_TypeD,
  T3vsT5 = Group_TypeC - Group_TypeE,
  T3vsT6 = Group_TypeC - Group_TypeF,
  T3vsT7 = Group_TypeC - Group_TypeG,
  T4vsT5 = Group_TypeD - Group_TypeE,
  T4vsT6 = Group_TypeD - Group_TypeF,
  T4vsT7 = Group_TypeD - Group_TypeG,
  T5vsT6 = Group_TypeE - Group_TypeF,
  T5vsT7 = Group_TypeE - Group_TypeG,
  T6vsT7 = Group_TypeF - Group_TypeG,
  levels=K27me3_design_type)

qlf_K27me3_type <- glmQLFTest(K27me3_fited_type, contrast=K27me3_Con_type)
qlf_K27me3_type <- data.frame(topTags(qlf_K27me3_type, 3000000))
DAE_K27me3_type <- qlf_K27me3_type[which(qlf_K27me3_type$FDR < 0.05), ]
dim(DAE_K27me3_type)

K27me3_Con_Type.PSychENCODE <- makeContrasts(T1vsT2 = Group_TypeA - Group_TypeB, levels=K27me3_design_type)
qlf_K27me3_Type.PSychENCODE <- glmQLFTest(K27me3_fited_type, contrast=K27me3_Con_Type.PSychENCODE)
qlf_K27me3_Type.PSychENCODE <- data.frame(topTags(qlf_K27me3_Type.PSychENCODE, 3000000))
dim(qlf_K27me3_Type.PSychENCODE)
DAE_K27me3_Type.PSychENCODE <- qlf_K27me3_Type.PSychENCODE[which(qlf_K27me3_Type.PSychENCODE$FDR < 0.05), ]
dim(DAE_K27me3_Type.PSychENCODE)

K27me3_Con_Type.CNSPC.GENSPC <- makeContrasts(T1vsT2 = Group_TypeC - Group_TypeD, levels=K27me3_design_type)
qlf_K27me3_Type.CNSPC.GENSPC <- glmQLFTest(K27me3_fited_type, contrast=K27me3_Con_Type.CNSPC.GENSPC)
qlf_K27me3_Type.CNSPC.GENSPC <- data.frame(topTags(qlf_K27me3_Type.CNSPC.GENSPC, 3000000))
DAE_K27me3_Type.CNSPC.GENSPC <- qlf_K27me3_Type.CNSPC.GENSPC[which(qlf_K27me3_Type.CNSPC.GENSPC$FDR < 0.05), ]
dim(DAE_K27me3_Type.CNSPC.GENSPC)

table(rownames(DAE_K27me3_Type.CNSPC.GENSPC) %in% rownames(DAE_K27me3_type))
table(rownames(DAE_K27me3_timepoints) %in% rownames(DAE_K27me3_type))

DAE_K27me3_timepoints <- data.frame(row.names(DAE_K27me3_timepoints)); colnames(DAE_K27me3_timepoints) <- "enhancer"  
DAE_K27me3_type <- data.frame(row.names(DAE_K27me3_type)); colnames(DAE_K27me3_type) <- "enhancer"
DAE_K27me3_Type.CNSPC.GENSPC <- data.frame(row.names(DAE_K27me3_Type.CNSPC.GENSPC)); colnames(DAE_K27me3_Type.CNSPC.GENSPC) <- "enhancer"
DAE_K27me3 <- unique(rbind(DAE_K27me3_timepoints[1], DAE_K27me3_type[1], DAE_K27me3_Type.CNSPC.GENSPC[1]))
dim(DAE_K27me3)

## Common DAEs in at least two epigenome marks

ATAC_DAE <- data.frame(DAE_ATAC [(DAE_ATAC$enhancer %in% DAE_DNase$enhancer|
                                    DAE_ATAC$enhancer %in% DAE_K27ac$enhancer|
                                    DAE_ATAC$enhancer %in% DAE_K4me1$enhancer|
                                    DAE_ATAC$enhancer %in% DAE_K4me2$enhancer|
                                    DAE_ATAC$enhancer %in% DAE_K4me3$enhancer|
                                    DAE_ATAC$enhancer %in% DAE_K27me3$enhancer),]) 
colnames(ATAC_DAE) <- 'Enhancers'
dim(ATAC_DAE)# Common ATC-DAE in at least two epigenome data 
dim(DAE_ATAC)# Number of initial defined DAE

DNase_DAE <- data.frame(DAE_DNase [(DAE_DNase$enhancer %in% DAE_ATAC$enhancer|
                                      DAE_DNase$enhancer %in% DAE_K27ac$enhancer|
                                      DAE_DNase$enhancer %in% DAE_K4me1$enhancer|
                                      DAE_DNase$enhancer %in% DAE_K4me2$enhancer|
                                      DAE_DNase$enhancer %in% DAE_K4me3$enhancer|
                                      DAE_DNase$enhancer %in% DAE_K27me3$enhancer),])
colnames(DNase_DAE)[1] <- 'Enhancers'
dim(DNase_DAE)
dim(DAE_DNase)

K27ac_DAE <- data.frame(DAE_K27ac [(DAE_K27ac$enhancer %in% DAE_ATAC$enhancer|
                                      DAE_K27ac$enhancer %in% DAE_DNase$enhancer|
                                      DAE_K27ac$enhancer %in% DAE_K4me1$enhancer|
                                      DAE_K27ac$enhancer %in% DAE_K4me2$enhancer|
                                      DAE_K27ac$enhancer %in% DAE_K4me3$enhancer|
                                      DAE_K27ac$enhancer %in% DAE_K27me3$enhancer),])
colnames(K27ac_DAE) <- 'Enhancers'
dim(K27ac_DAE)
dim(DAE_K27ac)

K4me1_DAE <- data.frame(DAE_K4me1 [(DAE_K4me1$enhancer %in% DAE_ATAC$enhancer|
                                      DAE_K4me1$enhancer %in% DAE_DNase$enhancer|
                                      DAE_K4me1$enhancer %in% DAE_K27ac$enhancer|
                                      DAE_K4me1$enhancer %in% DAE_K4me2$enhancer|
                                      DAE_K4me1$enhancer %in% DAE_K4me3$enhancer|
                                      DAE_K4me1$enhancer %in% DAE_K27me3$enhancer),])
colnames(K4me1_DAE) <- 'Enhancers'
dim(K4me1_DAE)
dim(DAE_K4me1)

K4me2_DAE <- data.frame(DAE_K4me2 [(DAE_K4me2$enhancer %in% DAE_ATAC$enhancer|
                                      DAE_K4me2$enhancer %in% DAE_DNase$enhancer|
                                      DAE_K4me2$enhancer %in% DAE_K27ac$enhancer|
                                      DAE_K4me2$enhancer %in% DAE_K4me1$enhancer|
                                      DAE_K4me2$enhancer %in% DAE_K4me3$enhancer|
                                      DAE_K4me2$enhancer %in% DAE_K27me3$enhancer),])
colnames(K4me2_DAE)<-'Enhancers'
dim(K4me2_DAE)
dim(DAE_K4me2)

K4me3_DAE <- data.frame(DAE_K4me3 [(DAE_K4me3$enhancer %in% DAE_ATAC$enhancer|
                                      DAE_K4me3$enhancer %in% DAE_DNase$enhancer|
                                      DAE_K4me3$enhancer %in% DAE_K27ac$enhancer|
                                      DAE_K4me3$enhancer %in% DAE_K4me1$enhancer|
                                      DAE_K4me3$enhancer %in% DAE_K4me2$enhancer|
                                      DAE_K4me3$enhancer %in% DAE_K27me3$enhancer),])
colnames(K4me3_DAE) <- 'Enhancers'
dim(K4me3_DAE)
dim(DAE_K4me3)

K27me3_DAE <- data.frame(DAE_K27me3 [(DAE_K27me3$enhancer %in% DAE_ATAC$enhancer|
                                        DAE_K27me3$enhancer %in% DAE_DNase$enhancer|
                                        DAE_K27me3$enhancer %in% DAE_K27ac$enhancer|
                                        DAE_K27me3$enhancer %in% DAE_K4me1$enhancer|
                                        DAE_K27me3$enhancer %in% DAE_K4me2$enhancer|
                                        DAE_K27me3$enhancer %in% DAE_K4me3$enhancer),])
colnames(K27me3_DAE) <- 'Enhancers'
dim(K27me3_DAE)
dim(DAE_K27me3)

## Merging defined DAEs
DAE_merge <- unique(rbind(ATAC_DAE[1], DNase_DAE[1], K27ac_DAE[1], K4me1_DAE[1], 
                          K4me2_DAE[1], K4me3_DAE[1], K27me3_DAE[1]))
dim(DAE_merge)

## nDAE regions
nDAE <- pCR_raw[! row.names(pCR_raw) %in% DAE_merge$Enhancers, ]
dim(nDAE)
DAE_merge <- DAE_merge %>% separate(Enhancers, c('chr','start','end'))

#write.table(nDAE[1:3], 'nDAE_162454_coordinate.bed', quote=F, sep="\t", row.names=F, col.names=T)
#write.table(DAE,'DAE_39709_coordinate.bed', quote=F, sep="\t", row.names=F, col.names=T)
