
######
###### Link DAE and nDAE to target genes
######

## Enhancers linked to the genes based on overlap between coordinates using bedtools intersect
## Data including enhancer coordinates,Gene type, GeneID, gene name for different gene prediction tools.
E_G_CP <- read.csv("Enhancer_CP.Genes.csv")
E_G_GZ <- read.csv("Enhancer_GZ.Genes.csv")

## selecting protein coding and LincRNA genes for downstream analysis
E_G_CP_pr.linc <- unique(E_G_CP[(E_G_CP$Gene.Type == 'protein_coding'|E_G_CP$Gene.Type =='lincRNA'),])
dim(unique(E_G_CP_pr.linc[c(1:3)])) # number of enhancers (DAEs/nDAEs) linked to the CP-protein coding and -LincRNA genes
dim(unique(E_G_CP_pr.linc[5]))      # number of CP-protein coding and -LincRNA genes linked to enhancers (DAEs/nDAEs) 

E_G_GZ_pr.linc <- unique(E_G_GZ[(E_G_GZ$Gene.Type == 'protein_coding'|E_G_GZ$Gene.Type =='lincRNA'),])
dim(unique(E_G_GZ_pr.linc[c(1:3)])) # number of enhancers (DAEs/nDAEs) linked to the GZ-protein coding and -LincRNA genes
dim(unique(E_G_GZ_pr.linc[5]))      # number of GZ-protein coding and -LincRNA genes linked to enhancers (DAEs/nDAEs)

## Separating DAE and nDAE to check how many enhancers linked to genes
### DAE regions linked to genes
DAE_G_CP_pr.linc <- unique(E_G_CP_pr.linc[(E_G_CP_pr.linc$Enhancer.type == 'DAE'),])
dim(unique(DAE_G_CP_pr.linc[1:3])) # Number of DAEs linked to the CP genes
dim(unique(DAE_G_CP_pr.linc[5]))   # Number of CP Genes linked to the DAEs
DAE_G_GZ_pr.linc <- unique(E_G_GZ_pr.linc[(E_G_GZ_pr.linc$Enhancer.type == 'DAE'),])
dim(unique(DAE_G_GZ_pr.linc[1:3])) # Number of DAEs linked to the GZ genes
dim(unique(DAE_G_GZ_pr.linc[5]))   # Number of GZ Genes linked to the DAEs

### nDAE regions linked to genes
nDAE_G_CP_pr.linc <- unique(E_G_CP_pr.linc[(E_G_CP_pr.linc$Enhancer.type == 'nDAE'),])
dim(unique(nDAE_G_CP_pr.linc[1:3]))# Number of nDAEs linked to the CP genes
dim(unique(nDAE_G_CP_pr.linc[5]))  # Number of CP Genes linked to the nDAEs
nDAE_G_GZ_pr.linc <- unique(E_G_GZ_pr.linc[(E_G_GZ_pr.linc$Enhancer.type == 'nDAE'),])
dim(unique(nDAE_G_GZ_pr.linc[1:3]))# Number of nDAEs linked to the GZ genes
dim(unique(nDAE_G_GZ_pr.linc[5]))  # Number of GZ Genes linked to the nDAEs

#### Enhancers interacting with gene Type
DAE_G_CP_pr<-DAE_G_CP_pr.linc[DAE_G_CP_pr.linc$Gene.Type =='protein_coding',]# separate DAE linked to protein coding CP-genes
dim(unique(DAE_G_CP_pr[1:3])) # Number of DAE that linked to CP protein coding genes
DAE.CP_linc<-DAE_G_CP_pr.linc[DAE_G_CP_pr.linc$Gene.Type =='lincRNA',]# separate DAE linked to lincRNA CP-genes
dim(unique(DAE.CP_linc[1:3]))# Number of DAE that linked to CP lincRNA genes
overlap_DAE.CP_pr.linc <- DAE_G_CP_pr[DAE_G_CP_pr$chr %in% DAE.CP_linc$chr & 
                                        DAE_G_CP_pr$start %in% DAE.CP_linc$start &
                                        DAE_G_CP_pr$end %in% DAE.CP_linc$end,]
dim(unique(overlap_DAE.CP_pr.linc[1:3]))# Number of shared DAE linked to both protein coding and lincRNA CP-genes

DAE_G_GZ_pr<-DAE_G_GZ_pr.linc[DAE_G_GZ_pr.linc$Gene.Type =='protein_coding',]
dim(unique(DAE_G_GZ_pr[1:3]))
DAE.GZ_linc<-DAE_G_GZ_pr.linc[DAE_G_GZ_pr.linc$Gene.Type =='lincRNA',]
dim(unique(DAE.GZ_linc[1:3]))
overlap_DAE.GZ_pr.linc <- DAE_G_GZ_pr[DAE_G_GZ_pr$chr %in% DAE.GZ_linc$chr & 
                                        DAE_G_GZ_pr$start %in% DAE.GZ_linc$start &
                                        DAE_G_GZ_pr$end %in% DAE.GZ_linc$end,]
dim(unique(overlap_DAE.GZ_pr.linc[1:3]))


nDAE_G_CP_pr<-nDAE_G_CP_pr.linc[nDAE_G_CP_pr.linc$Gene.Type =='protein_coding',]# separate nDAE linked to protein coding CP-genes
dim(unique(nDAE_G_CP_pr[1:3]))# Number of nDAE that linked to CP protein coding genes
nDAE.CP_linc<-nDAE_G_CP_pr.linc[nDAE_G_CP_pr.linc$Gene.Type =='lincRNA',]# separate nDAE linked to lincRNA CP-genes
dim(unique(nDAE.CP_linc[1:3]))# Number of nDAE that linked to CP lincRNA genes
overlap_nDAE.CP_pr.linc <- nDAE_G_CP_pr[nDAE_G_CP_pr$chr %in% nDAE.CP_linc$chr & 
                                          nDAE_G_CP_pr$start %in% nDAE.CP_linc$start &
                                          nDAE_G_CP_pr$end %in% nDAE.CP_linc$end,]
dim(unique(overlap_nDAE.CP_pr.linc[1:3]))# Number of shared nDAE linked to both protein coding and lincRNA CP-genes

nDAE_G_GZ_pr<-nDAE_G_GZ_pr.linc[nDAE_G_GZ_pr.linc$Gene.Type =='protein_coding',]
dim(unique(nDAE_G_GZ_pr[1:3]))
nDAE.GZ_linc<-nDAE_G_GZ_pr.linc[nDAE_G_GZ_pr.linc$Gene.Type =='lincRNA',]
dim(unique(nDAE.GZ_linc[1:3]))
overlap_nDAE.GZ_pr.linc <- nDAE_G_GZ_pr[nDAE_G_GZ_pr$chr %in% nDAE.GZ_linc$chr & 
                                          nDAE_G_GZ_pr$start %in% nDAE.GZ_linc$start &
                                          nDAE_G_GZ_pr$end %in% nDAE.GZ_linc$end,]
dim(unique(overlap_nDAE.GZ_pr.linc[1:3]))

##### Enhancer (DAEs and nDAEs) overlap between CP and GZ regions
DAE_CP_enhancer <- unique(DAE_G_CP_pr.linc[1:3])
DAE_GZ_enhancer <- unique(DAE_G_GZ_pr.linc[1:3])
table(DAE_CP_enhancer$chr %in% DAE_GZ_enhancer$chr & 
        DAE_CP_enhancer$start %in% DAE_GZ_enhancer$start &
        DAE_CP_enhancer$end %in% DAE_GZ_enhancer$end)

nDAE_CP_enhancer <- unique(nDAE_G_CP_pr.linc[1:3])
nDAE_GZ_enhancer <- unique(nDAE_G_GZ_pr.linc[1:3])
table(nDAE_CP_enhancer$chr %in% nDAE_GZ_enhancer$chr & 
        nDAE_CP_enhancer$start %in% nDAE_GZ_enhancer$start &
        nDAE_CP_enhancer$end %in% nDAE_GZ_enhancer$end)

##### CP genes interacting with enhancers (DAEs and nDAEs)
DAE_G.CP_genes_unique <- unique(DAE_G_CP_pr.linc[5])
nDAE_G.CP_genes_unique <- unique(nDAE_G_CP_pr.linc[5])
table(DAE_G.CP_genes_unique$GeneID %in% nDAE_G.CP_genes_unique$GeneID)

DAE_G.GZ_genes_unique <- unique(DAE_G_GZ_pr.linc[5])
nDAE_G.GZ_genes_unique <- unique(nDAE_G_GZ_pr.linc[5])
table(DAE_G.GZ_genes_unique$GeneID %in% nDAE_G.GZ_genes_unique$GeneID)

#### DAE and nDAE specific genes
specific_DAE_G.CP <- DAE_G_CP_pr.linc[!(DAE_G_CP_pr.linc$GeneID %in% nDAE_G_CP_pr.linc$GeneID), ]
dim(unique(specific_DAE_G.CP[5]))
specific_DAE_G.GZ <- DAE_G_GZ_pr.linc[!(DAE_G_GZ_pr.linc$GeneID %in% nDAE_G_GZ_pr.linc$GeneID), ]
dim(unique(specific_DAE_G.GZ[5]))

specific_nDAE_G.CP <- nDAE_G_CP_pr.linc[!(nDAE_G_CP_pr.linc$GeneID %in% DAE_G_CP_pr.linc$GeneID), ]
dim(unique(specific_nDAE_G.CP[5]))
specific_nDAE_G.GZ <- nDAE_G_GZ_pr.linc[!(nDAE_G_GZ_pr.linc$GeneID %in% DAE_G_GZ_pr.linc$GeneID), ]
dim(unique(specific_nDAE_G.GZ[5]))

##### Gene (DAEs and nDAEs) overlap between CP and GZ regions
DAE_CP_genes <- unique(DAE_G_CP_pr.linc[5])
DAE_GZ_genes <- unique(DAE_G_GZ_pr.linc[5])
table(DAE_CP_genes$GeneID %in% DAE_GZ_genes$GeneID)

nDAE_CP_genes <- unique(specific_nDAE_G.CP[5])
nDAE_GZ_genes <- unique(specific_nDAE_G.GZ[5])
table(nDAE_CP_genes$GeneID %in% nDAE_GZ_genes$GeneID)

## Making a table of DAEs and nDAEs target genes for downstream analysis
DAE.nDAE_CP_genes <- unique(rbind(DAE_G_CP_pr.linc, specific_nDAE_G.CP))
dim(unique(DAE.nDAE_CP_genes[5]))

DAE.nDAE_GZ_genes <- unique(rbind(DAE_G_GZ_pr.linc, specific_nDAE_G.GZ))
dim(unique(DAE.nDAE_GZ_genes[5]))

#write.table(DAE.nDAE_CP_genes), "Enhancer_Genes.CP.bed", quote=F, sep="\t", row.names=F, col.names=T)
#write.table(DAE.nDAE_GZ_genes), "Enhancer_Genes.GZ.bed", quote=F, sep="\t", row.names=F, col.names=T)

######
###### Expression
######

## Read DAEs and nDAEs linked to CP/GZ genes
## Data including enhancer coordinates,geneID, gene type, gene name (HiC) and enhacer type.
DAE.nDAE_genes <- read.delim("Enhancer_Genes.CP.bed")
#DAE.nDAE_genes <- read.delim("Enhancer_Genes.GZ.bed")

## Read brain related Expression data including geneID, gene name and expression level from ENCODE,Roadmap and AllenBrain
Expression <- read.csv('Expression_data.csv')

## Compare gene expression across different barin parts
## Extracting expression level for initial table "DAE.nDAE_genes"
DAE.nDAE_genes <- DAE.nDAE_genes %>% 
  mutate (`Frontal cortex` = Expression$Fetal.Frontal.Cortex [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Temporal lobe` = Expression$Fetal.Temporal.lobe [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Parietal lobe` = Expression$Fetal.Parietal.lobe [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Diencephalon lobe` = Expression$Fetal.Diencephalon.lobe [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Occipital lobe` = Expression$Fetal.Occipital.lobe [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Cerebellum` = Expression$Fetal.Cerebellum [match (GeneID, Expression$GeneID)]) %>%
  mutate (`Spinal Cord` = Expression$Fetal.Spinal.Cord [match (GeneID, Expression$GeneID)])

## separate expression level columns which added to "DAE.nDAE_genes" table for drawing plot
Expr_Brain.Tissues <- unique(DAE.nDAE_genes[c(4, 7:14)])
## Remove low regions and calculate log data
Expr_Brain.Tissues_melt <- melt(Expr_Brain.Tissues, id=c("GeneID","Enhancer.Type"), na.rm = T) %>%
  filter (value > 0) %>%
  mutate (log2_value = log2(value +1))

ggplot(Expr_Brain.Tissues_melt, aes(x=variable, y=log2_value, fill=Enhancer.Type, color=Enhancer.Type)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, outlier.size=0.6, alpha=0.8,notch = F) +
  scale_y_continuous("Log2(FPKM)") +
  theme_classic() +
  coord_cartesian(ylim=c(0, 9))+
  theme(axis.text.x=element_text(size=9, color="black", angle=90, hjust=1, vjust=0.3),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        legend.text=element_text(color="black", size=9), legend.position = "none") +
  scale_fill_manual(name="Putative Enhancer", values=c("gray", "black")) +
  scale_color_manual(name="Putative Enhancer", values = c("red", "red")) 

## Comparing expression level across Fetal brain time points and compare with adult brain
## To extract expression level the GeneID of initial table "DAE.nDAE_genes" matched with expression table
DAE.nDAE_genes <- DAE.nDAE_genes %>%
  mutate (`Fetal sources(Mean)` = Expression$Fetal.Brain_Mean [match (GeneID, Expression$GeneID)]) %>%
  mutate (`12PCW` = Expression$Fetal.Brain_Yan [match (GeneID, Expression$GeneID)]) %>%
  mutate (`15-17PCW` = Expression$Fetal.Brain_Ubieta [match (GeneID, Expression$GeneID)]) %>%
  mutate (`17PCW` = Expression$Fetal.Brain_Roadmap [match (GeneID, Expression$GeneID)]) %>%
  mutate (`81Y` = Expression$Adult.Brain_Roadmap [match (GeneID, Expression$GeneID)]) 

## Separating expression level columns which added to "DAE.nDAE_genes" table for drawing plot
Expr_Brain.TimePoints <- unique(DAE.nDAE_genes[c(4,7,15:19)])
## Separating regions that linked to Fetal and Adult samples to give them names as "Fetal Brain" and "Adult Brain" in order to group them in plot
Expr_Brain.TimePoints_melt_FB <- melt(Expr_Brain.TimePoints[c(1,2,4,5,6,3)], id=c("GeneID", "Enhancer.Type"), na.rm = T) %>%
  mutate (group = 'Fetal Brain')
Expr_Brain.TimePoints_melt_AB <- melt(Expr_Brain.TimePoints[c(1,2,7)], id=c("GeneID", "Enhancer.Type"), na.rm = T) %>%
  mutate (group = 'Adult Brain') %>%
  filter (value > 0)
Expr_Brain.TimePoints_melt <- na.omit(rbind(Expr_Brain.TimePoints_melt_FB, Expr_Brain.TimePoints_melt_AB)) %>%
  mutate (log2_value = log2(value+1))
## order x-axis
Expr_Brain.TimePoints_melt$variable <- factor(Expr_Brain.TimePoints_melt$variable,
                                              levels = c("12PCW", "15-17PCW", "17PCW", "81Y","Fetal sources(Mean)"))

ggplot(Expr_Brain.TimePoints_melt, aes(x=variable, y=log2_value, fill=Enhancer.Type, color=group)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("Log2(Normlized count)") +
  theme_classic() +
  coord_cartesian(ylim=c(0, 11)) +
  theme(axis.text.x=element_text(size=10, color='black', angle=90, vjust=0.3, hjust=1),
        axis.title=element_blank(), axis.text.y=element_text(size=10, color='black'),
        legend.position="none") +
  scale_fill_manual(name="Putative Enhancer", values=c("gray", "black")) +
  scale_color_manual(name="Expression Source", values=c("blue", "red"))

######
###### Number of Gene per Enhancer
######

# Read enhancer_gene tables and separate DAEA and nDAE interactions
## Read E_G.CP table

Enhancer_Genes.CP <- read.csv("C:/Users/g/Desktop/DAE.CP_genes.csv")
## Separating DAE-CP regions liked to genes
DAE_CP_genes <- Enhancer_Genes.CP[Enhancer_Genes.CP$variability == "DAE",(1:4)] %>% unique() 
## Read E_G.GZ table
Enhancer_Genes.GZ <- read.csv("C:/Users/g/Desktop/DAE.GZ_genes.csv")
## Separating DAE-GZ regions liked to genes
DAE_GZ_genes <- Enhancer_Genes.GZ[Enhancer_Genes.GZ$variability == "DAE",(1:4)] %>% unique()

######
###### Number of Enhancer per Gene
######
## For nDAE this part should be run for "specific_nDAE_G.GZ/CP"
#### Counting number of enhancer per gene for CP and GZ
EpG_CP <- data.frame(DAE_CP_genes %>% group_by(across(4)) %>%
                       summarize(EpG=n()) %>%
                       mutate(HiC="CP"))
EpG_GZ <- data.frame(DAE_GZ_genes %>% group_by(across(4)) %>%
                       summarize(EpG=n()) %>%
                       mutate(HiC="CP"))
## grouping enhancers with >5 genes together and count gene per enhancer for CP and GZ
EpG_CP$EpG[EpG_CP$EpG > 5] <- 5 
Table.EpG_CP <- data.frame(table(EpG_CP$EpG)) %>%
  mutate(type="CP")
EpG_GZ$EpG[EpG_GZ$EpG > 5] <- 5 
Table.EpG_GZ <- data.frame(table(EpG_GZ$EpG)) %>%
  mutate(type="GZ")
EpG_tables <- rbind(Table.EpG_CP, Table.EpG_GZ)

ggplot(EpG_tables, aes(x=Var1, y=Freq, group=type)) +
  geom_line(color="black", size=0.5) +
  geom_point(aes(color=type), size=1.5) + theme_classic() +
  theme(axis.text.x=element_text(size=10, color='black', angle=0, vjust=0.3, hjust=1),
        axis.title=element_blank(), axis.text.y=element_text(size=10, color='black')) +
  scale_color_manual(values=c("red", "blue"))

## Merging expression data with EpG tables to see how expression level is changing when a gene interaction with multiple enhancers
Expr_EpG <- merge(Expression[c(1,15)], EpG_CP[-3], 'GeneID') # Also the same for "EpG_GZ[-3]"
Expr_EpG_melt <- melt(Expr_EpG, id=c('GeneID',"EpG"), na.rm = T) %>%
  mutate(log2_value=log2(value + 1))
Expr_EpG_melt$EpG[Expr_EpG_melt$EpG > 5] <- 5

ggplot(Expr_EpG_melt, aes(x=EpG, y=log2_value, fill=as.factor(EpG))) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("Log2(FPKM)") +
  scale_x_continuous("DAE per gene (GZ)") +
  theme_classic() +
  coord_cartesian(ylim=c(0, 11)) +
  theme(axis.text.x=element_text(size=9, color="black", angle=0, hjust=0.3, vjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("gray90", "gray75", "gray60", "gray45", "gray30"))

#### Counting number of Gene per Enhancer
## Counting number of gene per enhancer for CP and GZ
GpE_CP <- data.frame(DAE_CP_genes %>% group_by(across(1:3)) %>%
                       summarize(GpE=n()) %>%
                       mutate(HiC="CP"))
GpE_GZ <- data.frame(DAE_GZ_genes %>% group_by(across(1:3)) %>%
                       summarize(GpE=n()) %>%
                       mutate(HiC="GZ"))
## Grouping enhancers with >5 genes together and count gene per enhancer for CP and GZ
GpE_CP$GpE[GpE_CP$GpE > 5] <- 5
Table.GpE_CP <- data.frame(table(GpE_CP$GpE)) %>%
  mutate(type="CP")
GpE_GZ$GpE[GpE_GZ$GpE > 5] <- 5
Table.GpE_GZ <- data.frame(table(GpE_GZ$GpE)) %>%
  mutate(type="GZ")
GpE_tables <- rbind(Table.GpE_CP, Table.GpE_GZ)

ggplot(GpE_tables, aes(x=Var1, y=Freq, group=type)) +
  geom_line(color="black", size=1, position = position_jitter(width=0.15, height=0.1)) +
  geom_point(aes(color=type), size=2, position = position_jitter(width=0.2, height=0.1)) + theme_classic() +
  theme(axis.text.x=element_text(size=10, color='black', angle=0, vjust=0.3, hjust=1),
        axis.title=element_blank(), axis.text.y=element_text(size=10, color='black')) +
  scale_color_manual(values=c("red", "blue"))

## Gene per enhancer enrichment(nCER,GC and conservation scores)
GpE_CP.GZ<-rbind(GpE_CP, GpE_GZ)# combine CP and GZ related data
GpE_CP.GZ_coor <- GpE_CP.GZ[1:3] %>% # select coordinates
  mutate (coordinate=paste0(chr, ":", start, "-", end))

### ncER
####nCER table is generated by intersecting all enhancers coordinates and nCER coordinates using bedtools
nCER_pCR <- read.delim('pCR_nCER.bed', header=F)
## separating related enhancer.CP and GZ coordinates from nCER table
GpE_CP.GZ_nCER <- nCER_pCR[(nCER_pCR$V1 %in% GpE_CP.GZ_coor$chr & nCER_pCR$V2 %in% GpE_CP.GZ_coor$start
                            & nCER_pCR$V3 %in% GpE_CP.GZ_coor$end), ] # separate related coordinates from nCER table
colnames(GpE_CP.GZ_nCER) <- c('chr', 'start', 'end', 'nCER')
GpE_CP.GZ_nCER <- merge(GpE_CP.GZ, GpE_CP.GZ_nCER, by=c('chr', 'start', 'end'))
head(GpE_CP.GZ_nCER)

### GC content
## Extracting GC content for enhancer.CP and GZ coordinates 
GpE_CP.GZ_GC <- data.frame(GCcontent(Hsapiens, GRanges(GpE_CP.GZ_coor$coordinate), as.prob=T))
colnames(GpE_CP.GZ_GC) <- "GC"
GpE_CP.GZ_GC <- cbind(GpE_CP.GZ, GpE_CP.GZ_GC)

### Conservation Score
## Extracting Conservation Score for enhancer.CP and GZ coordinates 
GpE_CP.GZ_phastcons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(GpE_CP.GZ_coor$coordinate)))
colnames(GpE_CP.GZ_phastcons)[6] <- "phastcons"
GpE_CP.GZ_phastcons <- cbind(GpE_CP.GZ, GpE_CP.GZ_phastcons[6])

#### Grouping and plot
## Group enhancers with >5 genes together
GpE_CP.GZ_nCER$GpE[GpE_CP.GZ_nCER$GpE > 5] <- 5
GpE_CP.GZ_GC$GpE[GpE_CP.GZ_GC$GpE > 5] <- 5
GpE_CP.GZ_phastcons$GpE[GpE_CP.GZ_phastcons$GpE > 5] <- 5

#### Drawing plot for each score we need to use the related score as input
ggplot(GpE_CP.GZ_nCER, aes(x=as.factor(GpE), y=nCER, fill=HiC)) +
  geom_boxplot(lwd=0.34, outlier.colour="black", outlier.shape=16, outlier.size=0.8, alpha=0.8) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.text.x=element_text(size=10, color="black", angle=0, hjust=0.3, vjust=1),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_manual(values=c("#FF0000", "#0000FF"))

######
###### DNA sequence features
######
## Comparing different subsets of enhancer.CP and GZ for nCER,GC and conservation scores
## Read related data
DAE <- read.delim('DAE_39709.bed')## read all DAE coordinates
nDAE <- read.delim('nDAE_162454.bed')## read all nDAE coordinates
## Separating DAE and nDAE regions interacting to genes (CP or GZ)
E_G.CP <- read.csv('Enhancer_CP.Genes.csv') ## read enhancers linked to the CP-genes
DAE_G.CP <- E_G.CP %>%
  filter(Enhancer.type == "DAE") %>%
  dplyr::select("chr", "start", "end", "GeneID") %>% unique()
nDAE_G.CP <- E_G.CP %>%
  filter(Enhancer.type == "nDAE") %>%
  dplyr::select("chr", "start", "end", "GeneID") %>% unique()

E_G.GZ <- read.csv('Enhancer_GZ.Genes.csv')## read enhancers linked to the GZ-genes
DAE_G.GZ <- E_G.GZ %>%
  filter(Enhancer.type == "DAE") %>%
  dplyr::select("chr", "start", "end", "GeneID") %>% unique()
nDAE_G.GZ <- E_G.GZ %>%
  filter(Enhancer.type == "nDAE") %>%
  dplyr::select("chr", "start", "end", "GeneID") %>% unique()

## Linking the separated CP and GZ data with OMIM to include those regions that linked to OMIM genes as a group in plot
OMIM <- read.delim('OMIM_data.bed')
DAE_G.CP_OMIM <- merge(DAE_G.CP, OMIM, 'GeneID')
DAE_G.CP_OMIM <- DAE_G.CP_OMIM[!DAE_G.CP_OMIM$MIM.description == '', 2:4] %>% unique()
DAE_G.GZ_OMIM <- merge(DAE_G.GZ, OMIM, 'GeneID')
DAE_G.GZ_OMIM <- DAE_G.GZ_OMIM[!DAE_G.GZ_OMIM$MIM.description == '', 2:4] %>% unique()
nDAE_G.CP_OMIM <- merge(nDAE_G.CP, OMIM, 'GeneID')
nDAE_G.CP_OMIM <- nDAE_G.CP_OMIM[!nDAE_G.CP_OMIM$MIM.description == '', 2:4] %>% unique()
nDAE_G.GZ_OMIM <- merge(nDAE_G.GZ, OMIM, 'GeneID')
nDAE_G.GZ_OMIM <- nDAE_G.GZ_OMIM[!nDAE_G.GZ_OMIM$MIM.description == '', 2:4] %>% unique()

### Making coordinate column
DAE <- DAE %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
nDAE <- nDAE %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
DAE_G.CP <- DAE_G.CP %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
nDAE_G.CP <- nDAE_G.CP %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
DAE_G.GZ <- DAE_G.GZ %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
nDAE_G.GZ <- nDAE_G.GZ %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
DAE_G.CP_OMIM <- DAE_G.CP_OMIM %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
nDAE_G.CP_OMIM <- nDAE_G.CP_OMIM %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
DAE_G.GZ_OMIM <- DAE_G.GZ_OMIM %>% mutate (coordinate=paste0(chr, ":", start, "-", end))
nDAE_G.GZ_OMIM <- nDAE_G.GZ_OMIM %>% mutate (coordinate=paste0(chr, ":", start, "-", end))

### Matching each enhancer data with features (nCER,GC content and conservation score)
### ncER
nCER_pCR <- read.delim('pCR_nCER.bed', header=F)
## Separating nCER score for DAEs
DAE_ncER <- nCER_pCR[(nCER_pCR$V1 %in% DAE$chr & nCER_pCR$V2 %in% DAE$start
                      & nCER_pCR$V3 %in% DAE$end),]
## Add name to the columns
colnames(DAE_ncER)[4] <- c('DAE')
DAE_ncER_melt <- reshape2::melt(DAE_ncER[4]) %>% 
  mutate(group="All Regions") %>%
  mutate(name="All Regions")
## Separating nCER score for nDAEs
nDAE_ncER <- nCER_pCR[(nCER_pCR$V1 %in% nDAE$chr & nCER_pCR$V2 %in% nDAE$start
                       & nCER_pCR$V3 %in% nDAE$end),]
colnames(nDAE_ncER)[4]<-c('nDAE')
nDAE_ncER_melt <- reshape2::melt(nDAE_ncER[4]) %>%
  mutate(group="All Regions") %>%
  mutate(name="All Regions")
## Separating nCER score for DAEs are linked to CP genes
DAE_G.CP_nCER <- nCER_pCR[(nCER_pCR$V1 %in% DAE_G.CP$chr & nCER_pCR$V2 %in% DAE_G.CP$start
                           & nCER_pCR$V3 %in% DAE_G.CP$end),]
colnames(DAE_G.CP_nCER)[4] <- c('DAE')
DAE_G.CP_nCER_melt <- reshape2::melt(DAE_G.CP_nCER[4]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")
## Separating nCER score for nDAEs are linked to CP genes
nDAE_G.CP_nCER <- nCER_pCR[(nCER_pCR$V1 %in% nDAE_G.CP$chr & nCER_pCR$V2 %in% nDAE_G.CP$start
                            & nCER_pCR$V3 %in% nDAE_G.CP$end),]
colnames(nDAE_G.CP_nCER)[4] <- c('nDAE')
nDAE_G.CP_nCER_melt <- reshape2::melt(nDAE_G.CP_nCER[4]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")
## Separating nCER score for DAEs are linked to GZ genes
DAE_G.GZ_nCER <- nCER_pCR[(nCER_pCR$V1 %in% DAE_G.GZ$chr & nCER_pCR$V2 %in% DAE_G.GZ$start
                           & nCER_pCR$V3 %in% DAE_G.GZ$end),]
colnames(DAE_G.GZ_nCER)[4] <- c('DAE')
DAE_G.GZ_nCER_melt <- reshape2::melt(DAE_G.GZ_nCER[4]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_GZ")
## Separating nCER score for nDAEs are linked to GZ genes
nDAE_G.GZ_nCER <- nCER_pCR[(nCER_pCR$V1 %in% nDAE_G.GZ$chr & nCER_pCR$V2 %in% nDAE_G.GZ$start
                            & nCER_pCR$V3 %in% nDAE_G.GZ$end),]
colnames(nDAE_G.GZ_nCER)[4] <- c('nDAE')
nDAE_G.GZ_nCER_melt <- reshape2::melt(nDAE_G.GZ_nCER[4]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_GZ")
## Separating nCER score for DAEs are linked to CP.OMIM genes
DAE_G.CP_OMIM_nCER <- nCER_pCR[(nCER_pCR$V1 %in% DAE_G.CP_OMIM$chr & nCER_pCR$V2 %in% DAE_G.CP_OMIM$start
                                & nCER_pCR$V3 %in% DAE_G.CP_OMIM$end),]
colnames(DAE_G.CP_OMIM_nCER)[4] <- c('DAE')
DAE_G.CP_OMIM_nCER_melt <- reshape2::melt(DAE_G.CP_OMIM_nCER[4]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")
## Separating nCER score for nDAEs are linked to CP.OMIM genes
nDAE_G.CP_OMIM_nCER <- nCER_pCR[(nCER_pCR$V1 %in% nDAE_G.CP_OMIM$chr & nCER_pCR$V2 %in% nDAE_G.CP_OMIM$start
                                 & nCER_pCR$V3 %in% nDAE_G.CP_OMIM$end),]
colnames(nDAE_G.CP_OMIM_nCER)[4] <- c('nDAE')
nDAE_G.CP_OMIM_nCER_melt <- reshape2::melt(nDAE_G.CP_OMIM_nCER[4]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")
## Separating nCER score for DAEs are linked to GZ.OMIM genes
DAE_G.GZ_OMIM_nCER <- nCER_pCR[(nCER_pCR$V1 %in% DAE_G.GZ_OMIM$chr & nCER_pCR$V2 %in% DAE_G.GZ_OMIM$start
                                & nCER_pCR$V3 %in% DAE_G.GZ_OMIM$end),]
colnames(DAE_G.GZ_OMIM_nCER)[4] <- c('DAE')
DAE_G.GZ_OMIM_nCER_melt <- reshape2::melt(DAE_G.GZ_OMIM_nCER[4]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_GZ.OMIM")
## Separating nCER score for nDAEs are linked to GZ.OMIM genes
nDAE_G.GZ_OMIM_nCER <- nCER_pCR[(nCER_pCR$V1 %in% nDAE_G.GZ_OMIM$chr & nCER_pCR$V2 %in% nDAE_G.GZ_OMIM$start
                                 & nCER_pCR$V3 %in% nDAE_G.GZ_OMIM$end),]
colnames(nDAE_G.GZ_OMIM_nCER)[4] <- c('nDAE')
nDAE_G.GZ_OMIM_nCER_melt <- reshape2::melt(nDAE_G.GZ_OMIM_nCER[4]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_GZ.OMIM")
## Combining all groups with nCER scores
nCER <- rbind(DAE_ncER_melt, nDAE_ncER_melt, 
              DAE_G.CP_nCER_melt, nDAE_G.CP_nCER_melt,
              DAE_G.GZ_nCER_melt, nDAE_G.GZ_nCER_melt,
              DAE_G.CP_OMIM_nCER_melt, nDAE_G.CP_OMIM_nCER_melt,
              DAE_G.GZ_OMIM_nCER_melt, nDAE_G.GZ_OMIM_nCER_melt)

ggplot(nCER, aes(x=name, y=value, fill=variable, color=group)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("ncER Percentile") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(),legend.position="none") +
  scale_fill_manual(name="Putative Enhancer", values=c("#E6E6E6", "#1A1A1A")) +
  scale_color_manual(name="Enriched Type", values=c("black", "#CC0000",'#0000CC'))


### GC content
## Adding GC score to all DAEs/nDAEs
DAE_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(DAE$coordinate)), as.prob=T))
## Adding column name
colnames(DAE_GC) <- "DAE"
DAE_GC_melt <- reshape2::melt(DAE_GC) %>%
  mutate(group="All Regions") %>%
  mutate(name="All Regions")
nDAE_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(nDAE$coordinate)), as.prob=T))
colnames (nDAE_GC) <- "nDAE"
nDAE_GC_melt <- reshape2::melt(nDAE_GC) %>%
  mutate(group="All Regions") %>%
  mutate(name="All Regions")
## Adding GC score to DAEs/nDAEs are linked to CP genes
DAE_G.CP_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(DAE_G.CP$coordinate)), as.prob=T))
colnames(DAE_G.CP_GC) <- "DAE"
DAE_G.CP_GC_melt <-reshape2::melt(DAE_G.CP_GC) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")
nDAE_G.CP_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(nDAE_G.CP$coordinate)), as.prob=T))
colnames(nDAE_G.CP_GC) <- "nDAE"
nDAE_G.CP_GC_melt <- reshape2::melt(nDAE_G.CP_GC) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")
## Adding GC score to DAEs/nDAEs are linked to GZ genes
DAE_G.GZ_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(DAE_G.GZ$coordinate)), as.prob=T))
colnames(DAE_G.GZ_GC) <- "DAE"
DAE_G.GZ_GC_melt <- reshape2::melt(DAE_G.GZ_GC) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_GZ")
nDAE_G.GZ_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(nDAE_G.GZ$coordinate)), as.prob=T))
colnames(nDAE_G.GZ_GC) <- "nDAE"
nDAE_G.GZ_GC_melt <- reshape2::melt(nDAE_G.GZ_GC) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_GZ")

## Adding GC score to DAEs/nDAEs are linked to CP.OMIM genes
DAE_G.CP_OMIM_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(DAE_G.CP_OMIM$coordinate)), as.prob=T))
colnames(DAE_G.CP_OMIM_GC) <- "DAE"
DAE_G.CP_OMIM_GC_melt <- reshape2::melt(DAE_G.CP_OMIM_GC) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")
nDAE_G.CP_OMIM_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(nDAE_G.CP_OMIM$coordinate)), as.prob=T))
colnames(nDAE_G.CP_OMIM_GC) <- "nDAE"
nDAE_G.CP_OMIM_GC_melt <- reshape2::melt(nDAE_G.CP_OMIM_GC) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")
## Adding GC score to DAEs/nDAEs are linked to GZ.OMIM genes
DAE_G.GZ_OMIM_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(DAE_G.GZ_OMIM$coordinate)), as.prob=T))
colnames(DAE_G.GZ_OMIM_GC) <- "DAE"
DAE_G.GZ_OMIM_GC_melt <- reshape2::melt(DAE_G.GZ_OMIM_GC) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_GZ.OMIM")
nDAE_G.GZ_OMIM_GC <- data.frame(GCcontent(Hsapiens, GRanges(unique(nDAE_G.GZ_OMIM$coordinate)), as.prob=T))
colnames(nDAE_G.GZ_OMIM_GC) <- "nDAE"
nDAE_G.GZ_OMIM_GC_melt <- reshape2::melt(nDAE_G.GZ_OMIM_GC) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_GZ.OMIM")

## Combining all groups with GC scores
GC <- rbind(DAE_GC_melt, nDAE_GC_melt, 
            DAE_G.CP_GC_melt, nDAE_G.CP_GC_melt,
            DAE_G.GZ_GC_melt, nDAE_G.GZ_GC_melt,
            DAE_G.CP_OMIM_GC_melt, nDAE_G.CP_OMIM_GC_melt,
            DAE_G.GZ_OMIM_GC_melt, nDAE_G.GZ_OMIM_GC_melt)

ggplot(GC, aes(x=name, y=value,fill= variable, color = group)) +
  geom_boxplot(lwd=0.35 ,outlier.colour="black", outlier.shape=16, outlier.size=0.9, alpha=0.9)+
  scale_y_continuous("GC Score")+
  theme_classic() +
  coord_cartesian(ylim = c(0, 1))+
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size=9,color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text())+
  scale_fill_manual(name= "Putative Enhancer", values = c("#E6E6E6", "#1A1A1A"))+
  scale_color_manual(name = "Enriched Type", values = c("black","#CC0000",'#0000CC'))

### Conservation Score
## Adding Conservation score to all DAEs/nDAEs
DAE_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(DAE$coordinate)))
colnames(DAE_phastCons)[6] <- "DAE"
DAE_phastCons_melt <- reshape2::melt(DAE_phastCons[6]) %>%
  mutate(group="All Regions") %>%
  mutate(name="All Regions")
nDAE_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(nDAE$coordinate)))
colnames(nDAE_phastCons)[6] <- "nDAE"
nDAE_phastCons_melt <- reshape2::melt(nDAE_phastCons[6]) %>%
  mutate(group = "All Regions") %>%
  mutate(name = "All Regions")

## Adding Conservation score to DAEs/nDAEs are linked to CP genes
DAE_G.CP_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(DAE_G.CP$coordinate))))
colnames(DAE_G.CP_phastCons)[6] <- "DAE"
DAE_G.CP_phastCons_melt <- reshape2::melt(DAE_G.CP_phastCons[6]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")
nDAE_G.CP_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(nDAE_G.CP$coordinate))))
colnames(nDAE_G.CP_phastCons)[6] <- "nDAE"
nDAE_G.CP_phastCons_melt <- reshape2::melt(nDAE_G.CP_phastCons[6]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_CP")

## Adding Conservation score to DAEs/nDAEs are linked to GZ genes
DAE_G.GZ_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(DAE_G.GZ$coordinate))))
colnames(DAE_G.GZ_phastCons)[6] <-"DAE"
DAE_G.GZ_phastCons_melt <- reshape2::melt(DAE_G.GZ_phastCons[6]) %>%
  mutate(group="Linked to genes") %>%
  mutate(name="HiC_GZ")
nDAE_G.GZ_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(nDAE_G.GZ$coordinate))))
colnames(nDAE_G.GZ_phastCons)[6] <- "nDAE"
nDAE_G.GZ_phastCons_melt <- reshape2::melt(nDAE_G.GZ_phastCons[6]) %>%
  mutate(group = "Linked to genes") %>%
  mutate(name = "HiC_GZ")

## Adding Conservation score to DAEs/nDAEs are linked to CP.OMIM genes
DAE_G.CP_OMIM_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(DAE_G.CP_OMIM$coordinate))))
colnames(DAE_G.CP_OMIM_phastCons)[6] <- "DAE"
DAE_G.CP_OMIM_phastCons_melt <- reshape2::melt(DAE_G.CP_OMIM_phastCons[6]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")
nDAE_G.CP_OMIM_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(nDAE_G.CP_OMIM$coordinate))))
colnames(nDAE_G.CP_OMIM_phastCons)[6] <- "nDAE"
nDAE_G.CP_OMIM_phastCons_melt <- reshape2::melt(nDAE_G.CP_OMIM_phastCons[6]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_CP.OMIM")

## Adding Conservation score to DAEs/nDAEs are linked to GZ.OMIM genes
DAE_G.GZ_OMIM_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(DAE_G.GZ_OMIM$coordinate))))
colnames(DAE_G.GZ_OMIM_phastCons)[6] <-"DAE"
DAE_G.GZ_OMIM_phastCons_melt <- reshape2::melt(DAE_G.GZ_OMIM_phastCons[6]) %>%
  mutate(group="Linked to OMIM genes") %>%
  mutate(name="HiC_GZ.OMIM")
nDAE_G.GZ_OMIM_phastCons <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(unique(nDAE_G.GZ_OMIM$coordinate))))
colnames(nDAE_G.GZ_OMIM_phastCons)[6] <- "nDAE"
nDAE_G.GZ_OMIM_phastCons_melt <- reshape2::melt(nDAE_G.GZ_OMIM_phastCons[6]) %>%
  mutate(group= "Linked to OMIM genes") %>%
  mutate(name= "HiC_GZ.OMIM")

## Combining all groups with Conservation scores
phastcons <- rbind(DAE_phastCons_melt, nDAE_phastCons_melt,
                   DAE_G.CP_phastCons_melt, nDAE_G.CP_phastCons_melt,
                   DAE_G.GZ_phastCons_melt, nDAE_G.GZ_phastCons_melt,
                   DAE_G.CP_OMIM_phastCons_melt, nDAE_G.CP_OMIM_phastCons_melt,
                   DAE_G.GZ_OMIM_phastCons_melt, nDAE_G.GZ_OMIM_phastCons_melt)

ggplot(phastcons, aes(x=name, y=value,fill=variable, color=group)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, outlier.size=0, alpha=0.8) +
  scale_y_continuous("Score") +
  theme_classic() +
  coord_cartesian(ylim=c(0, 0.6)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none") +
  scale_fill_manual(name="Putative Enhancer", values=c("#E6E6E6", "#1A1A1A")) +
  scale_color_manual(name="Enriched Type", values=c("black", "#CC0000", "#0000CC"))

### Add orion and CADD scores to enhacer groups
Orion_score <- read.csv('Orion.score.csv',header = T,sep = ',')
CADD_score <- read.csv('CADD.score.csv',header = T,sep = ',')

ggplot(CADD_score, aes(x=groups, y=CADD, fill= Enhancer_type, color=Enriched_type))+
  geom_boxplot(lwd=0.35 ,outlier.colour="black", outlier.shape=16, outlier.size=0.4, alpha=0.8) +
  theme_classic() +
  scale_fill_manual(values = c('#E6E6E6', '#1A1A1A')) + 
  scale_color_manual(values = c('black', '#CC0000','#0000CC')) +
  theme(axis.text.x = element_blank(), 
        axis.text.y =element_text(size=9,color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  ylim(0,15)

ggplot(Orion_score, aes(x=groups, y=Orion, fill= Enhancer_type, color=Enriched_type))+
  geom_boxplot(lwd=0.35 ,outlier.colour="black", outlier.shape=16, outlier.size=0.4, alpha=0.8) +
  theme_classic() +
  scale_fill_manual(values = c('#E6E6E6', '#1A1A1A')) + 
  scale_color_manual(values = c('black', '#CC0000','#0000CC')) +
  theme(axis.text.x = element_blank(), 
        axis.text.y =element_text(size=9,color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  ylim(-1,0.2)

######
###### PIL Score
######
## Read PLI (downloaded from Lek etal 2016) and enhancer_gene data
## Data including gene name and score
pLI <- read.csv('pLI.score.csv')

## Read DAEs and nDAEs linked to genes
E_G.CP <- read.delim('Enhancer_Genes.CP.bed')
E_G.GZ <- read.delim('Enhancer_Genes.GZ.bed')

## Merging PLI score with enhancer_gene data
E_G.CP_pLI <- merge(E_G.CP, pLI, by.x='GeneName', by.y='gene') %>%
  mutate(group='HiC_CP')
E_G.CP_pLI_melt <- reshape2::melt(E_G.CP_pLI[7:9]) %>% na.omit()

E_G.GZ_pLI<-merge(E_G.GZ,pLI,by.x='GeneName',by.y='gene') %>%
  mutate (group='HiC_GZ')
E_G.GZ_pLI_melt <- reshape2::melt(E_G.GZ_pLI[7:9]) %>% na.omit()

E_G.CP.GZ_pLI_melt <- rbind(E_G.CP_pLI_melt, E_G.GZ_pLI_melt) %>%
  mutate(Enhancer.Type=as.factor(Enhancer.Type)) %>%
  mutate(group=as.factor(group))

ggplot(E_G.CP.GZ_pLI_melt, aes(x=group, y=value, fill=Enhancer.Type, color=group)) +
  geom_boxplot(lwd=0.35 ,outlier.colour="black", outlier.shape=16, outlier.size=0.2, alpha=0.8) +
  scale_y_continuous("ncER Percentile") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle=90),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_manual(name="Putative Enhancer", values=c("gray", "black")) +
  scale_color_manual(name="Enriched Type", values=c("black", "black"))

######
###### LOF score
######
## Intersecting between each enhancer group and LOF regions are applied using bedtools intersect.
## Reading each enhancer_LoF data
## Data including enhancer coordinates, LoF coordinates, LoF type, LoF score.
DAE_LoF <- read.delim('DAE_39709_LoF.bed', header=F) %>%
  mutate(enhancer_type='All DAE')
DAE.CP_LoF <- read.delim('DAE_CP_LoF.bed', header=F) %>%
  mutate(enhancer_type='DAE linked to CP genesP')
DAE.GZ_LoF <- read.delim('DAE_GZ_LoF.bed', header=F) %>%
  mutate(enhancer_type='DAE linked to GZ genes')
nDAE_LoF <- read.delim('nDAE_162454_LoF.bed', header=F) %>%
  mutate(enhancer_type='All nDAE')
nDAE.CP_LoF <- read.delim('nDAE_CP_LoF.bed', header=F) %>%
  mutate (enhancer_type='nDAE linked to CP genes')
nDAE.GZ_LoF <- read.delim('nDAE_GZ_LoF.bed', header=F) %>%
  mutate (enhancer_type='nDAE linked to GZ genes')

enhancer_LoF <- rbind(DAE_LoF[7:9], DAE.CP_LoF[7:9], DAE.GZ_LoF[7:9],
                      nDAE_LoF[7:9], nDAE.CP_LoF[7:9], nDAE.GZ_LoF[7:9])

ggplot(enhancer_LoF, aes(x=V8, color=enhancer_type)) +
  geom_density(size=1) + theme_classic() +
  scale_x_continuous("LoF Score") +
  scale_color_manual(values=c('darkgray', 'black', 'red', 'green', 'orange', 'yellow')) +
  geom_vline(xintercept=0.5, colour="black", linetype="dashed") +
  theme(axis.text.x=element_text(size=9, color="black"),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_text(),
        axis.title.y=element_text())

######
###### TE enrichment
######
# Intersecting between each enhancer group and repeatmasker(downloaded from UCSC) are applied using bedtools intersect.
Repeatmasker <- unique(read.delim('RepeatMasker.bed', header=T))
## We can do it for other subsets as well (DAE_O.E.Ratio,DAE.CP_O.E.Ratio,DAE.GZ_O.E.Ratio)
repeat_enriched <- Repeatmasker[Repeatmasker$nDAE_O.E.Ratio >= 0.1,]# Threshold for enrichment
repeat_enriched <- unique(repeat_enriched[order(- repeat_enriched$nDAE_O.E.Ratio),])
repeat_enriched <- repeat_enriched[1:10,] ## select top 10 enriched TEs

ggplot(repeat_enriched, aes(y=reorder(TE_Name, nDAE_O.E.Ratio), x=nDAE_O.E.Ratio, nDAE_O.E.Ratio, fill=TE_class)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  scale_x_continuous(limits=c(0,12), breaks=seq(0, 12, by=1), expand=c(0, 0)) +
  theme(axis.text.x=element_text(size=10, color="black", angle=0, hjust=0.3, vjust=1.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_manual(values=c('Low_complexity'="darkblue", 'Simple_repeat'='darkgray',
                             'Unknown?'='darkorange', 'LTR'='darkgreen', 'DNA'='darkred',
                             'LINE'='black', 'Satellite'='purple'))

######
###### OMIM enrichment
######
## Each DAE_CP and -GZ are linked to OMIM phenotype using the same gene
DAE.CP_OMIM<-read.delim("DAE.CP_OMIM.phenotype.bed")
head(DAE.CP_OMIM)
DAE.CP_OMIM_top40 <- data.frame(DAE.CP_OMIM %>%
                                  group_by(across(phenotype)) %>%
                                  mutate(freq=n()) %>% 
                                  dplyr::select(c(5,6)) %>%
                                  unique() %>%
                                  mutate (type = "CP")) %>% arrange (-freq)
DAE.CP_OMIM_top40<-DAE.CP_OMIM_top40[1:40,]
DAE.GZ_OMIM<-read.delim("DAE.GZ_OMIM.phenotype.bed")

DAE.GZ_OMIM_top40 <- data.frame(DAE.GZ_OMIM %>%
                                  group_by(across(phenotype)) %>%
                                  mutate(freq=n()) %>% 
                                  dplyr::select(c(5,6)) %>%
                                  unique() %>%
                                  mutate (type = "GZ")) %>% arrange (-freq)
DAE.GZ_OMIM_top40<-DAE.GZ_OMIM_top40[1:40,]


DAE_OMIM<-merge(DAE.CP_OMIM_top40,DAE.GZ_OMIM_top40,'phenotype',all=T)

OMIM_CP <- DAE_OMIM %>% dplyr::select (c(1:3)) %>%
  mutate (freq.x = ifelse(is.na(freq.x), 0, freq.x)) %>%
  mutate (type.x = ifelse(is.na(type.x), "CP", type.x)) %>%
  arrange (freq.x) %>%
  mutate (ORDER = c(1:nrow(DAE_OMIM)))

OMIM_GZ <- DAE_OMIM %>% dplyr::select (c(1,4:5)) %>%
  mutate (freq.y = ifelse(is.na(freq.y), 0, freq.y)) %>%
  mutate (type.y = ifelse(is.na(type.y), "GZ", type.y)) %>%
  arrange (freq.y) %>%
  mutate (ORDER = OMIM_CP$ORDER [match (phenotype, OMIM_CP$phenotype)])

colnames(OMIM_GZ)[1:3] <- colnames(OMIM_CP)[1:3] <- c('phenotype','freq','type')

OMIM_CP.GZ<-rbind(OMIM_CP,OMIM_GZ)

ggplot(OMIM_CP.GZ, aes(y=reorder(phenotype, ORDER),x=freq, fill=type)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic() +
  scale_x_continuous(breaks=seq(0, 120, by = 30),expand = c(0, 0))+
  theme(axis.text.x = element_text(size=10,color = "black",angle=0, hjust = 0.3,vjust=1.5),
        axis.text.y =element_text(size=10,color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none")+
  scale_fill_manual(values = c("#FF0000","#0000FF"))

######
###### GWAS enrichment
######
# DAEs and nDAEs are intersected with brain related GWAS (downloaded from UCSC) using bedtools.
## data including enhancer coordinates, GWAS coordinates, GWAS phenotype
DAE_GWAS<-read.delim('DAE_GWAS.bed',header=F)
nDAE_GWAS<-read.delim('nDAE_GWAS.bed',header=F)

## Count number of DAE/nDAE enhancers for each GWAS phenotype
DAE.GWAS_group<- unique(data.frame(DAE_GWAS[c(1:3,7)] %>%
                                     group_by(across(V7)) %>%
                                     mutate(DAE.Freq=n())))
DAE.GWAS_group<-unique(DAE.GWAS_group[c(4,5)])
nDAE_GWAS_group<- unique(data.frame(nDAE_GWAS[c(1:3,7)] %>%
                                      group_by(across(V7)) %>%
                                      mutate(nDAE.Freq=n())))
nDAE_GWAS_group<-unique(nDAE_GWAS_group[c(4,5)])

## Merge DAE and nDAE GWAS together
DAE_nDAE_GWAS<-merge(DAE.GWAS_group,nDAE_GWAS_group,'V7',all=T)
## Calculate Odds Ratio
DAE_nDAE_GWAS[is.na(DAE_nDAE_GWAS)]<-0
DAE_nDAE_GWAS$DAE_Expect<- 39709- DAE_nDAE_GWAS$DAE.Freq
DAE_nDAE_GWAS$nDAE_Expect<- 162454- DAE_nDAE_GWAS$nDAE.Freq
## Add 0.5 value to all cells (based on the Haldane-Anscombe correction to adjust the odds ratio) to prevent getting INF value.
DAE_nDAE_GWAS$DAE.Freq<- DAE_nDAE_GWAS$DAE.Freq +0.5
DAE_nDAE_GWAS$nDAE.Freq<- DAE_nDAE_GWAS$nDAE.Freq +0.5
DAE_nDAE_GWAS$DAE_Expect<- DAE_nDAE_GWAS$DAE_Expect +0.5
DAE_nDAE_GWAS$nDAE_Expect<- DAE_nDAE_GWAS$nDAE_Expect +0.5
## calculate the OddsRatio
DAE_nDAE_GWAS$odd.Ratio<- log2((DAE_nDAE_GWAS$DAE.Freq * DAE_nDAE_GWAS$nDAE_Expect) / (DAE_nDAE_GWAS$nDAE.Freq * DAE_nDAE_GWAS$DAE_Expect))
DAE_nDAE_GWAS<-DAE_nDAE_GWAS[order(-DAE_nDAE_GWAS$odd.Ratio),]

ggplot(DAE_nDAE_GWAS[c(1:25),],aes(y=reorder(V7, odd.Ratio),x=odd.Ratio, V7)) +
  geom_point(size = 1,color='red') +
  theme_classic() +
  theme(axis.text.x = element_text(size=10,color = "black", hjust = 0.3,vjust =0.5),
        axis.text.y =element_text(size=9,color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none")+
  scale_color_manual(values = 'red')+
  geom_segment(aes(xend=0, yend=V7), color="gray")+
  scale_x_continuous(limits = c(0, 5),breaks = c(0,1,2,3,4,5),expand = c(0, 0.03))

######
###### Brandler and zhou enrichment
######

## Read data
## Data is prepared by intersecting between regions 
results <- read.delim('Brandler_results.csv',sep=',')
tabletext <- cbind(
  c("Regions","DAE", "DAE_CP", "DAE_GZ", "nDAE", "DAE"),
  c("pvalue","0.144", "0.003", "0.010", "0.268", "0.051"))## lable group and pvalue to plot

frdata <- results[,c("oddsratio","low_ci", "hi_ci")]
forestplot(tabletext,frdata, lwd.ci=1, boxsize=0.25, clip=c(0,12.2), zero=0, cex=0.9,
           is.summary = c(FALSE,rep(FALSE,8)), lty.ci = 1, xticks=c(seq(0,12,2)),lineheight=unit(0.4,'cm'),
           col=fpColors(box="black", lines="gray50", zero = 0),
           grid = structure(c(1), gp = gpar(lty = 2, col = "red")), 
           vertices=T, graph.pos = 2, colgap=unit(10,"mm"),
           txt_gp=fpTxtGp(label=gpar(cex=0.5),
                          ticks=gpar(cex=0.5),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1)))
### Read Zhou data
results <- read.delim('Zhou_result.csv',sep=',')
tabletext <- cbind(
  c("Regions","DAE", "DAE_CP_ASD", "DAE_GZ_ASD"),
  c("pvalue","0.157", "0.021", "0.746"))## lable group and pvalue to plot

frdata <- results[,c("oddsratio","low_ci", "hi_ci")]
forestplot(tabletext,frdata, lwd.ci=1, boxsize=0.25, clip=c(0,12.2), zero=0, cex=0.9,
           is.summary = c(FALSE,rep(FALSE,8)), lty.ci = 1, xticks=c(seq(0,12,2)),lineheight=unit(0.4,'cm'),
           col=fpColors(box="black", lines="gray50", zero = 0),
           grid = structure(c(1), gp = gpar(lty = 2, col = "red")), 
           vertices=T, graph.pos = 2, colgap=unit(10,"mm"),
           txt_gp=fpTxtGp(label=gpar(cex=0.5),
                          ticks=gpar(cex=0.5),
                          xlab=gpar(cex = 1),
                          title=gpar(cex = 1)))

######
###### Enhancer_Enhancer interaction
######
## Enhancer-Enhancer interaction table for CP
E_E_CP_data<-read.csv("E.DER_E.Non.DER_CP_enrich_OMIM.csv")[1:18,]
E_E_CP_data$No.Interaction[E_E_CP_data$No.Interaction > 5]<- 5
E_E_CP_data <- E_E_CP_data %>% 
  group_by(No.Interaction) %>%
  summarise(across(everything(), list(sum))) %>%
  mutate (ratio_OMIM_all =  DAE.linked.OMIM.gene_1 / DAE.linked.non.OMIM.gene_1,
          HiC = "CP") %>%
  dplyr::select (No.Interaction, ratio_OMIM_all, HiC)

## Enhancer-Enhancer interaction table for GZ
E_E_GZ_data<-read.csv("E.DER_E.Non.DER_GZ_enrich_OMIM.csv")[1:18,]
E_E_GZ_data$No.Interaction[E_E_GZ_data$No.Interaction > 5]<- 5
E_E_GZ_data <- E_E_GZ_data %>% 
  group_by(No.Interaction) %>%
  summarise(across(everything(), list(sum))) %>%
  mutate (ratio_OMIM_all =  DAE.linked.OMIM.gene_1 / DAE.linked.non.OMIM.gene_1,
          HiC = "GZ") %>%
  dplyr::select (No.Interaction, ratio_OMIM_all, HiC)

E_E_CP.GZ <- rbind (E_E_CP_data, E_E_GZ_data)

ggplot(data=E_E_CP.GZ, aes(x=No.Interaction, y=ratio_OMIM_all , group=HiC)) +
  geom_line(aes(color=HiC))+
  theme_classic()+
  scale_y_continuous("Ratio (OMIM/Non.OMIM genes)", expand = c(0, 0.1))+
  scale_x_continuous("Enhancer per DAE", expand = c(0, 0.1),breaks = 0:5)+
  coord_cartesian(ylim = c(0,1))+
  scale_color_manual(values=c("red", "blue")) + 
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'))


#####
##### Overlapping with BrainVar genes
#####
## Read Brain Var table that is downloaded from "Werling et al., 2020, Cell Reports"
brainvar <- read.csv("BrainVar.csv")
## Separating regions based on Trajectory groups
Falling <- brainvar %>% filter (Trajectory_group == "Falling")
Rising <- brainvar %>% filter (Trajectory_group == "Rising")
Non.transitional <- brainvar %>% filter (Trajectory_group == "Non-transitional")

dim(unique(Falling[1]))
dim(unique(Rising[1]))
dim(unique(Non.transitional[1]))

## Raed DAE and nDAE linked to CP genes
DAE.CP <- read.csv("DAE.CP_genes.csv")
nDAE.CP <- read.csv("specific_nDAE.CP_genes.csv")
## Merging DAE/nDAE and BrainVar data for Falling treajectory group
DAE.CP_Falling <- merge(DAE.CP, Falling,by.x="GeneID", by.y="EnsemblID")
nDAE.CP_Falling <- merge(nDAE.CP, Falling,by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.CP_Falling[1]))
dim(unique(nDAE.CP_Falling[1]))
## Calculating OddsRatio between DAE and nDAE for Falling treajectory group
GenesOddsRatio_DAE.nDAE_falling <- matrix(c(1092,787, 5946-1092, 4979-787), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_falling)
## Merging DAE/nDAE and BrainVar data for Rising treajectory group 
DAE.CP_Rising <- merge(DAE.CP, Rising, by.x="GeneID", by.y="EnsemblID")
nDAE.CP_Rising <- merge(nDAE.CP, Rising, by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.CP_Rising[1]))
dim(unique(nDAE.CP_Rising[1]))
## Calculating OddsRatio between DAE and nDAE for Rising treajectory group
GenesOddsRatio_DAE.nDAE_Rising <- matrix(c(1460,1175, 5946-1460, 4979-1175), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_Rising)
## Merging DAE/nDAE and BrainVar data for Non.transitional treajectory group
DAE.CP_Non.transitional <- merge(DAE.CP,Non.transitional, by.x="GeneID", by.y="EnsemblID")
nDAE.CP_Non.transitional <- merge(nDAE.CP,Non.transitional, by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.CP_Non.transitional[1]))
dim(unique(nDAE.CP_Non.transitional[1]))
## Calculating fisher.test between DAE and nDAE for Non.transitional treajectory group
GenesOddsRatio_DAE.nDAE_Non.transitional <- matrix(c(2129,1734, 5946-2129, 4979-1734), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_Non.transitional)
## making a table based on 95 percent confidence interval and odds ratio
OddsRatio_DAE.nDAE_CP <- data.frame(boxLabels=c("Falling", "Rising", "Constant"), 
                                    boxOdds=c(1.198289 , 1.053667 , 1.04382 ), 
                                    boxCILow=c(1.082422, 0.9638761 ,0.9639536 ), 
                                    boxCIHigh=c(1.326947, 1.1519268, 1.1303018),
                                    group="CP")

## Raed DAE and nDAE linked to GZ genes
DAE.GZ <- read.csv("DAE.GZ_genes.csv")
nDAE.GZ <- read.csv("specific_nDAE.GZ_genes.csv")
## Merging DAE/nDAE and BrainVar data for Falling treajectory group
DAE.GZ_Falling <- merge(DAE.GZ, Falling,by.x="GeneID", by.y="EnsemblID")
nDAE.GZ_Falling <- merge(nDAE.GZ, Falling,by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.GZ_Falling[1]))
dim(unique(nDAE.GZ_Falling[1]))
## Calculating OddsRatio between DAE and nDAE for Falling treajectory group
GenesOddsRatio_DAE.nDAE_Falling <- matrix(c(1128,821, 6085-1128, 5092-821), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_Falling)
## Merging DAE/nDAE and BrainVar data for Rising treajectory group
DAE.GZ_Rising <- merge(DAE.GZ, Rising, by.x="GeneID", by.y="EnsemblID")
nDAE.GZ_Rising <- merge(nDAE.GZ, Rising, by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.GZ_Rising[1]))
dim(unique(nDAE.GZ_Rising[1]))
## Calculating OddsRatio between DAE and nDAE for Rising treajectory group
GenesOddsRatio_DAE.nDAE_Rising <- matrix(c(1457,1199, 6085-1457, 5092-1199), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_Rising)
## Merging DAE/nDAE and BrainVar data for Non.transitional treajectory group
DAE.GZ_Non.transitional <- merge(DAE.GZ, Non.transitional, by.x="GeneID", by.y="EnsemblID")
nDAE.GZ_Non.transitional <- merge(nDAE.GZ, Non.transitional, by.x="GeneID", by.y="EnsemblID")
dim(unique(DAE.GZ_Non.transitional[1]))
dim(unique(nDAE.GZ_Non.transitional[1]))
## Calculating fisher.test between DAE and nDAE for Non.transitional treajectory group
GenesOddsRatio_DAE.nDAE_Non.transitional <- matrix(c(2209,1752, 6085-2209, 5092-1752), nrow=2)
fisher.test(GenesOddsRatio_DAE.nDAE_Non.transitional)
## making a table based on 95 percent confidence interval and odds ratio
OddsRatio_DAE.nDAE_GZ <- data.frame(boxLabels=c("Falling", "Rising", "Constant"), 
                                    boxOdds=c(1.183776, 1.022192, 1.086477), 
                                    boxCILow=c(1.071193, 0.9356652, 1.004262), 
                                    boxCIHigh=c(1.308504, 1.1168126, 1.175487),
                                    group="GZ")
## Combining DAE.nDAE_CP and DAE.nDAE_GZ tables together
OddsRatio_DAE.nDAE <- rbind(OddsRatio_DAE.nDAE_CP, OddsRatio_DAE.nDAE_GZ)

ggplot(OddsRatio_DAE.nDAE, aes(x=boxOdds, y=boxLabels, color=group)) + 
  geom_vline(aes(xintercept=1), size=0.1, linetype="dashed", color="darkgray") + 
  geom_errorbarh(aes(xmax=boxCIHigh, xmin=boxCILow), size=.5,
                 height=.05, position=position_dodge(.2)) +
  geom_point(size=1, position=position_dodge(.2)) +
  coord_trans(x=scales:::exp_trans(10)) +
  theme_classic() +
  ylab("") +
  xlab("Odds ratio\n(DAE/nDAE)") +
  theme(axis.text.x=element_text(color="black", angle=0, hjust=1, size=9),
        axis.text.y=element_text(color="black", angle=0, hjust=1, size=9),
        axis.title.x=element_text(size=9),legend.position = "none") +
  scale_color_manual(name="Target Gene", values=c("red","blue"))

#####
##### LDSC
#####
## Read output of LDSC function (from bash)
## Data including Name of traits and results of LDSC function
LDSC <- read.csv("LDSC.csv")
head(LDSC)
## Convert Z-score to p-value and adjuste significantly
dat <- LDSC %>%
  mutate(p=pnorm(Coefficient_z.score, lower.tail=FALSE)) %>%
  mutate(`-log10FDR`=-log10(p.adjust(p, "BH"))) %>% 
  dplyr::select(Category, Abbreviation, `-log10FDR`) %>% 
  spread(Abbreviation, `-log10FDR`) %>% 
  tibble::column_to_rownames("Category")

pheatmap(t(dat),
         col=colorRampPalette(brewer.pal(9, "OrRd"))(10), margins=c(5, 5),
         main="", cluster_rows=F, legend_labels="- log10 (q-value)",
         cluster_cols=F, angle_col=0,
         fontsize=12, fontsize_row=10, fontsize_col=10,
         border_color="grey60", number_color="black",
         display_numbers=T)

#####
##### HOMER
#####

library("readxl")
## Read each file for each bin group and select significant motifs
### "DAE_Homer.KnownMotif.xlsx" can be replaced to get plot for nDAE
bin1_5 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=1)
bin1_5 <- data.frame(bin1_5 %>% 
                       mutate(type="bin01_05") %>% 
                       filter(as.numeric(`P-value`) <= 1e-2) %>% 
                       dplyr::select(c(1:4,6,10)))
bin5_10 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=2)
bin5_10 <- data.frame(bin5_10 %>% 
                        mutate(type="bin05_10") %>% 
                        filter(as.numeric(`P-value`) <= 1e-2) %>% 
                        dplyr::select(c(1:4,6,10)))
bin10_15 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=3)
bin10_15 <- data.frame(bin10_15 %>% 
                         mutate(type="bin10_15") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin15_20 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=4)
bin15_20 <- data.frame(bin15_20 %>% 
                         mutate(type="bin15_20") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select (c(1:4,6,10)))
bin20_25 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=5)
bin20_25 <- data.frame(bin20_25 %>% 
                         mutate(type="bin20_25") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select (c(1:4,6,10)))
bin25_30 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=6)
bin25_30 <- data.frame(bin25_30 %>% 
                         mutate(type="bin25_30") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin30_35 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=7)
bin30_35 <- data.frame(bin30_35 %>% 
                         mutate(type="bin30_35") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin35_40 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=8)
bin35_40 <- data.frame(bin35_40 %>% 
                         mutate(type="bin35_40") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select (c(1:4,6,10)))
bin40_45 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=9)
bin40_45 <- data.frame(bin40_45 %>% 
                         mutate(type="bin40_45") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin45_50 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=10)
bin45_50 <- data.frame(bin45_50 %>% 
                         mutate(type="bin45_50") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin50_55 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=11)
bin50_55 <- data.frame(bin50_55 %>% 
                         mutate(type="bin50_55") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin55_60 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=12)
bin55_60 <- data.frame(bin55_60 %>% 
                         mutate(type="bin55_60") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin60_65 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=13)
bin60_65 <- data.frame(bin60_65 %>% 
                         mutate(type="bin60_65") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin65_70 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=14)
bin65_70 <- data.frame(bin65_70 %>% 
                         mutate(type="bin65_70") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin70_75 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=15)
bin70_75 <- data.frame(bin70_75 %>% 
                         mutate(type="bin70_75") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin75_80 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=16)
bin75_80 <- data.frame(bin75_80 %>% 
                         mutate(type="bin75_80") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin80_85 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=17)
bin80_85 <- data.frame(bin80_85 %>% 
                         mutate(type="bin80_85") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin85_90 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=18)
bin85_90 <- data.frame(bin85_90 %>% 
                         mutate(type="bin85_90") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin90_95 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=19)
bin90_95 <- data.frame(bin90_95 %>% 
                         mutate(type="bin90_95") %>% 
                         filter(as.numeric(`P-value`) <= 1e-2) %>% 
                         dplyr::select(c(1:4,6,10)))
bin95_100 <- read_excel("DAE_Homer.KnownMotif.xlsx", sheet=20)
bin95_100 <- data.frame(bin95_100 %>% 
                          mutate(type="bin95_100") %>% 
                          filter(as.numeric(`P-value`) <= 1e-2) %>% 
                          dplyr::select(c(1:4,6,10)))

## Change column name
for (files in ls(pattern="bin")) {
  FILE <- get(files)
  colnames(FILE)[5] <- "Target.Sequences.with.Motif"
  assign(files, FILE)
  rm(FILE)
} 
rm(files)

## Calculating number of motifs and number of target sequences with motif for bins
Number.motifs <- data.frame(bins=c("bin01_05","bin05_10","bin10_15","bin15_20","bin20_25","bin25_30","bin30_35","bin35_40",
                                   "bin40_45","bin45_50","bin50_55","bin55_60","bin60_65","bin65_70","bin70_75","bin75_80",
                                   "bin80_85","bin85_90","bin90_95","bin95_100"),
                            KnownMotifs=c(dim(bin1_5)[1],dim(bin5_10)[1],dim(bin10_15)[1],dim(bin15_20)[1],dim(bin20_25)[1],
                                          dim(bin25_30)[1],dim(bin30_35)[1],dim(bin35_40)[1],dim(bin40_45)[1],dim(bin45_50)[1],
                                          dim(bin50_55)[1],dim(bin55_60)[1],dim(bin60_65)[1],dim(bin65_70)[1],dim(bin70_75)[1],
                                          dim(bin75_80)[1],dim(bin80_85)[1],dim(bin85_90)[1],dim(bin90_95)[1],dim(bin95_100)[1]))

Target.seq <- data.frame(bins=c("bin01_05","bin05_10","bin10_15","bin15_20","bin20_25","bin25_30","bin30_35","bin35_40",
                                "bin40_45","bin45_50","bin50_55","bin55_60","bin60_65","bin65_70","bin70_75","bin75_80",
                                "bin80_85","bin85_90","bin90_95","bin95_100"),
                         KnownMotifs=c(sum(bin1_5[5]),sum(bin5_10[5]),sum(bin10_15[5]),sum(bin15_20[5]),sum(bin20_25[5]),
                                       sum(bin25_30[5]),sum(bin30_35[5]),sum(bin35_40[5]),sum(bin40_45[5]),sum(bin45_50[5]),
                                       sum(bin50_55[5]),sum(bin55_60[5]),sum(bin60_65[5]),sum(bin65_70[5]),sum(bin70_75[5]),
                                       sum(bin75_80[5]),sum(bin80_85[5]),sum(bin85_90[5]),sum(bin90_95[5]),sum(bin95_100[5])))

ggplot(Number.motifs, aes(x=bins, y=KnownMotifs )) +
  geom_bar(stat="identity", width=0.8, fill='darkgray') +
  theme_classic() +
  coord_cartesian(ylim=c(0, 170)) +    ## when plot for 'Number.motifs'
  scale_y_continuous(expand = c(0,0.1)) + 
  xlab("DAE Bins") + 
  ylab("Number of target sequences with motif") +
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(),
        axis.text.y=element_text(size=10, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=90, vjust=0.5, hjust=1))

ggplot(Target.seq, aes(x=bins, y=KnownMotifs )) +
  geom_bar(stat="identity", width=0.8, fill='darkgray') +
  theme_classic() +
  coord_cartesian(ylim = c(0, 4500000))+ ## when plot for 'Target.seq'
  scale_y_continuous(expand = c(0,0.1)) + 
  xlab("DAE Bins") + 
  ylab("Number of target sequences with motif") +
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(),
        axis.text.y=element_text(size=10, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=90, vjust=0.5, hjust=1))

## combining all motifs to select those ones that are not exist in all bins
Motif <- rbind(bin1_5,bin5_10,bin10_15,bin15_20,bin20_25,
               bin25_30,bin30_35,bin35_40,bin40_45,bin45_50,
               bin50_55,bin55_60,bin60_65,bin65_70,bin70_75,
               bin75_80,bin80_85,bin85_90,bin90_95,bin95_100)
## Selecting motif names that are not exist in all bins
NotCommonMotifs <- names(rowSums(table(Motif$Motif.Name, Motif$type))[rowSums(table(Motif$Motif.Name, Motif$type)) < 20])
## separating them from motif table
selected_motif <- Motif[Motif$Motif.Name %in% NotCommonMotifs,]
## convert log10pvalue to -log10pvale
selected_motif$Log.P.value <- - selected_motif$Log.P.value
## extracting name of the top 10 motifs to label in plot
name <- unique(selected_motif[selected_motif$type %in% "bin45_50" & selected_motif$Log.P.value > 250,])
order_selected_motif <- selected_motif[order(- selected_motif$Log.P.value),]
selected_motif$name <- gsub("/.*", "", name$Motif.Name)[match(selected_motif$Motif.Name, name$Motif.Name)]
selected_motif$name[selected_motif$type != "bin45_50"] <- "" 

library(ggplot2)
library(directlabels)
ggplot(selected_motif, aes(x=type, y=Log.P.value, group=Motif.Name, color=Motif.Name)) +
  geom_line(aes(color=Motif.Name), size=1) +
  theme_classic() +
  scale_y_continuous(expand=c(0, 0.1)) +
  coord_cartesian(ylim=c(0, 1100)) +
  xlab("DAE Bins") + ylab("-log 10 p-value") +
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(),
        axis.text.y=element_text(size=10, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=90, vjust=0.5, hjust=1),
        legend.position="none") +
  geom_dl(aes(label=name), method=list(dl.trans(x= x - 0.2), "top.points")) 

#####
##### Clustering DAE and nDAE for Adult K27ac
#####
## Data including DAE and nDAE coordinates and K27ac samples from fetal and adult brain
### Input needs to be replaced by other groups
H3K27ac <- read.delim("DAE_CP_H3K27ac_Fetal.Adult.bed")
#H3K27ac <- read.delim("DAE_GZ_H3K27ac_Fetal.Adult.bed")
#H3K27ac <- read.delim("nDAE_CP_H3K27ac_Fetal.Adult.bed")
#H3K27ac <- read.delim("nDAE_GZ_H3K27ac_Fetal.Adult.bed")

pam <- eclust(H3K27ac, "pam", k=4, hc_metric="spearman")## add k=5 for nDAEs groups
H3K27ac_pam <- cbind(H3K27ac, cluster=pam$clustering)

log_H3K27ac_pam <- log2(H3K27ac_pam[-c(15)] + 1)

## When we use nDAE regions the ,cluster5 = "yellow" will add to colours
my_colour <- list(cluster=c(cluster1="red", cluster2="blue",
                            cluster3="green", cluster4="purple"))
## the below colours will use when the Input is DAE_GZ
# my_colour <- list(cluster=c(cluster1="red", cluster2="green",
#                             cluster3="blue", cluster4="purple"))
pheatmap(log_H3K27ac_pam, annotation_colors=my_colour,
         color=colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100),
         show_rownames=F, cluster_cols=F, cluster_rows=T,
         annotation_legend=F, annotation_row=H3K27ac_pam[15])

#### epigenome pattern
H3K27ac_fetal.Adult <- H3K27ac_pam %>% group_by(cluster)%>% 
  summarise(
    CBC_119days = mean(CBC_119days),
    DFC_119days = mean(DFC_119days),
    CBC_133days = mean(CBC_133days),   
    DFC_133days = mean(DFC_133days),
    CBC_147days = mean(CBC_147days),
    DFC_147days = mean(DFC_147days),       
    CBC_268days = mean(CBC_268days),
    DFC_268days = mean(DFC_268days),
    CBC_23Y = mean(CBC_23Y),
    DFC_23Y = mean(DFC_23Y),
    CBC_30Y = mean(CBC_30Y),
    DFC_30Y = mean(DFC_30Y),     
    CBC_37Y  = mean(CBC_37Y),    
    DFC_37Y = mean(DFC_37Y))

H3K27ac_fetal.Adult_melt <- melt(H3K27ac_fetal.Adult, id="cluster")
H3K27ac_fetal.Adult_melt$log2 <- log2(H3K27ac_fetal.Adult_melt$value + 1)

ggplot(H3K27ac_fetal.Adult_melt, aes(x=variable, y=log2, group=cluster)) +
  geom_line(aes(color=factor(cluster)), size=0.3) +
  geom_point(aes(color=factor(cluster)), size=0.5) +
  theme_classic() + 
  ylab("log2 (Mean Normlized Count)") +
  coord_cartesian(ylim = c(0,6)) +
  theme(axis.text.x=element_text(size=9, color='black', angle=90, vjust=0.5, hjust=1),
        axis.title=element_blank(), axis.text.y=element_text(size=10, color='black'),
        legend.position="none")+
  scale_color_manual(values = c("red","blue","green","purple","yellow"))
## for DAE_GZ use :values = c("red","green","blue","purple","yellow")



