
######
###### Clustring DAEs across different epigenome data for three time points
######
## Read normalized epigenoem data for DAEs linked to CP or GZ genes
## Data including DAE coordinates and different normalized epigenome data
Input<-read.delim("DAE_CP_epigenome.Normalized.bed")
#Input<-read.delim("DAE_GZ_epigenome.Normalized.bed")

## Separating epigenome data based on time-points
range_8_12week <- Input[c(7:11,20:23,48:51)]
range_13_17week <- Input[c(1:6,12:18,20,21,24:29,36:45,52:65,70:83)]
range_18_moreweek <- Input[c(19,30:35,46,47,66:69,84,85)]

### Clustring epigenome data
pam_8_12week <- eclust(range_8_12week,"pam", k=2, "spearman") 
Epigenome_cluster_8_12week <- cbind(range_8_12week, cluster=pam_8_12week$clustering)

pam_13_17week <- eclust(range_13_17week,"pam", k=3, "spearman")
Epigenome_cluster_13_17week <- cbind(range_13_17week, cluster=pam_13_17week$clustering)

pam_18_moreweek <- eclust(range_18_moreweek,"pam", k=2, "spearman") 
Epigenome_cluster_18_moreweek <- cbind(range_18_moreweek, cluster=pam_18_moreweek$clustering)

## Average of epigenome data for samples with replicates for range_8_12week
Epigenome_cluster_8_12week$D_56.58 <- rowMeans(Epigenome_cluster_8_12week[c(1,2)])
Epigenome_cluster_8_12week$D_72.76 <- rowMeans(Epigenome_cluster_8_12week[c(3,4)])
Epigenome_cluster_8_12week$D_85 <- Epigenome_cluster_8_12week$DNase_day85
Epigenome_cluster_8_12week$ac_D_49 <- Epigenome_cluster_8_12week$H3K27ac_day49
Epigenome_cluster_8_12week$ac_D_60 <- Epigenome_cluster_8_12week$H3K27ac_day60
Epigenome_cluster_8_12week$ac_D_84 <- rowMeans(Epigenome_cluster_8_12week[c(8,9)])
Epigenome_cluster_8_12week$me2_D_49 <- Epigenome_cluster_8_12week$H3k4me2_day49
Epigenome_cluster_8_12week$me2_D_60 <- Epigenome_cluster_8_12week$H3k4me2_day60
Epigenome_cluster_8_12week$me2_D_84 <- rowMeans(Epigenome_cluster_8_12week[c(12,13)])

## Average of epigenome data for samples with replicates for range_13_17week
Epigenome_cluster_13_17week$ATAC_D_105 <- rowMeans(Epigenome_cluster_13_17week[c(1,2)])
Epigenome_cluster_13_17week$ATAC_D_112 <- rowMeans(Epigenome_cluster_13_17week[c(3,4)])
Epigenome_cluster_13_17week$ATAC_D_119 <- rowMeans(Epigenome_cluster_13_17week[c(5,6)])
Epigenome_cluster_13_17week$D_96 <- Epigenome_cluster_13_17week$DNase_day96
Epigenome_cluster_13_17week$D_101 <- Epigenome_cluster_13_17week$DNase_day101
Epigenome_cluster_13_17week$D_104 <- Epigenome_cluster_13_17week$DNase_day104
Epigenome_cluster_13_17week$D_105 <- Epigenome_cluster_13_17week$DNase_day105
Epigenome_cluster_13_17week$D_109 <- Epigenome_cluster_13_17week$DNase_day109
Epigenome_cluster_13_17week$D_117 <- Epigenome_cluster_13_17week$DNase_day117
Epigenome_cluster_13_17week$D_122 <- Epigenome_cluster_13_17week$DNase_day122
Epigenome_cluster_13_17week$ac_D_105 <- rowMeans(Epigenome_cluster_13_17week[c(22,23)])
Epigenome_cluster_13_17week$ac_D_109 <- rowMeans(Epigenome_cluster_13_17week[c(16,17)])
Epigenome_cluster_13_17week$ac_D_113 <- rowMeans(Epigenome_cluster_13_17week[c(18,19)])
Epigenome_cluster_13_17week$ac_D_123 <- rowMeans(Epigenome_cluster_13_17week[c(20,21)])
Epigenome_cluster_13_17week$me1_D_105 <- Epigenome_cluster_13_17week$H3K4me1_GENSPC_HuFNSC04_day105
Epigenome_cluster_13_17week$me1_D_119 <- rowMeans(Epigenome_cluster_13_17week[c(25,26,28:31)])
Epigenome_cluster_13_17week$me1_D_120 <- Epigenome_cluster_13_17week$H3K4me1_H22676_day120
Epigenome_cluster_13_17week$me3_D_105 <- rowMeans(Epigenome_cluster_13_17week[c(32:34)])
Epigenome_cluster_13_17week$me3_D_109 <- rowMeans(Epigenome_cluster_13_17week[c(35,36)])
Epigenome_cluster_13_17week$me3_D_113 <- rowMeans(Epigenome_cluster_13_17week[c(37,38)])
Epigenome_cluster_13_17week$me3_D_119 <- rowMeans(Epigenome_cluster_13_17week[c(39,42)])
Epigenome_cluster_13_17week$me3_D_122 <- Epigenome_cluster_13_17week$H3K4me3_H.22510_day122
Epigenome_cluster_13_17week$me3_D_123 <- rowMeans(Epigenome_cluster_13_17week[c(44,45)])
Epigenome_cluster_13_17week$me27_D_105 <- rowMeans(Epigenome_cluster_13_17week[c(46,47)])
Epigenome_cluster_13_17week$me27_D_109 <- rowMeans(Epigenome_cluster_13_17week[c(48,49)])
Epigenome_cluster_13_17week$me27_D_113 <- rowMeans(Epigenome_cluster_13_17week[c(50,51)])
Epigenome_cluster_13_17week$me27_D_119 <- rowMeans(Epigenome_cluster_13_17week[c(52:55)])
Epigenome_cluster_13_17week$me27_D_120 <- Epigenome_cluster_13_17week$H3K27me3_H22676_day120
Epigenome_cluster_13_17week$me27_D_122 <- Epigenome_cluster_13_17week$H3K27me3_H22510_day122
Epigenome_cluster_13_17week$me27_D_123 <- rowMeans(Epigenome_cluster_13_17week[c(58,59)])

## Average of epigenome data for samples with replicates for range_18_moreweek
Epigenome_cluster_18_moreweek$D_142 <- Epigenome_cluster_18_moreweek$DNase_day142
Epigenome_cluster_18_moreweek$ac_D_133 <- rowMeans(Epigenome_cluster_18_moreweek[c(2,3)])
Epigenome_cluster_18_moreweek$ac_D_147 <- rowMeans(Epigenome_cluster_18_moreweek[c(4,5)])
Epigenome_cluster_18_moreweek$ac_D_268 <- rowMeans(Epigenome_cluster_18_moreweek[c(6,7)])
Epigenome_cluster_18_moreweek$me1_D_140 <- rowMeans(Epigenome_cluster_18_moreweek[c(8,9)])
Epigenome_cluster_18_moreweek$me3_D_140 <- rowMeans(Epigenome_cluster_18_moreweek[c(10,11)])
Epigenome_cluster_18_moreweek$me3_D_147 <- rowMeans(Epigenome_cluster_18_moreweek[c(12,13)])
Epigenome_cluster_18_moreweek$me27_D_140 <- rowMeans(Epigenome_cluster_18_moreweek[c(14,15)])

## Selecting prepared columns to plot
Epigenome_cluster_Range_8_12week <- Epigenome_cluster_8_12week[14:23]#range_8_12week
Epigenome_cluster_Range_13_17week <- Epigenome_cluster_13_17week[60:90]#range_13_17week
Epigenome_cluster_Range_18_moreweek <- Epigenome_cluster_18_moreweek[16:24]#range_18_moreweek

## We need to use each time point as Input "Epigenome_cluster_Range_13_17week","Epigenome_cluster_Range_18_moreweek" separately
log_Epigenome_cluster <- log2(Epigenome_cluster_Range_8_12week[-1]+1)
my_colour <- list(cluster=c(cluster1 = "#CC0000", cluster2 = "#00CC00"))# for range_13_17week we need to cluster3 = "#0066CC"
pheatmap(log_Epigenome_cluster,
         annotation_colors=my_colour,
         color=colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100),
         show_rownames=F, clustering_distance_cols=as.dist(1 - cols.cor),
         cluster_cols=F, cluster_rows=T,
         fontsize_row=10, fontsize_col=9, treeheight_row=0,annotation_legend=F,
         annotation_row=Epigenome_cluster_Range_8_12week[1])## the "annotation_row=" needs to be changed based on time-point

### Epigenome pattern
## Group epigenome data for each defined cluster

###range_8_12week
Epigenome_pattern_8_12week <- Epigenome_cluster_8_12week %>%
  group_by(cluster) %>%
  summarise(
    D_56.58 = mean(D_56.58),
    D_72.76 = mean(D_72.76),
    D_85 = mean(D_85),
    ac_D_49 = mean(ac_D_49),
    ac_D_60 = mean(ac_D_60),
    ac_D_84 = mean(ac_D_84),
    me2_D_49= mean(me2_D_49),
    me2_D_60 = mean(me2_D_60),
    me2_D_84 = mean(me2_D_84))

###range_13_17week
Epigenome_pattern_13_17week <- Epigenome_cluster_13_17week %>%
  group_by(cluster) %>%
  summarise(
    ATAC_D_105 = mean(ATAC_D_105),
    ATAC_D_112 = mean(ATAC_D_112),
    ATAC_D_119 = mean(ATAC_D_119),
    D_96 = mean(D_96),
    D_101 = mean(D_101),
    D_104 = mean(D_104),
    D_105 = mean(D_105),
    D_109 = mean(D_109),
    D_117 = mean(D_117),
    D_122 = mean(D_122),
    ac_D_105 = mean(ac_D_105),
    ac_D_109 = mean(ac_D_109),
    ac_D_113 = mean(ac_D_113),
    ac_D_123 = mean(ac_D_123),
    me1_D_105 = mean(me1_D_105),
    me1_D_119 = mean(me1_D_119),
    me1_D_120 = mean(me1_D_120),
    me3_D_105 = mean(me3_D_105),
    me3_D_109 = mean(me3_D_109),
    me3_D_113 = mean(me3_D_113),
    me3_D_119 = mean(me3_D_119),
    me3_D_122 = mean(me3_D_122),
    me3_D_123 = mean(me3_D_123),
    me27_D_105 = mean(me27_D_105),
    me27_D_109 = mean(me27_D_109),
    me27_D_113 = mean(me27_D_113),
    me27_D_119 = mean(me27_D_119),
    me27_D_120 = mean(me27_D_120),
    me27_D_122 = mean(me27_D_122),
    me27_D_123 = mean(me27_D_123))

###range_18_moreweek
Epigenome_pattern_18_moreweek<- Epigenome_cluster_18_moreweek %>%
  group_by(cluster) %>%
  summarise(
    D_142 = mean(D_142),
    ac_D_133 = mean(ac_D_133),
    ac_D_147 = mean(ac_D_147),
    ac_D_268 = mean(ac_D_268),
    me1_D_140 = mean(me1_D_140),
    me3_D_140 = mean(me3_D_140),
    me3_D_147 = mean(me3_D_147),
    me27_D_140 = mean(me27_D_140))

## The Input needs to be changed based on time-points
epigenome_pattern <- melt(Epigenome_pattern_8_12week, id="cluster")
epigenome_pattern$log2 <- log2(epigenome_pattern$value)

ggplot(epigenome_pattern, aes(x=variable, y=log2, group=cluster)) +
  geom_line(aes(color=factor(cluster)), size=0.3) +
  geom_point(aes(color=factor(cluster)), size=0.5) +
  theme_classic() +
  ylab("log2 (Mean Normlized Count)") +
  coord_cartesian(ylim=c(2, 9)) +
  theme(axis.text.x=element_text(size=9, color='black', angle=90, vjust=0.5, hjust=1),
        axis.title=element_blank(), axis.text.y=element_text(size=10, color='black'),
        legend.position="none") +
  scale_color_manual(values=c("#CC0000", "#00CC00","#0066CC"))

### Expression pattern
## Read gene expression (from AllenBrain) data for DAEs and nDAEs linked to CP or GZ genes
Expression <- read.delim('DAE.nDAE_Gene.expression_AllenBrain.bed')

## Separating related DAEs_CP/GZ from expression data
### The Input needs to be changed based on time-points
Epigenome_cluster_8_12week$coordinate <- row.names(Epigenome_cluster_8_12week)
Epigenome_cluster_13_17week$coordinate <- row.names(Epigenome_cluster_13_17week)
Epigenome_cluster_18_moreweek$coordinate <- row.names(Epigenome_cluster_18_moreweek)

Expression_pattern_cluster_8_12week <- merge(Expression, Epigenome_cluster_8_12week, 'coordinate')
Expression_pattern_cluster_13_17week <- merge(Expression, Epigenome_cluster_13_17week, 'coordinate')
Expression_pattern_cluster_18_moreweek <- merge(Expression, Epigenome_cluster_18_moreweek, 'coordinate')

## Select time-points for plotting
Expression_pattern_cluster_8_12week <- unique(Expression_pattern_cluster_8_12week[, c(9:11,35)])#range_8_12week
Expression_pattern_cluster_13_17week <- unique(Expression_pattern_cluster_13_17week[,c(12:14,81)])#range_13_17week
Expression_pattern_cluster_18_moreweek <- unique(Expression_pattern_cluster_18_moreweek[,c(15:21,37)])#range_18_moreweek

## preparing extracted data to plot
### The Input needs to be changed based on time-points
Expression_pattern <- melt(Expression_pattern_cluster_8_12week, id="cluster", na.rm=T) 
Expression_pattern <- Expression_pattern[Expression_pattern$value > 0.5,]
Expression_pattern$log2 <- log2(Expression_pattern$value + 1)
Expression_pattern$cluster <- as.factor(Expression_pattern$cluster)
Expression_pattern$cluster <- factor(Expression_pattern$cluster,
                                     levels=c('1','3','2'), ordered=TRUE)## this is only for range_13_18 to order gene expression boxplot

ggplot(Expression_pattern, aes(x=variable, y=log2, fill=cluster)) +
  geom_boxplot(width=0.8, outlier.colour="black", outlier.shape=16, outlier.size=0.4, 
               alpha=0.7, position=position_dodge(width=0.9)) +
  scale_y_continuous("Log2(RPKM +1)") +
  theme_classic() + coord_cartesian(ylim=c(0,9)) +
  theme(axis.text.x=element_text(size=10, color='black', angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("#CC0000","#00CC00"))# For range_13_18 -> values=c("#CC0000","#0066CC", "#00CC00")


