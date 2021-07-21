
###
###  Distrbution of nCER,GC and conservation scores based on 10bps bin
###

## Read enhancers with 10 bps bin linked to nCER (for DAEs and nDAEs)
## Data including initial coordinate of enhancers, 10 bps bin coordinates and nCER score
ncER_10bin <- read.delim('DAE_ncER_10bin.bed')
#ncER_10bin <- read.delim('nDAE_ncER_10bin.bed')
## Counting number of bins based on initial enhancer length
ncER_10bin <- data.frame(ncER_10bin %>% 
                           group_by(across(1:3)) %>%
                           mutate(numberofbin=n()) %>%
                           mutate(count=seq(1:n())))
## Defining relative position between 1 and 100 for bins.
ncER_10bin$relative.Pos <- ifelse (ncER_10bin$count == 1, 1,
                                   ifelse (ncER_10bin$count == ncER_10bin$numberofbin, 100,
                                           (1 + ((ncER_10bin$count - 1)*(99/(ncER_10bin$numberofbin -1 ))))))

## Making coordinate column in order to match GC and conservation scores
ncER_10bin$coordinate <- paste0(ncER_10bin$chr.bin, ":", ncER_10bin$start.bin, "-", ncER_10bin$end.bin)
## Adding GC and Conservation scores to each bin
GC_10bin <- data.frame(GCcontent(Hsapiens,GRanges(ncER_10bin$coordinate), as.prob=T))
phastcons_10bin <- as.data.frame(gscores(phastCons100way.UCSC.hg19, GRanges(ncER_10bin$coordinate)))
ncER_GC_phastcons_10bin <- cbind(ncER_10bin, GC_10bin, phastcons_10bin[6])

## calculating mean of each score for bins with the same positions
BinToPlot <- data.frame(ncER_GC_phastcons_10bin %>%
                          group_by(relative.Pos) %>% 
                          summarize(nEnha = n(), 
                                   mean_ncER = mean(ncER),
                                    mean_GC = mean(C.G),
                                    mean_phastcons = mean(na.omit(default))))

## Drawing plot 
ggplot(data=BinToPlot) + aes(x=relative.Pos) +
  stat_smooth(aes(x=relative.Pos, y=mean_ncER), 
              method="lm", formula=y ~ poly(x, 21), se=F) +
  ylab("nCER percentile") +
  xlab("Relative position\nof each bin") +
  update_geom_defaults("smooth", list(size=.8)) +
  theme_classic() +
  coord_cartesian(ylim = c(50,65))+
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_text(),
        axis.title.y=element_text())

######
###### Enrichment of epigenome data across 100 bins
######
## Read epigenome data for each 10bps bin.
## Data including initial coordinate of enhancers,10 bps bin coordinates and epigenome data
epigenom_10bins <- read.delim('DAE_epigenome_10bins.bed', header=T)
## Counting the number of bins based on initial enhancer length
epigenom_10bins_group <- data.frame(epigenom_10bins %>%
                                      group_by(across(1:3)) %>%
                                      mutate(numberofbin = n()) %>%
                                      mutate(count = seq(1:n())))
## Defining relative position between 1 and 100 for bins.
epigenom_10bins_group$relative.Pos <- ifelse (epigenom_10bins_group$count == 1, 1,
                                            ifelse (epigenom_10bins_group$count == epigenom_10bins_group$numberofbin, 100,
                                                    (1 + ((epigenom_10bins_group$count - 1)*(99/(epigenom_10bins_group$numberofbin -1 ))))))
## calculating mean of each epigenome sample for bins with the same positions
BinToPlot <- data.frame(epigenom_10bins_group %>% group_by(relative.Pos)  %>%
                          summarize(nEnha = n (),
                                    ATAC_Cp.day105 = mean (ATAC_Cp.day105),
                                    ATAC_GZ.day105 = mean (ATAC_GZ.day105),
                                    ATAC_Cp.day112 = mean (ATAC_Cp.day112),
                                    ATAC_GZ.day112 = mean (ATAC_GZ.day112),
                                    ATAC_Cp.day119 = mean (ATAC_Cp.day119),
                                    ATAC_GZ.day119 = mean (ATAC_GZ.day119),
                                    DNase_day.85 = mean (DNase_day.85),
                                    DNase_96days = mean (DNase_96days),
                                    DNase_day.101 = mean (DNase_day.101),
                                    DNase_day.104 = mean (DNase_day.104),
                                    DNase_day.109 = mean (DNase_day.109),
                                    DNase_day.117 = mean (DNase_day.117),
                                    DNase_day.122 = mean (DNase_day.122),
                                    DNase_day.142 = mean (DNase_day.142),
                                    H3K27ac_day49 = mean (H3K27ac_day49),
                                    H3K27ac_day60 = mean (H3K27ac_day60),
                                    H3K27ac_O_day84 = mean (H3K27ac_O_day84),
                                    H3K27ac_F_day84 = mean (H3K27ac_F_day84),
                                    H3K27ac_CP.day109 = mean (H3K27ac_CP.day109),
                                    H3K27ac_GZ.day109 = mean (H3K27ac_GZ.day109),
                                    H3K27ac_CP.day113 = mean (H3K27ac_CP.day113),
                                    H3K27ac_GZ.day113 = mean (H3K27ac_GZ.day113),
                                    H3K27ac_CBC.day119 = mean (H3K27ac_CBC.day119),
                                    H3K27ac_DFC.day119 = mean (H3K27ac_DFC.day119),
                                    H3K27ac_CP.day123 = mean (H3K27ac_CP.day123),
                                    H3K27ac_GZ.day123 = mean (H3K27ac_GZ.day123),
                                    H3K27ac_CBC.day133 = mean (H3K27ac_CBC.day133),
                                    H3K27ac_DFC.day133 = mean (H3K27ac_DFC.day133),
                                    H3K27ac_CBC.day268 = mean (H3K27ac_CBC.day268),
                                    H3K27ac_DFC.day268 = mean (H3K27ac_DFC.day268),
                                    H3K27ac_CBC.day147 = mean (H3K27ac_CBC.day147),
                                    H3K27ac_DFC.day147 = mean (H3K27ac_DFC.day147),
                                    H3K4me1_day119.120 = mean (H3K4me1_day119.120),
                                    H3K4me2_day49 = mean (H3K4me2_day49),
                                    H3K4me2_day60= mean (H3K4me2_day60),
                                    H3K4me2_O_day84= mean (H3K4me2_O_day84),
                                    H3K4me2_F_day84= mean (H3K4me2_F_day84),
                                    K4me3_CP.day109= mean (K4me3_CP.day109),
                                    K4me3_GZ.day109= mean (K4me3_GZ.day109),
                                    K4me3_CP.day113= mean (K4me3_CP.day113),
                                    K4me3_GZ.day113= mean (K4me3_GZ.day113),
                                    K4me3_CP.day123= mean (K4me3_CP.day123),
                                    K4me3_GZ.day123= mean (K4me3_GZ.day123),
                                    K4me3_CBC.day147= mean (K4me3_CBC.day147),
                                    K4me3_DFC.day147= mean (K4me3_DFC.day147),
                                    K27me3_CP.day109= mean (K27me3_CP.day109),
                                    K27me3_GZ.day109= mean (K27me3_GZ.day109),
                                    K27me3_CP.day113= mean (K27me3_CP.day113),
                                    K27me3_GZ.day113= mean (K27me3_GZ.day113),
                                    K27me3_CP.day123= mean (K27me3_CP.day123),
                                    K27me3_GZ.day123= mean (K27me3_GZ.day123)))

## Drawing plot for each epigenome mark
ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_CP.day109), color='#9999FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_GZ.day109), color='#0000CC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_CP.day113), color='#FF9999', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_GZ.day113), color='#FF0000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_CP.day123), color='#33FF33', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K27me3_GZ.day123), color='#00FF00', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_CP.day109), color='#9999FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_GZ.day109), color='#0000CC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_CP.day113), color='#FF9999', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_GZ.day113), color='#FF0000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_CP.day123), color='#33FF33', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=K4me3_GZ.day123), color='#00FF00', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  theme(axis.text.x= element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10,color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=H3K4me2_day49), color='#CCCC00', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K4me2_day60), color='#CCCC00',  method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K4me2_O_day84), color='#CCCCCC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K4me2_F_day84), color='#666666', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  scale_y_continuous(breaks=c(0.3, 0.8, 1.3, 1.8)) +
  coord_cartesian(ylim=c(0.3, 2)) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=H3K4me1_day119.120), color='#006600', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  scale_y_continuous(breaks=c(0.15,0.18,0.21,0.24,0.27)) +
  coord_cartesian(ylim=c(0.15, 0.27))+
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_day49), color='#CCCC00', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_day60), color='#CCCC00', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_O_day84), color='#CCCCCC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_F_day84), color='#666666', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CP.day109), color='#9999FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_GZ.day109), color='#0000CC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CP.day113), color='#FF9999', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_GZ.day113), color='#FF0000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CBC.day119), color='#66FF66', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_DFC.day119), color='#99FF99', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CP.day123), color='#33FF33', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_GZ.day123), color='#00FF00', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CBC.day133), color='#00CCCC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_DFC.day133), color='#66FFFF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CBC.day268), color='#FF8000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_DFC.day268), color='#FFB366', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_CBC.day147), color='#FF00FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=H3K27ac_DFC.day147), color='#FF99FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2)) +
  coord_cartesian(ylim=c(0, 2.5)) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3,vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_Cp.day105), color='#9999FF', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_GZ.day105), color='#0000CC', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_Cp.day112), color='#FF9999', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_GZ.day112), color='#FF0000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_Cp.day119), color='#33FF33', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=ATAC_GZ.day119), color='#00FF00', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  coord_cartesian(ylim=c(0.5, 2.1)) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3, vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")

ggplot(data=BinToPlot, aes(x=relative.Pos)) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.85), color='#000000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_96days), color='#000000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.101), color='#000000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.104), color='#000000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.109), color='#000066', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.117), color='#660000', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.122), color='#006600', method="lm", formula=y ~ poly(x, 21), se=F) +
  stat_smooth(aes(x=relative.Pos, y=DNase_day.142), color='#660066', method="lm", formula=y ~ poly(x, 21), se=F) +
  theme_classic() + update_geom_defaults("smooth", list(size=.7)) +
  scale_y_continuous(breaks=c(0.2,0.5,0.8,1.1)) +
  coord_cartesian(ylim=c(0.2, 1.4)) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.3,vjust=0.5),
        axis.text.y=element_text(size=10, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")
