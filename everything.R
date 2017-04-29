setwd('~/FS2/16s/V4/')

library(vegan)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(reshape2)


set.seed(7777777)

### functions ###

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

#### end functions  ######


#### reading in the data ####

tax <- read.table('V4.final.taxonomy', header = TRUE, stringsAsFactors = FALSE)
shared <- read.table('V4.final.shared', header=TRUE, stringsAsFactors = FALSE)
rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

### NTC stuff ###

ntc <- shared[rownames(shared) == 'NTC',]
rowSums(ntc) # only 71 sequences in my NTC, not bad I'd say
ntc[,which(ntc >1)]
# this could be due to index read errors or low level contamination, either way there are very few sequences so I feel fine with this
### end NTC ###

shared <- shared[rownames(shared) != 'NTC',]

# metadata whatnot #

meta <- as.data.frame(rownames(shared))
colnames(meta) <- 'group'
meta$experiment <- NA
meta$tissue <- NA
meta$day <- NA
meta$pig_num <- NA
meta$treatment <- NA


FS2.cont <- c(67,68,69,70,71,72,73,81,82,83,84,85,86,87)
FS2.RPS <- c(74,75,76,77,78,79,80,90,91,92,93,94,95,96)

FS4.cont <- c(19,4,13,22,25,2)
FS4.tet <- c(21,1,27,10,6,5)
FS4.RPS <- c(11,18,24,26,28,23)
FS4.in <- c(15,3,16,7,9,30)
FS4.but <- c(14,20,29,8,12,17)

FS5.car <- c(38,39,43)
FS5.SCID <- c(40,41,42,44)

meta$experiment[grep('X2', meta$group)] <- 'FS2'
meta$experiment[grep('X3', meta$group)] <- 'FS4'
meta$experiment[grep('X5', meta$group)] <- 'FS5'

meta$pig_num <- factor(as.numeric(gsub('X.P([0-9]+).*', '\\1', meta$group)))
meta$day <- factor(as.numeric(gsub('.*D([0-9]+)', '\\1', meta$group)))
meta$tissue <- gsub('X[0-9]P[0-9]+.*([A-Za-z])D[0-9]+', '\\1', meta$group)

meta$tissue[meta$tissue == 'F'] <- 'feces'
meta$tissue[meta$tissue == 'C'] <- 'colon'
meta$tissue[meta$tissue %in% c('i', 'I')] <- 'ileum'
meta$tissue[meta$tissue %in% c('X', 'E')] <- 'cecum'
meta$tissue[meta$tissue == 'U'] <- 'cec_cont_RNA'
meta$tissue[meta$tissue == 'M'] <- 'cec_cont_DNA'

meta$treatment[meta$experiment == 'FS2' & meta$pig_num %in% FS2.cont] <- 'control'
meta$treatment[meta$experiment == 'FS2' & meta$pig_num %in% FS2.RPS] <- 'RPS'

meta$treatment[meta$experiment == 'FS4' & meta$pig_num %in% FS4.cont] <- 'control'
meta$treatment[meta$experiment == 'FS4' & meta$pig_num %in% FS4.but] <- 'butyrate'
meta$treatment[meta$experiment == 'FS4' & meta$pig_num %in% FS4.in] <- 'inulin'
meta$treatment[meta$experiment == 'FS4' & meta$pig_num %in% FS4.RPS] <- 'RPS'
meta$treatment[meta$experiment == 'FS4' & meta$pig_num %in% FS4.tet] <- 'tet'

meta$treatment[meta$experiment == 'FS5' & meta$pig_num %in% FS5.car] <- 'carrier'
meta$treatment[meta$experiment == 'FS5' & meta$pig_num %in% FS5.SCID] <- 'SCID'


meta$num_seqs <- rowSums(shared)
meta$day[meta$day == 42] <- 41

meta[meta$pig_num %in% c(83,87,93) & meta$tissue != 'feces', ]$day <- 41 # this is because ELI misnamed some of the samples. Jen noticed and is a hero.

meta$design <- paste(meta$experiment, meta$tissue, meta$day, meta$treatment, sep = '_')
meta$treatXday <- paste(meta$treatment, meta$day, sep = '_day_')
design <- meta[,c(1,8)]

FS2.accnos <- (meta$group[meta$experiment == 'FS2'])

rownames(shared) == meta$group # just checking...


# writing out meta and shared for correlation #

meta_for_corr <- filter(meta, experiment == 'FS2')
meta_for_corr$sample <- paste(meta_for_corr$tissue, meta_for_corr$day, meta_for_corr$pig_num, sep = '_')
#write.table(meta_for_corr, '~/FS2/correlation/but_corr/16S_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
shared2 <- shared[meta$experiment == 'FS2',]
shared2$group <- rownames(shared2)
write.table(shared2, '~/FS2/correlation/but_corr/16S_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)



#write.table(design, file = 'V4.design', quote = FALSE, sep = '\t', row.names = FALSE)
#write.table(FS2.accnos, file = 'FS2.accnos', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
#write.table(meta, file = 'V4.metadata.txt', quote = FALSE, sep = '\t', row.names = FALSE)
# I wrote these out so I can load them later or in other scripts so I dont have to run this chunk constantly

######### SHITS FUCKED UP HERE ##############
######### IS SHIT REALLY FUCKED UP HERE?? ###########
################# alpha diversity calcs ##########

FS2 <- shared[meta$experiment == "FS2",]
FS2.meta <- meta[meta$experiment == 'FS2',]
FS2 <- FS2[FS2.meta$tissue != 'cec_cont_DNA',]
FS2.meta <- FS2.meta[FS2.meta$tissue != 'cec_cont_DNA',]

#write.table(FS2, '~/FS2/correlation/but_corr/FS216S.shared', quote = FALSE, sep = '\t', row.names = TRUE)
colnames(FS2)[2070]

FS2 <- FS2[FS2.meta$day %in% c(0:21),]
FS2.meta <- FS2.meta[FS2.meta$day %in% c(0:21),]

FS2.meta$design <- as.factor(FS2.meta$design )

FecalTime <- FS2[FS2.meta$tissue == 'feces',]

FecalTime.meta <- FS2.meta[FS2.meta$tissue == 'feces',]

FecalTime <- FecalTime[rowSums(FecalTime) > 4200,]
FecalTime.meta <- FecalTime.meta[FecalTime.meta$group %in% rownames(FecalTime),]
FecalTime <- rrarefy(FecalTime, min(rowSums(FecalTime)))
rowSums(FecalTime)
FecalTime.meta$invsimpson <- diversity(FecalTime, index = 'invsimpson')
FecalTime.meta$shannon <- diversity(FecalTime, index = 'shannon')

#filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
#  ggplot() + geom_boxplot(aes(day, shannon, fill = treatment)) +
#  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time')


### USE THIS FIG FOR FECES ALPHA DIV ###
p1 <- filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time') +
  annotate('segment', x=4.7, xend=5.2, y=30, yend=30, size=.3) + annotate('text', x=4.9, y=31, label='* p=0.044', size=4)
###

## nifty little thing to do wilcoxon tests on alpha diversity at each day between the two treatments

fecal_alpha_wilcox_results <- FecalTime.meta %>% group_by(day) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(day, Wilcox = wilco$p.value)




Tissue21 <- FS2[FS2.meta$day ==21 & FS2.meta$tissue %in% c('ileum', 'cecum', 'colon', 'cec_cont_RNA'),]
Tissue21.meta <- FS2.meta[FS2.meta$day ==21 & FS2.meta$tissue %in% c('ileum', 'cecum', 'colon', 'cec_cont_RNA'),]

Tissue21 <- Tissue21[rowSums(Tissue21) > 1000,]  # removes samples with less than 1000 seqs
Tissue21 <- rrarefy(Tissue21, min(rowSums(Tissue21)))
Tissue21.meta <- Tissue21.meta[Tissue21.meta$group %in% rownames(Tissue21),]
rowSums(Tissue21)
Tissue21.meta$invsimpson <- diversity(Tissue21, index = 'invsimpson')


# wilcoxon tests for tissues
tissue_alpha_wilcox_results <- Tissue21.meta %>% group_by(tissue) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(tissue, Wilcox = wilco$p.value)

# use this figure for alpha diversity of tissues
p2 <- ggplot(Tissue21.meta) + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') +
  ggtitle('Alpha diversity of tissues') + annotate('text', x=4, y=36.5, label= '* p = 0.013 *')+
  annotate('segment', x=3.8, xend=4.2, y=35, yend=35)



### ###

FS2.meta$design <- as.factor(FS2.meta$design )


### ###

FS2 <- shared[meta$experiment == "FS2" & meta$day == 21,]
FS2.meta <- meta[meta$experiment == 'FS2'& meta$day == 21,]


colon.shared <- FS2[FS2.meta$tissue == 'colon', ]
colon.meta <- FS2.meta[FS2.meta$tissue == 'colon', ]
min(rowSums(colon.shared)) # 17890 seqs per sample for the colon tissue
colon.shared <- as.data.frame(rrarefy(colon.shared ,min(rowSums(colon.shared))))
#boxplot(diversity(colon.shared, index = 'invsimpson')~colon.meta$treatment)
#wilcox.test((diversity(colon.shared, index = 'invsimpson')~colon.meta$treatment))


# wrote these so I could use mothur's lefse implementation

#colon.shared <- cbind(label= '0.03', Group = rownames(colon.shared), numOtus=2069, colon.shared)
#write.table(colon.shared, 'colon.shared', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

cecum.shared <- FS2[FS2.meta$tissue == 'cecum', ]
cecum.meta <- FS2.meta[FS2.meta$tissue == 'cecum', ]
min(rowSums(cecum.shared)) # 31441 seqs per sample in the colon mucosa
cecum.shared <- as.data.frame(rrarefy(cecum.shared ,min(rowSums(cecum.shared))))
#boxplot(diversity(cecum.shared, index = 'invsimpson')~cecum.meta$treatment)
#wilcox.test((diversity(cecum.shared, index = 'invsimpson')~cecum.meta$treatment))
#cecum.shared <- cbind(label= '0.03', Group = rownames(cecum.shared), numOtus=2069, cecum.shared)
#write.table(cecum.shared, 'cecum.shared', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)


ileum.shared <- FS2[FS2.meta$tissue == 'ileum', ]
ileum.meta <- FS2.meta[FS2.meta$tissue == 'ileum', ]
ileum.meta <- ileum.meta[rowSums(ileum.shared) > 1000,]
ileum.shared <- ileum.shared[rowSums(ileum.shared) > 1000,]

min(rowSums(ileum.shared)) # only 674 per sequence
ileum.shared <- as.data.frame(rrarefy(ileum.shared ,min(rowSums(ileum.shared))))
#boxplot(diversity(ileum.shared, index = 'invsimpson')~ileum.meta$treatment)
#wilcox.test((diversity(ileum.shared, index = 'invsimpson')~ileum.meta$treatment))
#ileum.shared <- cbind(label= '0.03', Group = rownames(ileum.shared), numOtus=2069, ileum.shared)
#write.table(ileum.shared, 'ileum.shared', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

feces.shared <- FS2[FS2.meta$tissue == 'feces', ]
feces.meta <- FS2.meta[FS2.meta$tissue == 'feces', ]
feces.meta <- feces.meta[rowSums(feces.shared) > 2500,]
feces.shared <- feces.shared[rowSums(feces.shared) > 2500,]


feces.shared <- as.data.frame(rrarefy(feces.shared ,min(rowSums(feces.shared))))
#boxplot(diversity(feces.shared, index = 'invsimpson')~feces.meta$treatment)
#wilcox.test((diversity(feces.shared, index = 'invsimpson')~feces.meta$treatment))

#feces.shared <- cbind(label= '0.03', Group = rownames(feces.shared), numOtus=2069, feces.shared)
#write.table(feces.shared, 'feces.shared', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)


cec_cont_RNA.shared <- FS2[FS2.meta$tissue == 'cec_cont_RNA', ]
cec_cont_RNA.meta <- FS2.meta[FS2.meta$tissue == 'cec_cont_RNA', ]

#min(rowSums(cec_cont_RNA.shared)) # 21805 seqs per sample in cec cont RNA
cec_cont_RNA.shared <- as.data.frame(rrarefy(cec_cont_RNA.shared ,min(rowSums(cec_cont_RNA.shared))))

#boxplot(diversity(cec_cont_RNA.shared, index = 'invsimpson')~cec_cont_RNA.meta$treatment)
#wilcox.test((diversity(cec_cont_RNA.shared, index = 'invsimpson')~cec_cont_RNA.meta$treatment))



### START ORDS ####


#### Feces Ordination ####
# these ordinations use the OTU tables generated above, they only contain one type of tissue at d21 and are rarefied to the min number of reads in that tissue category
FS2.bray.feces <- vegdist(feces.shared, method = 'bray')
FS2.mds.feces <- metaMDS(FS2.bray.feces, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.feces$stress

FS2.nmds.feces <-as.data.frame(FS2.mds.feces$points)
FS2.nmds.feces$group <- rownames(FS2.nmds.feces)

FS2.metanmds.feces <- merge(meta, FS2.nmds.feces, by = 'group')
FS2.metanmds.feces$design <- factor(FS2.metanmds.feces$design)
FS2.metanmds.feces$treatment <- factor(FS2.metanmds.feces$treatment)



FS2.metanmds.feces$treatmentXday <- paste(FS2.metanmds.feces$treatment, FS2.metanmds.feces$day, sep = ' day ')
FS2.metanmds.feces$treatmentXday <- factor(FS2.metanmds.feces$treatmentXday)


ord <- ordiellipse(FS2.mds.feces,FS2.metanmds.feces$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.feces <- aggregate(FS2.metanmds.feces[,10:11], list(group=FS2.metanmds.feces$treatment), mean)

df_ell <- data.frame()
for (d in levels(FS2.metanmds.feces$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds.feces[FS2.metanmds.feces$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

FS2.metanmds.feces$centroidX <- NA
FS2.metanmds.feces$centroidY <- NA


FS2.metanmds.feces[FS2.metanmds.feces$treatment == 'control',]$centroidX <- NMDS.mean.feces$MDS1[1]
FS2.metanmds.feces[FS2.metanmds.feces$treatment == 'RPS',]$centroidX <- NMDS.mean.feces$MDS1[2]

FS2.metanmds.feces[FS2.metanmds.feces$treatment == 'control',]$centroidY <- NMDS.mean.feces$MDS2[1]
FS2.metanmds.feces[FS2.metanmds.feces$treatment == 'RPS',]$centroidY <- NMDS.mean.feces$MDS2[2]


### Statistical Tests ###

adonis.feces <- adonis(feces.shared~feces.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.feces, group = FS2.metanmds.feces$treatment)
permutest(BDISP, permutations = 9999)

### ###



p.feces <- ggplot(FS2.metanmds.feces, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.feces, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Fecal beta-diversity: NMDS ordination using Bray-Curtis distances',
          subtitle = 'Group similarity (PERMANOVA): p = 1e-05, F = 5.66, R2 = 0.18\nGroup dispersion (PERMDISP2): p = 0.0015, F = 13.12') + 
  geom_segment(data = FS2.metanmds.feces,
               aes(x=FS2.metanmds.feces$MDS1,
                   xend=FS2.metanmds.feces$centroidX,
                   y=FS2.metanmds.feces$MDS2,
                   yend=FS2.metanmds.feces$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") + 
  labs(caption = 'Ordination stress = 0.14')

#p.feces

####### ileum ord #####

FS2.bray.ileum <- vegdist(ileum.shared, method = 'bray')
FS2.mds.ileum <- metaMDS(FS2.bray.ileum, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.ileum$stress

FS2.nmds.ileum <-as.data.frame(FS2.mds.ileum$points)
FS2.nmds.ileum$group <- rownames(FS2.nmds.ileum)

FS2.metanmds.ileum <- merge(meta, FS2.nmds.ileum, by = 'group')
FS2.metanmds.ileum$design <- factor(FS2.metanmds.ileum$design)
FS2.metanmds.ileum$treatment <- factor(FS2.metanmds.ileum$treatment)



FS2.metanmds.ileum$treatmentXday <- paste(FS2.metanmds.ileum$treatment, FS2.metanmds.ileum$day, sep = ' day ')
FS2.metanmds.ileum$treatmentXday <- factor(FS2.metanmds.ileum$treatmentXday)


ord <- ordiellipse(FS2.mds.ileum,FS2.metanmds.ileum$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.ileum <- aggregate(FS2.metanmds.ileum[,10:11], list(group=FS2.metanmds.ileum$treatment), mean)


df_ell <- data.frame()
for (d in levels(FS2.metanmds.ileum$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds.ileum[FS2.metanmds.ileum$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

FS2.metanmds.ileum$centroidX <- NA
FS2.metanmds.ileum$centroidY <- NA


FS2.metanmds.ileum[FS2.metanmds.ileum$treatment == 'control',]$centroidX <- NMDS.mean.ileum$MDS1[1]
FS2.metanmds.ileum[FS2.metanmds.ileum$treatment == 'RPS',]$centroidX <- NMDS.mean.ileum$MDS1[2]

FS2.metanmds.ileum[FS2.metanmds.ileum$treatment == 'control',]$centroidY <- NMDS.mean.ileum$MDS2[1]
FS2.metanmds.ileum[FS2.metanmds.ileum$treatment == 'RPS',]$centroidY <- NMDS.mean.ileum$MDS2[2]

### Statistical Tests ###

BDISP <- betadisper(FS2.bray.ileum, group = FS2.metanmds.ileum$treatment)
permutest(BDISP, permutations = 99999)
adonis(ileum.shared~ileum.meta$treatment, permutations = 99999)

### ###

p.ileum <- ggplot(FS2.metanmds.ileum, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.ileum, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Ileum mucosa bacterial beta-diversity: NMDS ordination using Bray-Curtis distances',
          subtitle = 'Group similarity (PERMANOVA) p = 0.002, F = 2.62, R2 = 0.19\nGroup dispersion (PERMDISP2): p = 0.55, F = 0.374') + 
  geom_segment(data = FS2.metanmds.ileum,
               aes(x=FS2.metanmds.ileum$MDS1,
                   xend=FS2.metanmds.ileum$centroidX,
                   y=FS2.metanmds.ileum$MDS2,
                   yend=FS2.metanmds.ileum$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") +
  labs(caption = 'Ordination stress = 0.12')


#p.ileum


###### cecum ord #######

FS2.bray.cecum <- vegdist(cecum.shared, method = 'bray')
FS2.mds.cecum <- metaMDS(FS2.bray.cecum, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.cecum$stress

FS2.nmds.cecum <-as.data.frame(FS2.mds.cecum$points)
FS2.nmds.cecum$group <- rownames(FS2.nmds.cecum)

FS2.metanmds.cecum <- merge(meta, FS2.nmds.cecum, by = 'group')
FS2.metanmds.cecum$design <- factor(FS2.metanmds.cecum$design)
FS2.metanmds.cecum$treatment <- factor(FS2.metanmds.cecum$treatment)



FS2.metanmds.cecum$treatmentXday <- paste(FS2.metanmds.cecum$treatment, FS2.metanmds.cecum$day, sep = ' day ')
FS2.metanmds.cecum$treatmentXday <- factor(FS2.metanmds.cecum$treatmentXday)


ord <- ordiellipse(FS2.mds.cecum,FS2.metanmds.cecum$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.cecum <- aggregate(FS2.metanmds.cecum[,10:11], list(group=FS2.metanmds.cecum$treatment), mean)

df_ell <- data.frame()
for (d in levels(FS2.metanmds.cecum$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds.cecum[FS2.metanmds.cecum$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

FS2.metanmds.cecum$centroidX <- NA
FS2.metanmds.cecum$centroidY <- NA


FS2.metanmds.cecum[FS2.metanmds.cecum$treatment == 'control',]$centroidX <- NMDS.mean.cecum$MDS1[1]
FS2.metanmds.cecum[FS2.metanmds.cecum$treatment == 'RPS',]$centroidX <- NMDS.mean.cecum$MDS1[2]

FS2.metanmds.cecum[FS2.metanmds.cecum$treatment == 'control',]$centroidY <- NMDS.mean.cecum$MDS2[1]
FS2.metanmds.cecum[FS2.metanmds.cecum$treatment == 'RPS',]$centroidY <- NMDS.mean.cecum$MDS2[2]


### Statistical Tests ###

adonis(cecum.shared~cecum.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.cecum, group = FS2.metanmds.cecum$treatment)
permutest(BDISP, permutations = 9999)

### ###

p.cecum <- ggplot(FS2.metanmds.cecum, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.cecum, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Cecal mucosa bacterial beta-diversity: NMDS ordination using Bray-Curtis distances',
          subtitle = 'Group similarity (PERMANOVA): p = 0.0077, F = 2.67, R2 = 0.18\nGroup dispersion (PERMDISP2): p = 0.003, F = 12.93') + 
  geom_segment(data = FS2.metanmds.cecum,
               aes(x=FS2.metanmds.cecum$MDS1,
                   xend=FS2.metanmds.cecum$centroidX,
                   y=FS2.metanmds.cecum$MDS2,
                   yend=FS2.metanmds.cecum$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") +
  labs(caption = 'Ordination stress = 0.12')
#geom_text_repel(data=data.frame(),aes(x = Inf, y = -Inf), label=paste('stress = ', round(FS2.mds.cecum$stress, 3), sep = ''), segment.size = 0)
#p.cecum

######## colon ord #####

FS2.bray.colon <- vegdist(colon.shared, method = 'bray')
FS2.mds.colon <- metaMDS(FS2.bray.colon, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.colon$stress

FS2.nmds.colon <-as.data.frame(FS2.mds.colon$points)
FS2.nmds.colon$group <- rownames(FS2.nmds.colon)

FS2.metanmds.colon <- merge(meta, FS2.nmds.colon, by = 'group')
FS2.metanmds.colon$design <- factor(FS2.metanmds.colon$design)
FS2.metanmds.colon$treatment <- factor(FS2.metanmds.colon$treatment)



FS2.metanmds.colon$treatmentXday <- paste(FS2.metanmds.colon$treatment, FS2.metanmds.colon$day, sep = ' day ')
FS2.metanmds.colon$treatmentXday <- factor(FS2.metanmds.colon$treatmentXday)


ord <- ordiellipse(FS2.mds.colon,FS2.metanmds.colon$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.colon <- aggregate(FS2.metanmds.colon[,10:11], list(group=FS2.metanmds.colon$treatment), mean)


df_ell <- data.frame()
for (d in levels(FS2.metanmds.colon$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds.colon[FS2.metanmds.colon$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

FS2.metanmds.colon$centroidX <- NA
FS2.metanmds.colon$centroidY <- NA

FS2.metanmds.colon[FS2.metanmds.colon$treatment == 'control',]$centroidX <- NMDS.mean.colon$MDS1[1]
FS2.metanmds.colon[FS2.metanmds.colon$treatment == 'RPS',]$centroidX <- NMDS.mean.colon$MDS1[2]

FS2.metanmds.colon[FS2.metanmds.colon$treatment == 'control',]$centroidY <- NMDS.mean.colon$MDS2[1]
FS2.metanmds.colon[FS2.metanmds.colon$treatment == 'RPS',]$centroidY <- NMDS.mean.colon$MDS2[2]

### Statistical tests ###

adonis(colon.shared~colon.meta$treatment, permutations = 99999)
BDISP <- betadisper(FS2.bray.colon, group = FS2.metanmds.colon$treatment)
permutest(BDISP, permutations = 99999)

### ###

p.colon <- ggplot(FS2.metanmds.colon, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.colon, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Colonic mucosa bacterial beta-diversity: NMDS ordination using Bray-Curtis distances', 
          subtitle = 'Group similarity (PERMANOVA): p = 0.171, F = 1.36, R2 = 0.10\n Group dispesion (PERMDISP2): p = 0.133, F = 2.65') + 
  geom_segment(data = FS2.metanmds.colon,
               aes(x=FS2.metanmds.colon$MDS1,
                   xend=FS2.metanmds.colon$centroidX,
                   y=FS2.metanmds.colon$MDS2,
                   yend=FS2.metanmds.colon$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") +
  labs(caption = 'Ordination stress = 0.112')

#p.colon

######## cec_cont_RNA ord #####

FS2.bray.cec_cont_RNA <- vegdist(cec_cont_RNA.shared, method = 'bray')
FS2.mds.cec_cont_RNA <- metaMDS(FS2.bray.cec_cont_RNA, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.cec_cont_RNA$stress

FS2.nmds.cec_cont_RNA <-as.data.frame(FS2.mds.cec_cont_RNA$points)
FS2.nmds.cec_cont_RNA$group <- rownames(FS2.nmds.cec_cont_RNA)

FS2.metanmds.cec_cont_RNA <- merge(meta, FS2.nmds.cec_cont_RNA, by = 'group')
FS2.metanmds.cec_cont_RNA$design <- factor(FS2.metanmds.cec_cont_RNA$design)
FS2.metanmds.cec_cont_RNA$treatment <- factor(FS2.metanmds.cec_cont_RNA$treatment)



#FS2.metanmds.cec_cont_RNA$treatmentXday <- paste(FS2.metanmds.cec_cont_RNA$treatment, FS2.metanmds.cec_cont_RNA$day, sep = ' day ')
#FS2.metanmds.cec_cont_RNA$treatmentXday <- factor(FS2.metanmds.cec_cont_RNA$treatmentXday)


ord <- ordiellipse(FS2.mds.cec_cont_RNA,FS2.metanmds.cec_cont_RNA$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.cec_cont_RNA <- aggregate(FS2.metanmds.cec_cont_RNA[,10:11], list(group=FS2.metanmds.cec_cont_RNA$treatment), mean)

df_ell <- data.frame()
for (d in levels(FS2.metanmds.cec_cont_RNA$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds.cec_cont_RNA[FS2.metanmds.cec_cont_RNA$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

FS2.metanmds.cec_cont_RNA$centroidX <- NA
FS2.metanmds.cec_cont_RNA$centroidY <- NA


FS2.metanmds.cec_cont_RNA[FS2.metanmds.cec_cont_RNA$treatment == 'control',]$centroidX <- NMDS.mean.cec_cont_RNA$MDS1[1]
FS2.metanmds.cec_cont_RNA[FS2.metanmds.cec_cont_RNA$treatment == 'RPS',]$centroidX <- NMDS.mean.cec_cont_RNA$MDS1[2]

FS2.metanmds.cec_cont_RNA[FS2.metanmds.cec_cont_RNA$treatment == 'control',]$centroidY <- NMDS.mean.cec_cont_RNA$MDS2[1]
FS2.metanmds.cec_cont_RNA[FS2.metanmds.cec_cont_RNA$treatment == 'RPS',]$centroidY <- NMDS.mean.cec_cont_RNA$MDS2[2]


### Statistical Tests ###

adonis(cec_cont_RNA.shared~cec_cont_RNA.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.cec_cont_RNA, group = FS2.metanmds.cec_cont_RNA$treatment)
permutest(BDISP, permutations = 9999)
### ###



p.cec_cont_RNA <- ggplot(FS2.metanmds.cec_cont_RNA, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.cec_cont_RNA, aes(color=treatment), size=2) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Cecal contents beta-diversity: NMDS ordination using Bray-Curtis distances',
          subtitle = 'Group similarity (PERMANOVA): p = 0.061, F = 1.88, R2 = 0.14\nGroup dispersion (PERMDISP2): p = 0.909, F = 0.0156') + 
  geom_segment(data = FS2.metanmds.cec_cont_RNA,
               aes(x=FS2.metanmds.cec_cont_RNA$MDS1,
                   xend=FS2.metanmds.cec_cont_RNA$centroidX,
                   y=FS2.metanmds.cec_cont_RNA$MDS2,
                   yend=FS2.metanmds.cec_cont_RNA$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") + 
  labs(caption = 'Ordination stress = 0.11')

#p.cec_cont_RNA

######################### Deseq2 stuff ######################################################################





otu <- import_mothur(mothur_shared_file = 'V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'V4.final.taxonomy')
meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]



FS2 <- phyloseq(otu, taxo)
FS2 <- merge_phyloseq(FS2, phy_meta)                       # combines the metadata with this phyloseq object
colnames(tax_table(FS2)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS2 <- subset_samples(FS2, experiment == 'FS2')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.
FS2 <- prune_taxa(taxa_sums(FS2) > 10, FS2)        # removes OTUs that occur less than 10 times globally

FS2.genus <- tax_glom(FS2, taxrank = "Genus")

######## D0 #########


FS2.D0 <- subset_samples(FS2.genus, day == 0)

sample_sums(FS2.D0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)

rowSums(FS2.D0@otu_table)

FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)

FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")


res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0 = res.D0[which(res.D0$padj < .05), ]

# No sigdiff genera at D0
# sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS2.D0)[rownames(sigtab.D0), ], "matrix"))
# format(sigtab.D0$padj, scientific = TRUE)
# sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
# sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "RPS", "Control")
# 
# 
# deseq.D0 <- ggplot(sigtab.D0, aes(x=reorder(rownames(sigtab.D0), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0), y=log2FoldChange+.6, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
#                                              axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#                                              axis.title.x=element_text(size = 10),
#                                              axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
# deseq.D0


########## D12 ###############
##########     ###############

FS2.D12 <- subset_samples(FS2.genus, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")


res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D12 = res.D12[which(res.D12$padj < .05), ]
sigtab.D12 = cbind(as(sigtab.D12, "data.frame"), as(tax_table(FS2.D12)[rownames(sigtab.D12), ], "matrix"))
format(sigtab.D12$padj, scientific = TRUE)
sigtab.D12$newp <- format(round(sigtab.D12$padj, digits = 3), scientific = TRUE)
sigtab.D12$Treatment <- ifelse(sigtab.D12$log2FoldChange >=0, "RPS", "Control")


deseq.D12 <- ggplot(sigtab.D12, aes(x=reorder(rownames(sigtab.D12), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D12), y=log2FoldChange+.6, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
#deseq.D12

# 4 genera enriched in RPS pigs, none in control, not a good figure....

################### D15 ####################

#FS2.D15 <- subset_samples(FS2, day == 15)
FS2.D15 <- subset_samples(FS2.genus, day == 15)
sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)


FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")


res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D15 = res.D15[which(res.D15$padj < .05), ]
# sigtab.D15 = cbind(as(sigtab.D15, "data.frame"), as(tax_table(FS2.D15)[rownames(sigtab.D15), ], "matrix"))
# format(sigtab.D15$padj, scientific = TRUE)
# sigtab.D15$newp <- format(round(sigtab.D15$padj, digits = 3), scientific = TRUE)
# sigtab.D15$Treatment <- ifelse(sigtab.D15$log2FoldChange >=0, "RPS", "Control")
# 
# 
# deseq.D15 <- ggplot(sigtab.D15, aes(x=reorder(rownames(sigtab.D15), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D15), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
#                                              axis.text.y=element_blank(), 
#                                              axis.title.x=element_text(size = 10),
#                                              axis.title.y=element_blank(),
#                                              axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: feces', subtitle = 'Day 15')+ coord_flip()
# deseq.D15


# NO DIFF ABUND GENERA AT D15
# There are 3 diff abund OTUs though.....
################### D19 ####################


FS2.D19 <- subset_samples(FS2.genus, day == 19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")



res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D19 = res.D19[which(res.D19$padj < .05), ]
sigtab.D19 = cbind(as(sigtab.D19, "data.frame"), as(tax_table(FS2.D19)[rownames(sigtab.D19), ], "matrix"))
format(sigtab.D19$padj, scientific = TRUE)
sigtab.D19$newp <- format(round(sigtab.D19$padj, digits = 3), scientific = TRUE)
sigtab.D19$Treatment <- ifelse(sigtab.D19$log2FoldChange >=0, "RPS", "Control")


#tax_tab <- as.data.frame(FS2@tax_table)


#levels(tax_tab$Phy)
deseq.D19 <- ggplot(sigtab.D19, aes(x=reorder(rownames(sigtab.D19), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D19), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: feces', subtitle = 'Day 19')+ coord_flip()
#deseq.D19



# not a bad fig, but maybe better for supplement, and use ordinations to show timepoints prior to D21
#### D21 should be here ####


FS2.D21 <- subset_samples(FS2.genus, day == 21 & tissue == 'feces')

########## CHANGED THIS TO LUMP ALL FECAL TIMEPOINTS IS THIS APPROPRIATE?????  ####################


FS2.D21 <- subset_samples(FS2.genus, day %in% c(21) & tissue == 'feces')




FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)
FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")



res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21 = res.D21[which(res.D21$padj < .05), ]
sigtab.D21 = cbind(as(sigtab.D21, "data.frame"), as(tax_table(FS2.D21)[rownames(sigtab.D21), ], "matrix"))
format(sigtab.D21$padj, scientific = TRUE)
sigtab.D21$newp <- format(round(sigtab.D21$padj, digits = 3), scientific = TRUE)
sigtab.D21$Treatment <- ifelse(sigtab.D21$log2FoldChange >=0, "RPS", "Control")


deseq.D21 <- ggplot(sigtab.D21, aes(x=reorder(rownames(sigtab.D21), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: feces', subtitle = 'Day 21')+ coord_flip()
#deseq.D21





############ This is where I break down the comparisons before running DeSeq2 ###################
############ Reading the DeSeq2 manual, I learned this may be necessary when there are strong dispersion effects 
############ in some groups and not others, lumping everything can lear to an overestimation of the dispersion in some way.
##### This is all D21 tissue stuff #####

FS2.cecum <- subset_samples(FS2.genus, tissue == 'cecum' & day ==21)
FS2.ileum <- subset_samples(FS2.genus, tissue == 'ileum'& day ==21)
FS2.colon <- subset_samples(FS2.genus, tissue == 'colon'& day ==21)
FS2.cec_cont_RNA <- subset_samples(FS2.genus, tissue == 'cec_cont_RNA'& day ==21)

FS2.cecum.De <- phyloseq_to_deseq2(FS2.cecum, ~ treatment)
FS2.colon.De <- phyloseq_to_deseq2(FS2.colon, ~ treatment)
FS2.ileum.De <- phyloseq_to_deseq2(FS2.ileum, ~ treatment)
FS2.cec_cont_RNA.De <- phyloseq_to_deseq2(FS2.cec_cont_RNA, ~ treatment)




FS2.cecum.De <- DESeq(FS2.cecum.De, test = "Wald", fitType = "parametric")
FS2.colon.De <- DESeq(FS2.colon.De, test = "Wald", fitType = "parametric")
FS2.ileum.De <- DESeq(FS2.ileum.De, test = "Wald", fitType = "parametric")
FS2.cec_cont_RNA.De <- DESeq(FS2.cec_cont_RNA.De, test = "Wald", fitType = "parametric")




res.cecum = results(FS2.cecum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.colon <- results(FS2.colon.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.ileum <- results(FS2.ileum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.cec_cont_RNA <- results(FS2.cec_cont_RNA.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')

alpha = 0.05

################ ##################

sigtab.ileum = res.ileum[which(res.ileum$padj < alpha), ]
sigtab.ileum = cbind(as(sigtab.ileum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.ileum), ], "matrix"))
format(sigtab.ileum$padj, scientific = TRUE)
sigtab.ileum$newp <- format(round(sigtab.ileum$padj, digits = 3), scientific = TRUE)
sigtab.ileum$Treatment <- ifelse(sigtab.ileum$log2FoldChange >=0, "RPS", "Control")

deseq.ileum <- ggplot(sigtab.ileum, aes(x=reorder(rownames(sigtab.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.ileum), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank()) +
  ggtitle('Differentially abundant genera: ileal mucosa')+ coord_flip()
# deseq.ileum

sigtab.cecum = res.cecum[which(res.cecum$padj < alpha), ]
sigtab.cecum = cbind(as(sigtab.cecum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cecum), ], "matrix"))
format(sigtab.cecum$padj, scientific = TRUE)
sigtab.cecum$newp <- format(round(sigtab.cecum$padj, digits = 3), scientific = TRUE)
sigtab.cecum$Treatment <- ifelse(sigtab.cecum$log2FoldChange >=0, "RPS", "Control")

deseq.cecum <- ggplot(sigtab.cecum, aes(x=reorder(rownames(sigtab.cecum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.cecum), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: cecal mucosa')+ coord_flip()
# deseq.cecum


sigtab.colon = res.colon[which(res.colon$padj < alpha), ]
sigtab.colon = cbind(as(sigtab.colon, "data.frame"), as(tax_table(FS2.genus)[rownames(sigtab.colon), ], "matrix"))
format(sigtab.colon$padj, scientific = TRUE)
sigtab.colon$newp <- format(round(sigtab.colon$padj, digits = 3), scientific = TRUE)
sigtab.colon$Treatment <- ifelse(sigtab.colon$log2FoldChange >=0, "RPS", "Control")

deseq.colon <- ggplot(sigtab.colon, aes(x=reorder(rownames(sigtab.colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.colon), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: colonic mucosa')+ coord_flip()
# deseq.colon


sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < alpha), ]
sigtab.cec_cont_RNA = cbind(as(sigtab.cec_cont_RNA, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cec_cont_RNA), ], "matrix"))
format(sigtab.cec_cont_RNA$padj, scientific = TRUE)
sigtab.cec_cont_RNA$newp <- format(round(sigtab.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.cec_cont_RNA$Treatment <- ifelse(sigtab.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

deseq.cec_cont_RNA <- ggplot(sigtab.cec_cont_RNA, aes(x=reorder(rownames(sigtab.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.cec_cont_RNA), y=0, label = Genus), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: cecal contents (RNA)')+ coord_flip()
# deseq.cec_cont_RNA


########################## Now doing everything over just at OTU level this time   ######################
# this is for supplementary table stuff #


otu <- import_mothur(mothur_shared_file = 'V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'V4.final.taxonomy')
meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]



FS2 <- phyloseq(otu, taxo)
FS2 <- merge_phyloseq(FS2, phy_meta)                       # combines the metadata with this phyloseq object
colnames(tax_table(FS2)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS2 <- subset_samples(FS2, experiment == 'FS2')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.
FS2 <- prune_taxa(taxa_sums(FS2) > 10, FS2)        # removes OTUs that occur less than 10 times globally


alpha = 0.3 # Changing the alpha cutoff to 0.2 so that each timepoint can have a table

FS2.D0 <- subset_samples(FS2, day == 0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)
FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)
FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")

res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0 = res.D0[which(res.D0$padj < alpha), ]

sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS2.D0)[rownames(sigtab.D0), ], "matrix"))
format(sigtab.D0$padj, scientific = TRUE)
sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "RPS", "Control")

write.table(sigtab.D0, file = './OTUs/D0_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

########## D12 ###############
##########     ###############

FS2.D12 <- subset_samples(FS2, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")


res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D12 = res.D12[which(res.D12$padj < alpha), ]
sigtab.D12 = cbind(as(sigtab.D12, "data.frame"), as(tax_table(FS2.D12)[rownames(sigtab.D12), ], "matrix"))
format(sigtab.D12$padj, scientific = TRUE)
sigtab.D12$newp <- format(round(sigtab.D12$padj, digits = 3), scientific = TRUE)
sigtab.D12$Treatment <- ifelse(sigtab.D12$log2FoldChange >=0, "RPS", "Control")


write.table(sigtab.D12, file = './OTUs/D12_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')


################### D15 ####################

FS2.D15 <- subset_samples(FS2, day == 15)
sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)


FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")


res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D15 = res.D15[which(res.D15$padj < alpha), ]
sigtab.D15 = cbind(as(sigtab.D15, "data.frame"), as(tax_table(FS2.D15)[rownames(sigtab.D15), ], "matrix"))
format(sigtab.D15$padj, scientific = TRUE)
sigtab.D15$newp <- format(round(sigtab.D15$padj, digits = 3), scientific = TRUE)
sigtab.D15$Treatment <- ifelse(sigtab.D15$log2FoldChange >=0, "RPS", "Control")


write.table(sigtab.D15, file = './OTUs/D15_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

################### D19 ####################


FS2.D19 <- subset_samples(FS2, day == 19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")



res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D19 = res.D19[which(res.D19$padj < alpha), ]
sigtab.D19 = cbind(as(sigtab.D19, "data.frame"), as(tax_table(FS2.D19)[rownames(sigtab.D19), ], "matrix"))
format(sigtab.D19$padj, scientific = TRUE)
sigtab.D19$newp <- format(round(sigtab.D19$padj, digits = 3), scientific = TRUE)
sigtab.D19$Treatment <- ifelse(sigtab.D19$log2FoldChange >=0, "RPS", "Control")


write.table(sigtab.D19, file = './OTUs/D19_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')


#### D21 #######


FS2.D21 <- subset_samples(FS2, day == 21 & tissue == 'feces')

FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)
FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")



res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21 = res.D21[which(res.D21$padj < 0.05), ]
sigtab.D21 = cbind(as(sigtab.D21, "data.frame"), as(tax_table(FS2.D21)[rownames(sigtab.D21), ], "matrix"))
format(sigtab.D21$padj, scientific = TRUE)
sigtab.D21$newp <- format(round(sigtab.D21$padj, digits = 3), scientific = TRUE)
sigtab.D21$Treatment <- ifelse(sigtab.D21$log2FoldChange >=0, "RPS", "Control")


write.table(sigtab.D21, file = './OTUs/D21_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')


##### This is all D21 tissue stuff #####

FS2.cecum <- subset_samples(FS2, tissue == 'cecum' & day ==21)
FS2.ileum <- subset_samples(FS2, tissue == 'ileum'& day ==21)
FS2.colon <- subset_samples(FS2, tissue == 'colon'& day ==21)
FS2.cec_cont_RNA <- subset_samples(FS2, tissue == 'cec_cont_RNA'& day ==21)

FS2.cecum.De <- phyloseq_to_deseq2(FS2.cecum, ~ treatment)
FS2.colon.De <- phyloseq_to_deseq2(FS2.colon, ~ treatment)
FS2.ileum.De <- phyloseq_to_deseq2(FS2.ileum, ~ treatment)
FS2.cec_cont_RNA.De <- phyloseq_to_deseq2(FS2.cec_cont_RNA, ~ treatment)

FS2.cecum.De <- DESeq(FS2.cecum.De, test = "Wald", fitType = "parametric")
FS2.colon.De <- DESeq(FS2.colon.De, test = "Wald", fitType = "parametric")
FS2.ileum.De <- DESeq(FS2.ileum.De, test = "Wald", fitType = "parametric")
FS2.cec_cont_RNA.De <- DESeq(FS2.cec_cont_RNA.De, test = "Wald", fitType = "parametric")

res.cecum = results(FS2.cecum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.colon <- results(FS2.colon.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.ileum <- results(FS2.ileum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.cec_cont_RNA <- results(FS2.cec_cont_RNA.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')

sigtab.ileum = res.ileum[which(res.ileum$padj < alpha), ]
sigtab.ileum = cbind(as(sigtab.ileum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.ileum), ], "matrix"))
format(sigtab.ileum$padj, scientific = TRUE)
sigtab.ileum$newp <- format(round(sigtab.ileum$padj, digits = 3), scientific = TRUE)
sigtab.ileum$Treatment <- ifelse(sigtab.ileum$log2FoldChange >=0, "RPS", "Control")

sigtab.cecum = res.cecum[which(res.cecum$padj < alpha), ]
sigtab.cecum = cbind(as(sigtab.cecum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cecum), ], "matrix"))
format(sigtab.cecum$padj, scientific = TRUE)
sigtab.cecum$newp <- format(round(sigtab.cecum$padj, digits = 3), scientific = TRUE)
sigtab.cecum$Treatment <- ifelse(sigtab.cecum$log2FoldChange >=0, "RPS", "Control")

sigtab.colon = res.colon[which(res.colon$padj < alpha), ]
sigtab.colon = cbind(as(sigtab.colon, "data.frame"), as(tax_table(FS2)[rownames(sigtab.colon), ], "matrix"))
format(sigtab.colon$padj, scientific = TRUE)
sigtab.colon$newp <- format(round(sigtab.colon$padj, digits = 3), scientific = TRUE)
sigtab.colon$Treatment <- ifelse(sigtab.colon$log2FoldChange >=0, "RPS", "Control")

sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < alpha), ]
sigtab.cec_cont_RNA = cbind(as(sigtab.cec_cont_RNA, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cec_cont_RNA), ], "matrix"))
format(sigtab.cec_cont_RNA$padj, scientific = TRUE)
sigtab.cec_cont_RNA$newp <- format(round(sigtab.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.cec_cont_RNA$Treatment <- ifelse(sigtab.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

write.table(sigtab.cec_cont_RNA, file = './OTUs/cec_cont_RNA_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
write.table(sigtab.colon, file = './OTUs/colon_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
write.table(sigtab.cecum, file = './OTUs/cecum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
write.table(sigtab.ileum, file = './OTUs/ileum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')




################ end OTU level stuff  ##############


##### most stuff below here is everything rarrefied to 4200 sequences? ####



FS2.4200 <- read.table("V4.FS2.4200.shared", header = TRUE)
rownames(FS2.4200) <- FS2.4200$Group
meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)

meta <- meta[meta$group %in% rownames(FS2.4200),]
meta <- meta[meta$design != "FS2_cec_cont_DNA_21_RPS",]

FS2.4200 <- FS2.4200[rownames(FS2.4200) %in% meta$group,]



FS2.4200 <- FS2.4200[,-c(1,2,3)]


#####  #####


############### REPLACE WITH EXTRACT MOTHUR TAX ##################

tax <- read.table('V4.final.taxonomy', header = TRUE, stringsAsFactors = FALSE)
tna <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
tax[tna] <- NA
rowcount <- 1

for (i in strsplit(tax$Taxonomy, split = ';')){
  levelcount <- 1
  for (x in i){
    tax[rowcount,3+levelcount] <- x
    levelcount <- levelcount + 1
    #print(x)
  }
  rowcount <- rowcount+1
}



timepo <- FS2.4200[,-1]/4200

timepo$group <- rownames(timepo)
timepo <- merge(meta, timepo, by = 'group')
timepo$day <- as.numeric(timepo$day)


timepo2 <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(tissue == 'feces', day <  22)
timepo2$otu2 <- factor(timepo2$otu)

tiss <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(day <  22)



############# ELIMINATE INDIVIDUAL OTU OBJECTS #############

OTU87 <- tiss %>% filter(otu == 'Otu00087' & day ==21)
p.otu87.tiss <- ggplot(OTU87, aes(x=tissue, y=value))  +
  geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  ggtitle('Relative proportion of OTU87 in each tissue a Day 21') +
  scale_fill_brewer(palette = 'Dark2')

OTU15 <- tiss %>% filter(otu == 'Otu00015' & day ==21)
p.otu15.tiss <- ggplot(OTU15, aes(x=tissue, y=value))  +  
  geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  ggtitle('Relative proportion of OTU15 in each tissue a Day 21') +
  scale_fill_brewer(palette = 'Dark2') + geom_jitter() 


OTU15 <- tiss %>% filter(otu == 'Otu00015' & day ==21)
ggplot(OTU15, aes(x=tissue, y=value))  +  
  geom_jitter(aes(color=treatment)) +
  ggtitle('Relative proportion of OTU15 in each tissue a Day 21') +
  scale_color_brewer(palette = 'Dark2') + geom_text(aes(label = pig_num))





# OTU40 <- tiss %>% filter(otu == 'Otu00040' & day ==21)
# ggplot(OTU40, aes(x=tissue, y=value))  +
#   geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
#   ggtitle('Relative proportion of OTU40 in each tissue a Day 21') +
#   scale_fill_brewer(palette = 'Dark2')

# OTU225 <- tiss %>% filter(otu == 'Otu00225' & day ==21)
# ggplot(OTU225, aes(x=tissue, y=value))  +
#   geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
#   ggtitle('Relative proportion of OTU225 in each tissue a Day 21') +
#   scale_fill_brewer(palette = 'Dark2')

# OTU325 <- tiss %>% filter(otu == 'Otu00325' & day ==21)
# ggplot(OTU325, aes(x=tissue, y=value))  +
#   geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
#   ggtitle('Relative proportion of OTU325 in each tissue a Day 21') +
#   scale_fill_brewer(palette = 'Dark2')


# C_SS_1 <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(39,85,223,322)])
# ggplot(C_SS_1, aes(x=day, y=value))  +
#   #geom_smooth(aes(colour=treatment, linetype=treatment), method = 'loess', alpha=.1) +
#   geom_vline(xintercept = 12, colour='black')+
#   geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + 
#   facet_wrap(~otu, scales = 'free')

######  OTU 87 OVER TIME  ###########
timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(85)]) %>% 
  ggplot(aes(x=day, y=value))  +
  #geom_smooth(aes(colour=treatment, linetype=treatment), method = 'loess', alpha=.1) +
  geom_vline(xintercept = 12, colour='black')+
  geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + ggtitle('Relative proportion of OTU87 in feces over time') + scale_fill_brewer(palette = 'Dark2')


timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(39)]) %>% 
  ggplot(aes(x=day, y=value))  +
  #geom_smooth(aes(colour=treatment, linetype=treatment), method = 'loess', alpha=.1) +
  geom_vline(xintercept = 12, colour='black')+
  geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + ggtitle('Relative proportion of OTU40 in feces over time') + scale_fill_brewer(palette = 'Dark2')

# OTU255 <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(223)])
# ggplot(OTU255, aes(x=day, y=value))  +
#   #geom_smooth(aes(colour=treatment, linetype=treatment), method = 'loess', alpha=.1) +
#   geom_vline(xintercept = 12, colour='black')+
#   geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + ggtitle('Relative proportion of OTU255 in feces over time') + scale_fill_brewer(palette = 'Dark2')

# OTU325 <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(322)])
# ggplot(OTU325, aes(x=day, y=value))  +
#   #geom_smooth(aes(colour=treatment, linetype=treatment), method = 'loess', alpha=.1) +
#   geom_vline(xintercept = 12, colour='black')+
#   geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + ggtitle('Relative proportion of OTU325 in feces over time') + scale_fill_brewer(palette = 'Dark2')
# 

#####  #####

FS2.bray2 <- vegdist(FS2.4200, method = 'bray')

### dispersion stuff, this has its own metathing so hopefully later steps still work ###

dispers <- betadisper(FS2.bray2, group = meta$design)
pdispers <- permutest(dispers, pairwise = TRUE)
pdispers$pairwise
dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group
metadisp <- merge(meta, dispersdf, by = 'group')

dispgroups <- summarize(group_by(metadisp, design), average_dist=mean(dispers.distances))

dispgroups <- unique(inner_join(dispgroups, meta[,-c(1,2,5,7,9)]))

metadisp %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) + 
  geom_boxplot() + scale_fill_brewer(palette = 'Dark2') + ylim(c(.15,.7)) + 
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")


# maybe this is useful?  polish a little bit

fecaldisptime <-  dispgroups %>% filter(tissue == 'feces' & day < 22)

ggplot(fecaldisptime, aes(x=day, y=average_dist, color = treatment)) + 
  geom_vline(xintercept = 12)+ geom_point(size = 2.5) +  geom_line() + ylim(c(.2,.5))+
  ggtitle('Community Variability (Dispersion)',
          subtitle = "Vegan's betadisper(): how much variability is there in a group's community structure?") + 
  ylab("Average distance to the group median") + scale_color_brewer(palette = "Dark2")


# I think I should include this  figure in the supplement,  Could I  fit a model to this?

#################
# these ordinations contain all samples and are rarefied to 4200 seqs per sample
FS2.mds <- metaMDS(FS2.bray2, k = 2,trymax = 100, autotransform = FALSE)
FS2.mds$stress

FS2.nmds <-as.data.frame(FS2.mds$points)
FS2.nmds$group <- rownames(FS2.nmds)



FS2.metanmds <- merge(meta, FS2.nmds, by = 'group')
FS2.metanmds$design <- factor(FS2.metanmds$design)
FS2.metanmds$treatment <- factor(FS2.metanmds$treatment)
FS2.metanmds$tissue <- factor(FS2.metanmds$tissue)


FS2.metanmds$treatmentXday <- paste(FS2.metanmds$treatment, FS2.metanmds$day, sep = ' day ')
FS2.metanmds$treatmentXday <- factor(FS2.metanmds$treatmentXday)

######
FS2.nmds$group == FS2.metanmds$group
FS2.metanmds$group <- as.character(FS2.metanmds$group)
FS2.metanmds <- FS2.metanmds[match(FS2.nmds$group,FS2.metanmds$group),] # ok I think this works as long as $group isnt a factor...
FS2.nmds$group == FS2.metanmds$group
############
# this is where the elipse weirdness happens, sometimes the .mds and .metanmds have different ordering....


ord <- ordiellipse(FS2.mds,FS2.metanmds$design, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean <- aggregate(FS2.metanmds[,10:11], list(group=FS2.metanmds$design), mean)



df_ell <- data.frame()
for (d in levels(FS2.metanmds$design)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds[FS2.metanmds$design == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}


############ this one is for all timepoints and tissues ###########
############## I think I should do this except for day 21 only..... ##################


df_ell$treatment <- NA
df_ell$treatment[grep('RPS', df_ell$group)] <- 'RPS'
df_ell$treatment[grep('control', df_ell$group)] <- 'control'
df_ell$tissue <- gsub('FS2_(.*)_[0-9]+_.*', '\\1', df_ell$group)
df_ell$day <-gsub('FS2_(.*)_([0-9]+)_.*', '\\2', df_ell$group)


colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')
FS2.metanmds <- merge(FS2.metanmds, NMDS.mean , by='design')

df_ell$treatmentXday <- paste(df_ell$treatment, df_ell$day, sep = ' ')
FS2.metanmds$day <- factor(FS2.metanmds$day)

p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day %in% c(0,21)), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.7, alpha = .75) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, color=treatment), size = .3) + 
  geom_path(data = subset(df_ell, day %in% c(0, 21) & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment), size=1) +
  ggtitle('Feces community similarity', subtitle = 'Day 0 vs Day 21\nPERMANOVA: control: F=10.8, p=0.0001; RPS: F=11.5, p=0.0001') + 
  scale_color_brewer(palette = 'Dark2') + labs(caption = 'Ordination stress = 0.165')

p

# adonis stats info for D0 vs D21 for each group.  Do I need to lump?

# FS2_feces_0_control vs FS2_feces_21_control 10.7868015 0.30141843 0.0001000     0.0435
# FS2_feces_21_control vs FS2_feces_0_RPS 11.5118295 0.31529040 0.0001000     0.0435

##########  ##########

p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 12), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==12  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 12') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p


p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 15), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==15  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 15') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p


p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 19), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==19  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 19') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p

p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 21), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==21  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 21') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p


### dont know if I need these daily ordinations above...  #####

groupxday <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxday <- filter(groupxday, day %in% c(0:22) & tissue == 'feces')


FS2.metanmds %>% filter(tissue == 'feces' & day %in% c(0:21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_path(aes(group=treatment, color=treatment), size = 2, alpha=.7) +
  geom_point(aes(color=treatment), size=7) +
  geom_text(data=groupxday, aes(x=groupX, y=groupY, label=day))+
  scale_color_brewer(palette = 'Dark2')  + ggtitle('Fecal community structure over time', subtitle = 'Treatment group centroids') + ylab('NMDS2') + xlab('NMDS1')

# the above one is nice for showing progression over time maybe
##### #####

groupxtiss <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxtiss <- filter(groupxtiss, day == 21)

FS2.metanmds %>% filter(day %in% c(21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_point(aes(color=treatment), size=8) +
  geom_path(data = subset(df_ell, day == 21), aes(x=NMDS1, y=NMDS2, group=group, color=treatment), size=1) +
  geom_text(data=groupxtiss, aes(x=groupX, y=groupY, label=tissue))+
  scale_color_brewer(palette = 'Dark2') + labs(x='NMDS 1', y= 'NMDS 2', caption = 'Ordination stress: 0.16') +
  ggtitle('Bray-curtis based NMDS ordination of tissues at day 21') 



# I really like this one.  I think I will keep it.
# do I need stats here?

###################################### ###########################################

#PERMANOVA with Adonis

x <- FS2.4200
factors <- meta$design


PW.Adonis <- pairwise.adonis(x,factors,sim.method="bray",p.adjust.m = "bonferroni")

# write.table(PW.Adonis,"Adonis-Results.csv",sep=",")


fecVSfec <- grep('.*_feces_.* vs .*_feces_.*' ,PW.Adonis$pairs)
colVScol <- grep('.*_colon_.* vs .*_colon_.*' ,PW.Adonis$pairs)
cecVScec <- grep('.*_cecum_.* vs .*_cecum_.*' ,PW.Adonis$pairs)
rVSr <- grep('.*_RNA_.* vs .*_RNA_.*' ,PW.Adonis$pairs)
ilVSil <- grep('.*_ileum_.* vs .*_ileum_.*' ,PW.Adonis$pairs)

good <- PW.Adonis[c(fecVSfec, colVScol, cecVScec, rVSr, ilVSil),]

# FS2_feces_0_control vs FS2_feces_21_control 10.7868015 0.30141843 0.0001000     0.0435
# FS2_feces_21_control vs FS2_feces_0_RPS 11.5118295 0.31529040 0.0001000     0.0435


fecsametime <- grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", good$pairs)
ilsametime <- grep("FS2_ileum_([0-9]+)_control vs FS2_ileum_\\1_RPS", good$pairs)
colsametime <- grep("FS2_colon_([0-9]+)_control vs FS2_colon_\\1_RPS", good$pairs)
cecsametime <- grep("FS2_cecum_([0-9]+)_control vs FS2_cecum_\\1_RPS", good$pairs)

good2 <- good[c(fecsametime, ilsametime, colsametime, cecsametime),]
#write.table(good2,"Adonis-Results.csv",sep=",") 

pwadon <- read.table("Adonis-Results.csv",sep=",")
just_poop <- pwadon[grep('.*_feces_.* vs .*_feces_.*', pwadon$pairs),]
poop_time <- just_poop[grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", just_poop$pairs),]
poop_time$day <- c(0,12,15,19,21,28,35,42)
poop_time$p.value <- round(poop_time$p.value, 4)

filter(poop_time, day <22)%>%
  ggplot(aes(x=day, y=F.Model)) +
  geom_line()+geom_point(size=2)+ geom_vline(xintercept = 12, color='purple') + geom_text(aes(y = F.Model + .2, label = paste('p =', p.value))) +
  ggtitle('Dissimilarity of fecal microbiota over time', subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint, how different are the two diets at each timepoint? ') + 
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%')

# this is nice...maybe add dispersion info on this too?

###  relative to D0 for each treatment  ###

reltoD0 <- grep('FS2_feces_0_control vs FS2_feces_[0-9]+_control', just_poop$pairs)
reltoD0 <- c(reltoD0, grep('FS2_feces_0_RPS vs FS2_feces_[0-9]+_RPS', just_poop$pairs))

R2D0 <- just_poop[reltoD0,]
R2D0 <- add_row(.data = R2D0, pairs = 'FS2_feces_0_control vs FS2_feces_0_control', F.Model = 0, p.value = 0, R2 = 0, p.adjusted = 0)
R2D0 <- add_row(.data = R2D0, pairs = 'FS2_feces_0_RPS vs FS2_feces_0_RPS', F.Model = 0, p.value = 0, R2 = 0, p.adjusted = 0)
R2D0$treatment <- ifelse(grepl('control', R2D0$pairs), yes = 'control', no = 'RPS')
R2D0$day <- gsub('.* vs FS2_feces_([0-9]+)_.*', '\\1', R2D0$pairs)
R2D0$day <- as.numeric(R2D0$day)


filter(R2D0, day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=F.Model, color = treatment)) + geom_point(size=2.5) + geom_line() + geom_vline(xintercept = 12)+
  ggtitle('Dissimilarity of fecal microbiota compared to Day 0 over time',
          subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint,\nhow different are the two diets from their original community structures at each timepoint? ') + 
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%') + scale_color_brewer(palette = 'Dark2')

###### Dont think I need this one.... ##########


########################################################################################################
########################################################################################################
########################################################################################################


####### phylum level grouping #######
setwd('~/FS2/16s/V4/')
meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)


otu <- import_mothur(mothur_shared_file = 'V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'V4.final.taxonomy')


phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]



FS2 <- phyloseq(otu, taxo)                                    # generates the initial phyloseq object for FS1
FS2 <- merge_phyloseq(FS2, phy_meta)  
colnames(tax_table(FS2)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.


FS2 <- subset_samples(FS2, experiment =='FS2')
FS2.phylum <- tax_glom(FS2, 'Phylum')
#FS2.phylum@tax_table


FS2.phylum.OTU <- as.data.frame(t(FS2.phylum@otu_table))
rowSums(FS2.phylum.OTU)
FS2.phylum.relabund <- FS2.phylum.OTU/rowSums(FS2.phylum.OTU)
rowSums(FS2.phylum.relabund)
meta <- meta[meta$group %in% rownames(FS2.phylum.relabund),]
#boxplot(FS2.phylum.relabund$Otu00001~meta$design)
FS2.phylum@tax_table@.Data[,2]
colnames(FS2.phylum.relabund) <- FS2.phylum@tax_table@.Data[,2]


FS2.phylum.relabund$group <- rownames(FS2.phylum.relabund)

FS2.phylum.relabund <- merge(FS2.phylum.relabund, meta, by = 'group')

colnames(FS2.phylum.relabund)
FS2.phylum.relabund$day
FS2.phylum.relabund <- melt(FS2.phylum.relabund, id.vars = c(1, 22:29) )

FS2.phylum.relabund$pig_num <- factor(FS2.phylum.relabund$pig_num)
FS2.phylum.relabund$day <- factor(FS2.phylum.relabund$day)

goodphy <- levels(FS2.phylum.relabund$variable)[c(1:8,13)]

#FS2.phylum


filter(FS2.phylum.relabund, tissue == 'feces' & day %in% c(0,21) & variable %in% goodphy) %>%
  ggplot(aes(x=day, y=value, fill=treatment))+ geom_boxplot(aes(group=design)) +
  scale_fill_brewer(palette = 'Dark2') + expand_limits(y=0) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  ggtitle('Weaning effects on major phyla in fecal microbiota') + 
  ylab('proportion of total community')

# boom!  use it next to ordination showing D0 to D21 shift

D21phy <- levels(FS2.phylum.relabund$variable)[c(1,2,3,4,5,7,8,10,12)]

filter(FS2.phylum.relabund, tissue == 'feces' & day == 21 & variable %in% D21phy) %>% 
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Feces') + 
  ylab('proportion of total community')


filter(FS2.phylum.relabund, tissue == 'cecum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal mucosa') +
  ylab('proportion of total community')

filter(FS2.phylum.relabund, tissue == 'colon' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Colonic mucosa') +
  ylab('proportion of total community')


filter(FS2.phylum.relabund, tissue == 'ileum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Ileal mucosa') + 
  ylab('proportion of total community')

filter(FS2.phylum.relabund, tissue == 'cec_cont_RNA' & day == 21 & variable %in% D21phy) %>% 
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal contents (RNA)') + 
  ylab('proportion of total community')

# These are good, but maybe its too much... can I concentrate this at all?



############### BUTNET ###############

ptm <- proc.time()




setwd('~/FS2/correlation/but_corr/')
system('ls')
library(reshape2)
library(ggplot2)
library(geomnet)
library(tidyr)
library(dplyr)
library(ccrepe)
library(vegan)
#library(igraph)
library(Hmisc)
############# ##############


rcorr_to_ggnet <- function(rcorr.list, pcut=0.05, spearcut=0.6){
  pval <- as.data.frame(rcorr.list$P)
  sim <- as.data.frame(rcorr.list$r)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  
  return(sigcor)
  
  
}

ccrepe_to_ggnet <- function(ccrepe.obj, pcut=0.01, spearcut=0.6){
  pval <- as.data.frame(ccrepe.obj$p.values)
  sim <- as.data.frame(ccrepe.obj$sim.score)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spearman')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spearman > spearcut,]
  
  sigcor2 <- sigcor[,c(2,1,3,4)]
  colnames(sigcor2) <- c('from', 'to', 'pval', 'spearman')
  sigglies <- rbind(sigcor, sigcor2)
  return(sigglies)
  
  
}

ccrepe_to_ggnet2 <- function(ccrepe.obj, pcut=0.05, spearcut=0.6){
  pval <- as.data.frame(ccrepe.obj$p.values)
  sim <- as.data.frame(ccrepe.obj$sim.score)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to', factorsAsStrings=TRUE)
  sim <- melt(sim, id.vars = 'to', factorsAsStrings=TRUE)
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  sigcor$from <- as.character(sigcor$from)
  sigcor$to <- as.character(sigcor$to)
  return(sigcor)
  
}


# I could make a generic function that accepted a matrix of pvalues and a matrix of correlation coeficcients instead of an
# object/list from rcorr or ccrepe...

extract_mothur_tax <- function(filename){
  
  tax <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  tna <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
  tax[tna] <- NA
  rowcount <- 1
  
  for (i in strsplit(tax$Taxonomy, split = ';')){
    levelcount <- 1
    for (x in i){
      tax[rowcount,3+levelcount] <- x
      levelcount <- levelcount + 1
      print(x)
    }
    rowcount <- rowcount+1
  }
  return(tax)
}

otu_tax_labels <- function(mothur_tax){
  tax <- mothur_tax
  swip <- tax$OTU
  names(swip) <- paste(tax$Genus, tax$OTU, sep = ':')
  names(swip)[grepl("NA", names(swip))] <- paste(tax$Family[grepl("NA", names(swip))], tax$OTU[grepl("NA", names(swip))], sep = ':')
  swap <- names(swip)
  names(swap) <- swip
  return(swap)
}





prune_graph <- function(fortified.edgelist, node.dataframe, min.vert = 1){
  #hopefully this will take a fortified edgelist correlation graph from geomnet, convert it to an igraph object, 
  #make sure it's undirected? is this necessary?
  #decompose the graph, throwing away subnetworks that have less than min.vert vertexes
  #recompose the graph and convert it back to the fortified edgelist format
  library(igraph)
  filtered <- graph_from_data_frame(fortified.edgelist, directed = FALSE)
  filtered <- as.undirected(filtered, mode= "collapse")
  filtered <- delete.vertices(filtered, 'NA')
  filtered <- decompose(filtered, mode = "weak", max.comps = NA, min.vertices = min.vert)
  temp <- make_empty_graph(n=0, directed=FALSE)
  
  for (x in 1:length(filtered)){
    temp <- temp %du% filtered[x]
  }
  
  filtered <- as_long_data_frame(temp)
  rm(temp)
  
  filtered$from <- filtered[,length(colnames(filtered))-1]
  filtered$to <- filtered[,length(colnames(filtered))]
  
  filtered <- filtered[,1:4]
  filtered <- fortify(as.edgedf(filtered), node.dataframe)
  return(filtered)
}




gather_nodes <- function(x, typ=NA){
  x <- as.data.frame(colnames(x))
  colnames(x)[1] <- 'node'
  x$type <- typ
  return(x)
  
}


##### enriched features ####

# setwd("target_dir/")

file_list <- list.files('~/FS2/16s/V4/OTUs', full.names = TRUE)
gsub('.*/(.*)_OTUs.txt', '\\1', file_list)

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
    dataset$file <- gsub('.*/(.*)_OTUs.txt', '\\1', file)
    dataset$node <- rownames(dataset)
    
    
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    temp_dataset$file <- gsub('.*/(.*)_OTUs.txt', '\\1', file)
    temp_dataset$node <- rownames(temp_dataset)
    
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}


########## but taxonomy ############



D21 <- read.table('~/FS2/16s/V4/OTUs/D21_OTUs.txt', header = TRUE, stringsAsFactors = FALSE)


xx <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
xx <- unlist(strsplit(xx, split = ' '))



blast <- read.table('../../but_amp/from_sci/butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
shared2 <- read.table('../../but_amp/from_sci/but2.shared', header = TRUE, stringsAsFactors = FALSE)

colnames(blast) <- xx

blast$tit1 <- gsub('.* \\[(.*)\\]', '\\1', blast$salltitles)
#blast$tit1 <- paste(blast$tit1, blast$pident, sep = ' ')
blast$tit1

blast$otu <- colnames(shared2)[-c(1,2,3)]
blast$gen <- gsub('(\\[?[A-Za-z]+\\]?).*', '\\1', blast$tit1)
blast$lab <- paste(blast$gen, blast$otu, sep = ' ')


#hist(rowSums(shared2[,-c(1,2,3)]), breaks = 100)

#sort(rowSums(shared2[,-c(1,2,3)]))
#rowSums(shared2[,-c(1,2,3)]) > 2000

#shared2$Group == meta$group


#### this chunk creates a named vector that you can use to swap out the generic otuxxx labels for the blast results

butswip <- blast$otu
names(butswip) <- blast$lab
butswap <- names(butswip)
names(butswap) <- butswip

butswap[blast$otu]
#

################# Read in data ####################

tax <- extract_mothur_tax('~/FS2/16s/V4/V4.final.taxonomy')
tax <- tax[,-c(2,3)]
swap <- otu_tax_labels(tax)
# colnames(cecbact) <- swap[colnames(cecbact)]

# 16S data #

meta <- read.table('~/FS2/correlation/but_corr/16S_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE) 

otu <- read.table('~/FS2/correlation/but_corr/16S_shared_forcorr.txt', header = TRUE) %>% filter(group %in% meta$group)



rownames(otu) <- meta$sample
# glycoside hydrolase data #
glyc <- read.table('../glycoside_hydrolase_matrix.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
glyc$sample <- paste(glyc$tissue, glyc$day, glyc$pig, sep = '_')

# vfa data #



##########################  VFA ANALYSIS  ####################

vfa.cec <- read.table('~/FS2/correlation/but_corr/cec_redo_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfa.fec <- read.table('~/FS2/correlation/but_corr/fecalplasma_R.txt', header = TRUE, stringsAsFactors = FALSE) #%>% 
#filter(tissue == 'feces')

vfas <- rbind(vfa.cec, vfa.fec)
vfas[,2:16] <- vfas[,2:16]*3
vfas$total <- rowSums(vfas[,2:16])
vfas$BCFA <- vfas$isobutyrate + vfas$isovalerate
colnames(vfas)[8] <- 'lactate'


vfas$sample <- paste(vfas$tissue, 21, vfas$pig_num, sep = '_')
vfas$treatment <- ifelse(vfas$pig_num %in% c(67,68,69,70,71,72,73,81,82,83,84,85,86,87), 'control', 'RPS')
vfas$design <- paste(vfas$tissue, vfas$treatment, sep = '_')
colnames(vfas)
#vfas.melt$variable %in% c('heptanoate', 'lactate1', )
vfas.melt <- melt(data = vfas, id.vars = c(1,17, 20:22)) 


vfas.melt %>% filter(tissue == 'feces') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
vfas.melt %>% filter(tissue == 'portal') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
vfas.melt %>% filter(tissue == 'cecum') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')

vfas.melt %>% filter(variable %in% c('propionate', 'butyrate', 'lactate', 'total') & tissue != 'portal') %>%
  ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  geom_boxplot() + facet_wrap(~variable, scales = 'free') + scale_fill_brewer(palette = 'Dark2') + ylab('Concentration (mM)') + ggtitle('VFA concentrations at day 21')
# but data #

but <- read.table('~/FS2/correlation/but_corr/but_shared_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
butmeta <- read.table('~/FS2/correlation/but_corr/but_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
but$group == butmeta$group
butmeta$sample <- paste(butmeta$tissue, butmeta$day, butmeta$pig_num, sep = '_')
rownames(but) <- butmeta$sample
#rownames(otu)

########## getting everything in order ############


meta <- filter(meta, sample %in% butmeta$sample)
meta <- meta[match(butmeta$sample, meta$sample),]

otu <- otu[otu$group %in% meta$group,]
otu <- otu[match(meta$group, otu$group),]
otu <- na.omit(otu)

otu <- otu[,which(names(otu) != 'group')]
but <- but[,which(names(but) != 'group')]




sort(rowSums(otu))
otu <- otu[rowSums(otu)>1000,]
meta <- meta[which(meta$sample %in% rownames(otu)),]
butmeta <- butmeta[which(butmeta$sample %in% meta$sample),]
butmeta$sample == meta$sample

but <- but[which(rownames(but) %in% butmeta$sample),]

rownames(but) == butmeta$sample

otu <- otu/rowSums(otu)
but <- but/rowSums(but)

colSums(otu)
colSums(but)
# CHANGE DIVIDED BY TO MATCH NUMBER OF SAMPLES


otu <- otu[,((colSums(otu)/186)*100) >0.01] # removes 16S otus with less than 0.01% average abundance
but <- but[,((colSums(but)/186)*100) >0.01] # removes but otus with less than 0.01% average abundance

colnames(otu) <- swap[colnames(otu)]
colnames(but) <- butswap[colnames(but)]


######## glyc wrangle  ########
glyc$sample <- gsub('cec_cont_(.*)', 'cec_cont_RNA_\\1', glyc$sample)

rownames(glyc) <- glyc$sample
#glyc.feces <- filter(glyc, tissue =='feces')
#rownames(glyc) <- glyc$sample
glyc <- glyc[,c(2:8)]
glyc <- glyc/rowSums(glyc)  # COMMENT THIS OUT IF YOU DON'T WANT RELATIVE ABUND FOR GLYC


rownames(glyc) == rownames(but)

rownames(otu) %in% rownames(glyc)
otu <- otu[which(rownames(otu) %in% rownames(glyc)),]
but <- but[which(rownames(but) %in% rownames(glyc)),]
glyc <- glyc[which(rownames(glyc) %in% rownames(but)),]
glyc <- glyc[match(rownames(otu),rownames(glyc)),]


rownames(but) == rownames(otu)
rownames(glyc) == rownames(otu)

####### vfa wrangle ########

vfa.cec$sample <- paste('cec_cont_RNA', 21, vfa.cec$pig_num, sep = '_')
rownames(vfa.cec) <- vfa.cec$sample

vfa.cec <- vfa.cec[,c(2:11,13,14,15)]
vfa.cec <- vfa.cec/rowSums(vfa.cec)  # COMMENT OUT IF DONT WANT RELABUND

# ############################# D21 feces ########################
# 
# meta.feces <- filter(meta, tissue =='feces' & day == 21)        # subset metadata to just d21 feces
# butmeta.feces <- filter(butmeta, tissue =='feces' & day == 21)
# vfa.feces <- filter(vfas, tissue == 'feces') %>% select(ends_with('e'))
# rownames(vfa.feces) <- vfa.feces$sample
# vfas.feces <- vfa.feces[,c(1:8,10:13)]
# 
# meta.feces$sample == butmeta.feces$sample                       # make sure all samples are in the same order
# 
# otu.feces <- otu[rownames(otu) %in% meta.feces$sample,]             # subset data to just d21 feces
# but.feces <- but[rownames(but) %in% butmeta.feces$sample,]
# 
# otu.feces <- otu.feces[rownames(otu.feces) %in% rownames(vfas.feces),]
# but.feces <- but.feces[rownames(but.feces) %in% rownames(vfas.feces),]
# 
# vfas.feces <- vfas.feces[match(rownames(vfas.feces), rownames(otu.feces)),]
# 
# rownames(otu.feces) == rownames(but.feces)
# rownames(vfas.feces) == rownames(but.feces)
# 
# 
# 
# 
# #otu.feces <- otu.feces[,-which(names(otu.feces) == 'group')]   # remove the group column from the otu tables
# #but.feces <- but.feces[,-which(names(but.feces) == 'group')]
# 
# #colnames(otu.feces) <- swap[colnames(otu.feces)]
# #colnames(but.feces) <- butswap[colnames(but.feces)]
# 
# 
# but_16S_D21fec <- ccrepe(otu.feces, but.feces, verbose = TRUE, min.subj = 7)
# 
# #####
# 
# 
# #rcorr(otu.feces, vfas.feces)
# #vfa_16S_D21fec <- ccrepe(otu.feces, vfas.feces, verbose = TRUE, min.subj = 7)
# #vfa_16S_D21fec <- ccrepe(otu.feces, vfas.feces, verbose = TRUE, min.subj = 7)
# 
# 
# 
# 
# D21fec.sigs <- ccrepe_to_ggnet2(ccrepe.obj = but_16S_D21fec, pcut = 0.05, spearcut = 0.7)
# 
# #D21vfa_otu.sigs <- ccrepe_to_ggnet2(ccrepe.obj = vfa_16S_D21fec, pcut = 0.05, spearcut = 0.6)
# 
# # grep('Otu00087', D21fec.sigs$to)
# # D21fec.sigs[153,]
# 
# nodes <- rbind(gather_nodes(otu.feces, '16S'),
#                gather_nodes(but.feces, 'but'))
# all
# 
# all <- rbind(D21fec.sigs)
# D21 <- fortify(as.edgedf(all), nodes)
# 
# 
# set.seed(12812)
# 
# p <- ggplot(D21, aes(from_id = from_id, to_id = to_id, label=from, color=type)) +
#   geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)),
#            aes(color = type, label = from_id),
#            linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
#            repel = TRUE, fontsize=3.5, singletons = FALSE,labelcolour="black",
#            labelgeom = 'text') +
#   theme_net()
# p
# 
# 
# 

### all timepoints only but and 16s ###

#otu <- otu[,-which(names(otu) == 'group')]
#but <- but[,-which(names(but) == 'group')]




###############################

rownames(but) == rownames(otu) # ok good

rownames(but)

but.cec <- but[grep('cec', rownames(but)),]
otu.cec <- otu[grep('cec', rownames(otu)),]
glyc.cec <- glyc[grep('cec', rownames(glyc)),]
vfa.cec <- vfa.cec[rownames(vfa.cec) %in% rownames(otu.cec),]

rownames(but.cec) == rownames(otu.cec) # ok good
rownames(glyc.cec) == rownames(otu.cec) # ok good
rownames(vfa.cec) == rownames(otu.cec) # ok good


but_16S_cec <- ccrepe(but.cec, otu.cec, verbose = TRUE, min.subj = 3)
cec_16S <- ccrepe(otu.cec, verbose = TRUE, min.subj = 3)
cec_but <- ccrepe(but.cec, verbose = TRUE, min.subj = 3)
vfa_glyc <- ccrepe(vfa.cec, glyc.cec, verbose = TRUE, min.subj = 3)
vfa_16S <- ccrepe(vfa.cec, otu.cec, verbose = TRUE, min.subj = 3)
vfa_but <- ccrepe(vfa.cec, but.cec, verbose = TRUE, min.subj = 3)
glyc_16S <- ccrepe(glyc.cec, otu.cec, verbose = TRUE, min.subj = 3)
glyc_but <- ccrepe(glyc.cec, but.cec, verbose = TRUE, min.subj = 3)


but_16S_cec$sim.score

### rcorr() stuff ###
#glyc_vfa_cec <- rcorr(as.matrix(vfa.cec), as.matrix(glyc.cec))
#glyc_16S_cec <- rcorr(as.matrix(glyc.cec), as.matrix(otu.cec))
#glyc_but_cec <- rcorr(as.matrix(glyc.cec), as.matrix(but.cec))

#vfa_but_cec <- rcorr(as.matrix(vfa.cec), as.matrix(but.cec))
#vfa_otu_cec <- rcorr(as.matrix(vfa.cec), as.matrix(otu.cec))

#glyc_vfa_cec.sigs <- rcorr_to_ggnet(glyc_vfa_cec, pcut = 0.05, spearcut = 0.60)
#glyc_16S_cec.sigs <- rcorr_to_ggnet(glyc_16S_cec, pcut = 0.05, spearcut = 0.60)
#glyc_but_cec.sigs <- rcorr_to_ggnet(glyc_but_cec, pcut = 0.05, spearcut = 0.60)
#vfa_but_cec.sigs <- rcorr_to_ggnet(vfa_but_cec, pcut = 0.05, spearcut = 0.60)
#vfa_otu_cec.sigs <- rcorr_to_ggnet(vfa_otu_cec, pcut = 0.05, spearcut = 0.60)


#glyc_16S_cec.sigs <- glyc_16S_cec.sigs[!(grepl('Otu', glyc_16S_cec.sigs$from) & grepl('Otu', glyc_16S_cec.sigs$to)),]
#glyc_but_cec.sigs <- glyc_but_cec.sigs[!(grepl('Otu', glyc_but_cec.sigs$from) & grepl('Otu', glyc_but_cec.sigs$to)),]
#vfa_but_cec.sigs <- vfa_but_cec.sigs[!(grepl('Otu', vfa_but_cec.sigs$from) & grepl('Otu', vfa_but_cec.sigs$to)),]
#vfa_otu_cec.sigs <- vfa_otu_cec.sigs[!(grepl('Otu', vfa_otu_cec.sigs$from) & grepl('Otu', vfa_otu_cec.sigs$to)),]
###

but_16S_cec.sigs <- ccrepe_to_ggnet2(ccrepe.obj = but_16S_cec, pcut = 0.05, spearcut = 0.6)
cec_16S.sigs <- ccrepe_to_ggnet2(ccrepe.obj = cec_16S, pcut = 0.05, spearcut = 0.6)
cec_but.sigs <- ccrepe_to_ggnet2(ccrepe.obj = cec_but, pcut = 0.05, spearcut = 0.6)
vfa_glyc.sigs <- ccrepe_to_ggnet2(ccrepe.obj = vfa_glyc, pcut = 0.05, spearcut = 0.6)
vfa_16S.sigs <- ccrepe_to_ggnet2(ccrepe.obj = vfa_16S, pcut = 0.05, spearcut = 0.6)
vfa_but.sigs <- ccrepe_to_ggnet2(ccrepe.obj = vfa_but, pcut = 0.05, spearcut = 0.6)
vfa_16S.sigs <- ccrepe_to_ggnet2(ccrepe.obj = vfa_16S, pcut = 0.05, spearcut = 0.6)
glyc_16S.sigs <- ccrepe_to_ggnet2(ccrepe.obj = glyc_16S, pcut = 0.05, spearcut = 0.6)
glyc_but.sigs <- ccrepe_to_ggnet2(ccrepe.obj = glyc_but, pcut = 0.05, spearcut = 0.6)
#rownames(as.matrix(otu)) == rownames(as.matrix(but))

#good_nodes <- c(but_16S_cec.sigs$from, but_16S_cec.sigs$to)

#cec_16S.sigs <- cec_16S.sigs[cec_16S.sigs$from %in% good_nodes & cec_16S.sigs$to %in% good_nodes,]
#cec_but.sigs <- cec_but.sigs[cec_but.sigs$from %in% good_nodes & cec_but.sigs$to %in% good_nodes,]

#glyc_16S_cec.sigs <- glyc_16S_cec.sigs[!(grepl('Otu', glyc_16S_cec.sigs$from) & grepl('Otu', glyc_16S_cec.sigs$to)),]


nodes <- rbind(gather_nodes(otu.cec, '16S'), 
               gather_nodes(but.cec, 'but'),
               gather_nodes(vfa.cec, 'vfa'),
               gather_nodes(glyc.cec, 'hydrolase'))


all <- rbind(cec_but.sigs, but_16S_cec.sigs)
all <- fortify(as.edgedf(all), nodes)



filtered3 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 3)
filtered4 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 4)
filtered5 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 5)
filtered6 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 6)
filtered7 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 7)


p <- ggplot(all, aes(from_id = from_id, to_id = to_id, label=from, color=type)) + 
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
           aes(color = type, label = from_id),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
           labelgeom = 'text') +
  theme_net()

p








colnames(otu)[84]

grep('Otu00087', colnames(otu))
otu[,84]
meta$treatment

boxplot(otu[,84]~meta$treatment)



net <- ggplot_build(p)
net <- as.data.frame(ggplot_build(p)$data)


enrich <- dataset %>% select(newp, Treatment, file, node)


uenrich <- unique(enrich)
colnames(uenrich)[4] <- 'from'

uenrich$from <- swap[uenrich$from]
net2 <- merge(uenrich, net, by = 'from', all = TRUE)

net2 <- net2[!is.na(net2$x),]
net2$type <- NA
net2$type[grep('Otu[0-9][0-9][0-9]', net2$from)] <- 'but'
net2$type[grep('Otu[0-9][0-9][0-9][0-9][0-9]', net2$from)] <- '16S'


ggplot(net2, aes(x=x, y=y, xend=xend, yend=yend)) +  geom_segment() +
  geom_point(aes(color = Treatment), size=12, alpha=.3) +
  geom_point(aes(color=type), size=3) + geom_text(aes(label=label), size=3)

net2$label <- butswap[net2$label]
net2$label[is.na(net2$label)] <- swap[net2$from[is.na(net2$label)] ]


net2[157,]




proc.time() - ptm


############################# IMCOR2 ####################


ptm <- proc.time()

setwd('~/FS2/correlation/')

library(ggplot2)
library(vegan)
library(ccrepe)
library(reshape2)
library(geomnet)
library(Hmisc)
#library(tidyverse)
library(ggrepel)
library(dplyr)




# functions ---------------------------------------------------------------


prune_graph <- function(fortified.edgelist, node.dataframe, min.vert = 1){
  #hopefully this will take a fortified edgelist correlation graph from geomnet, convert it to an igraph object, 
  #make sure it's undirected? is this necessary?
  #decompose the graph, throwing away subnetworks that have less than min.vert vertexes
  #recompose the graph and convert it back to the fortified edgelist format
  library(igraph)
  filtered <- graph_from_data_frame(fortified.edgelist, directed = FALSE)
  filtered <- as.undirected(filtered, mode= "collapse")
  filtered <- delete.vertices(filtered, 'NA')
  filtered <- decompose(filtered, mode = "weak", max.comps = NA, min.vertices = min.vert)
  temp <- make_empty_graph(n=0, directed=FALSE)
  
  for (x in 1:length(filtered)){
    temp <- temp %du% filtered[x]
  }
  
  filtered <- as_long_data_frame(temp)
  rm(temp)
  
  filtered$from <- filtered[,length(colnames(filtered))-1]
  filtered$to <- filtered[,length(colnames(filtered))]
  
  filtered <- filtered[,1:4]
  filtered <- fortify(as.edgedf(filtered), node.dataframe)
  return(filtered)
}



gather_nodes <- function(x, typ=NA){
  x <- as.data.frame(colnames(x))
  colnames(x)[1] <- 'node'
  x$type <- typ
  return(x)
  
}

rcorr_to_ggnet <- function(rcorr.list, pcut=0.05, spearcut=0.6){
  pval <- as.data.frame(rcorr.list$P)
  sim <- as.data.frame(rcorr.list$r)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  
  return(sigcor)
  
  
}

ccrepe_to_ggnet <- function(ccrepe.obj, pcut=0.01, spearcut=0.6){
  pval <- as.data.frame(ccrepe.obj$p.values)
  sim <- as.data.frame(ccrepe.obj$sim.score)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spearman')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spearman > spearcut,]
  
  sigcor2 <- sigcor[,c(2,1,3,4)]
  colnames(sigcor2) <- c('from', 'to', 'pval', 'spearman')
  sigglies <- rbind(sigcor, sigcor2)
  return(sigglies)
  
  
}

ccrepe_to_ggnet2 <- function(ccrepe.obj, pcut=0.01, spearcut=0.6){
  pval <- as.data.frame(ccrepe.obj$p.values)
  sim <- as.data.frame(ccrepe.obj$sim.score)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  
  return(sigcor)
  
}


# I could make a generic function that accepted a matrix of pvalues and a matrix of correlation coeficcients instead of an
# object/list from rcorr or ccrepe...

extract_mothur_tax <- function(filename){
  
  tax <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  tna <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
  tax[tna] <- NA
  rowcount <- 1
  
  for (i in strsplit(tax$Taxonomy, split = ';')){
    levelcount <- 1
    for (x in i){
      tax[rowcount,3+levelcount] <- x
      levelcount <- levelcount + 1
      print(x)
    }
    rowcount <- rowcount+1
  }
  return(tax)
}

otu_tax_labels <- function(mothur_tax){
  tax <- mothur_tax
  swip <- tax$OTU
  names(swip) <- paste(tax$Genus, tax$OTU, sep = ':')
  names(swip)[grepl("NA", names(swip))] <- paste(tax$Family[grepl("NA", names(swip))], tax$OTU[grepl("NA", names(swip))], sep = ':')
  swap <- names(swip)
  names(swap) <- swip
  return(swap)
}



# read in data ------------------------------------------------------------


# this has cd3- population in it as well #
allflow <- read.table('FS2_all_flow.txt', header = TRUE, as.is = TRUE, check.names = FALSE)
allflow <- filter(allflow, tissue =='cec')
rownames(allflow) <- allflow$pig_num
allflow <- allflow[,-c(1,34:36)]
allflow <- allflow/rowSums(allflow)
sum((colSums(allflow)/14)*100 >.1)


(colSums(allflow[order(colSums(allflow))])/14)*100
allflow <- allflow[,(colSums(allflow)/14)*100 >.1] # removes cell types with less than 0.05% abundance (% total live single cells)

flowr <- allflow
#

# uncomment for only cd3+ cells
# flowc <- read.table("Flow_abs_counts.txt", header = TRUE, comment.char = '', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
# rownames(flowc) <- flowc$Pig
# flowc <- flowc[-15,-c(1,2,3)]
# flowr <- flowc/rowSums(flowc)
# rowSums(flowr)  #checking rows sum to 1
# 
# colnames(flowr)

#

meta <- read.table('~/FS2/16s/V4/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)
meta <- meta[meta$day == 21 & meta$tissue == 'cecum',]

#

# misc #
##################### NEED TO CHANGE THIS SECTION HERE #################
misc <- read.table('miscforcorr.txt', header = TRUE, stringsAsFactors = FALSE, sep = ',', as.is = TRUE, check.names = FALSE)

rownames(misc) <- misc$pig_num
misc <- misc[,-c(1,3,4)]


misc2 <- misc
misc2$treatment <- c(rep('control', 7), rep('RPS', 7))

boxplot(misc2$`ng IgA/mg dry contents`~misc2$treatment)
wilcox.test(misc2$`ng IgA/mg dry contents`~misc2$treatment)



misc2$totmuc <- misc2$NacBDglu.rate + misc2$Bdgalacto.rate + misc2$aLfuco.rate




boxplot(misc2$totmuc~misc2$treatment)

wilcox.test(misc2$totmuc~misc2$treatment)




#


tax <- extract_mothur_tax('~/FS2/16s/V4/V4.final.taxonomy')
swap <- otu_tax_labels(tax)

vfas <- read.table('vfas.txt', header = TRUE, sep = '\t')
vfas <- filter(vfas, day == 21 & location == 'cecum')
rownames(vfas) <- vfas$number
colnames(vfas)[24] <- 'lactate'


cecbact <- read.table('~/FS2/16s/V4/cecum.shared', header = TRUE)
rownames(cecbact) <- meta$pig_num
cecbact <- cecbact[,-c(1,2,3)]
rowSums(cecbact)  #31441 seqs per sample
cecbact <- cecbact/rowSums(cecbact)  # converts to relative abundance
cecbact <- cecbact[,colSums(cecbact)>.0001] #removes OTUs with less than 0.01% abundance across all cecal tissue samples
colnames(cecbact) <- swap[colnames(cecbact)]

cecbactcont <- read.table('~/FS2/16s/V4/cec_cont_RNA.shared', header = TRUE)
rownames(cecbactcont) <- meta$pig_num
cecbactcont <- cecbactcont[,-c(1,2,3)]
rowSums(cecbactcont)  #21805 seqs per sample
cecbactcont <- cecbactcont/rowSums(cecbactcont)  # converts to relative abundance
cecbactcont <- cecbactcont[,colSums(cecbactcont)>.0001] #removes OTUs with less than 0.01% abundance across all cecal tissue samples
colnames(cecbactcont) <- swap[colnames(cecbactcont)]


qPCR <- read.table('Cec_tissue_forcorr.txt', header = TRUE)
qPCR <- aggregate(x=qPCR, data=qPCR, by=list(qPCR$pig_num), FUN=mean)
rownames(qPCR) <- qPCR$pig_num
qPCR <- qPCR[,-c(1,2)]

qPCR$ACTb/qPCR$DEFb1
deltaCT <- qPCR - qPCR$ACTb


deltaCT.m <- as.matrix(deltaCT)
deltaCT.m <- 1/deltaCT.m[,-1]
vfas.m <- as.matrix(vfas[,-c(1:6,12,13, 17, 22,23,25)])
flowrm <- as.matrix(flowr)
cecbactm <- as.matrix(cecbact)
cecbactcontm <- as.matrix(cecbactcont)

deltaCT.m[,8]

misc2$mucexpr <- deltaCT.m[,8]

misc2$mucdeg <- misc2$totmuc / misc2$mucexpr

boxplot(misc2$mucdeg~misc2$treatment)
wilcox.test(misc2$mucdeg~misc2$treatment)

boxplot(misc2$mucexpr~misc2$treatment)
wilcox.test(misc2$mucexpr~misc2$treatment)

miscm <- as.matrix(misc)
miscm <- miscm[,c(4,5,9,10,11,12,13)]

# Correlation calculations ------------------------------------------------


TxB <- ccrepe(x = cecbact, y = flowr, min.subj = 7, verbose = TRUE)
bac <- ccrepe(x = cecbact, min.subj = 7, verbose = TRUE, compare.within.x = TRUE)
im <- ccrepe(x = flowr, min.subj = 7, verbose = TRUE)


min(na.omit(TxB$q.values))
vfaVSflow <- rcorr(vfas.m, flowrm)
vfaVSbact <- rcorr(vfas.m, cecbactm)
qPCRvsFlow <- rcorr(deltaCT.m, flowrm)
qPCRvsBac <- rcorr(deltaCT.m, cecbactm)
miscVSflow <- rcorr(miscm, flowrm)
miscVSqPCR <- rcorr(miscm, deltaCT.m)
miscVSbac <- rcorr(miscm, cecbactm)
miscVSbaccont <- rcorr(miscm, cecbactcontm)
vfaVSmisc <- rcorr(vfas.m, miscm)

vfaVSqPCR <- rcorr(vfas.m, deltaCT.m)





# Gathering node data -----------------------------------------------------

nodes <- rbind(gather_nodes(flowrm, 'T-cell'), 
               gather_nodes(cecbact, '16s'), 
               gather_nodes(vfas.m, 'VFA'),
               gather_nodes(deltaCT.m, 'mRNA'),
               gather_nodes(misc), gather_nodes(cecbactcont, '16S cont'))

nodes$type[grep('CD3-', nodes$node)] <- 'CD3neg'
nodes$type[grep('CD3\\+', nodes$node)] <- 'CD3pos'
# Convert to geomnet format -----------------------------------------------

TxB.sigs <- ccrepe_to_ggnet2(TxB, spearcut = 0.6)
im.sigs <- ccrepe_to_ggnet2(im, spearcut = 0.6)

vf.sig <- rcorr_to_ggnet(vfaVSflow, pcut = 0.05)
qPCR_flow.sig <- rcorr_to_ggnet(qPCRvsFlow, pcut = 0.05)
misc_qPCR.sig <- rcorr_to_ggnet(miscVSqPCR, pcut = 0.05)
misc_flow.sig <- rcorr_to_ggnet(miscVSflow, pcut = 0.05)
vfa_misc.sig <- rcorr_to_ggnet(vfaVSmisc, pcut = 0.05)
vfa_qPCR.sig <- rcorr_to_ggnet(vfaVSqPCR, pcut = 0.05)

# these below need to be pruned #

vfbac.sig <- rcorr_to_ggnet(vfaVSbact, pcut = 0.05)
misc_baccont.sig <- rcorr_to_ggnet(miscVSbaccont, pcut = 0.05)
qPCR_bact.sig <- rcorr_to_ggnet(qPCRvsBac, pcut = 0.05)
misc_bac.sig <- rcorr_to_ggnet(miscVSbac, pcut = 0.05)
bac.sigs <- ccrepe_to_ggnet2(bac, pcut = 0.05, spearcut = 0.6)

psdsd <- miscVSbac$P




#miscVSbaccont
# plot(misc$log.ng.IgA.prot, flowrm[,5])
# cor(misc$AvDailyGain,misc$log.ng.IgA.prot)

# removing stuff from bac.sigs and vf.sigs --------------------------------

# removes tcell-tcell correlations calculated by rcorr()

vf.sig <- vf.sig[!(grepl('CD25', vf.sig$from) & grepl('CD25', vf.sig$to)),] 
misc_flow.sig <- misc_flow.sig[!(grepl('CD25', misc_flow.sig$from) & grepl('CD25', misc_flow.sig$to)),] 
misc_bac.sig<- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),] 
misc_baccont.sig <- misc_baccont.sig[!(grepl('Otu', misc_baccont.sig$from) & grepl('Otu', misc_baccont.sig$to)),] 

vfbac.sig <- vfbac.sig[!(grepl('Otu', vfbac.sig$from) & grepl('Otu', vfbac.sig$to)),] 
grep('unclassified', vfbac.sig$from)
grep('unclassified', vfbac.sig$to)

misc_baccont.sig <- misc_baccont.sig[!(grepl('Otu', misc_baccont.sig$from) & grepl('Otu', misc_baccont.sig$to)),] 
qPCR_bact.sig <- qPCR_bact.sig[!(grepl('Otu', qPCR_bact.sig$from) & grepl('Otu', qPCR_bact.sig$to)),] 
misc_bac.sig <- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),] 





##### ######
# this makes it so that bac to bac correlations are limited to those OTUs which also correlate with other features




bac.sigs <- rbind(bac.sigs[bac.sigs$from %in% TxB.sigs$to,],
                  bac.sigs[bac.sigs$to %in% TxB.sigs$to,],
                  bac.sigs[bac.sigs$from %in% vfbac.sig$to,],
                  bac.sigs[bac.sigs$to %in% vfbac.sig$from,])

qPCR_flow.sig <- qPCR_flow.sig[!(grepl('CD25', qPCR_flow.sig$from) & grepl('CD25', qPCR_flow.sig$to)),]

qPCR_flow.sig


qPCR_bact.sig <- qPCR_bact.sig[!(grepl('Otu', qPCR_bact.sig$from) & grepl('Otu', qPCR_bact.sig$to)),]


all <- rbind(TxB.sigs,
             vf.sig,
             im.sigs,
             qPCR_flow.sig,
             misc_flow.sig,
             misc_bac.sig,
             vfa_misc.sig, 
             vfa_qPCR.sig)

all <- fortify(as.edgedf(all), nodes)

all$from_id


all$label <- gsub('(.*):Otu[0-9]+', '\\1', all$from_id)

# all$label <- gsub('(CD3./CD4./CD8a./)(FoxP3./CD25.)', '\\1\n\\2', all$label)
# 
# allnet <- network(all)
# get.neighborhood(allnet)


# this removes subnetworks with less than 3 nodes

remove <- c(grep('IL.10', all$from_id),
            grep('TGF.beta', all$from_id),
            grep('IL1b', all$from_id),
            grep('oxalate', all$from_id),
            grep('formate', all$from_id),
            grep('succinate', all$from_id),
            grep('MUC2', all$from_id)
)


############# NEED TO CHANGE ENZYME DATA AND REMOVE HISTO FROM MISC ###########
###############################################################################
############# TRY AND ADD BUTS IN ############## ################ #############



filtered <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 5)

#all <- all[-remove,]




# library(sna)
# library(igraph)
#library(network)
# graphhy <- graph_from_data_frame(all)
# V(graphhy)[17]
# 
# g=delete.vertices(graphhy,which(degree(graphhy)==1))
# 
# V(g)
# neighborhood(g, order = 3, nodes = V(g)[14])
# get.edgelist(g)
# cliques(g)

# Plotting ----------------------------------------------------------------

set.seed(12812)

p <- ggplot(filtered, aes(from_id = from_id, to_id = to_id, label=from_id, color=type)) + 
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
           aes(color = type),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = FALSE, fontsize=3.5, singletons = FALSE,labelcolour="black",
           labelgeom = 'label') +
  theme_net()
p
#list(c('cell.pointpointrad', 3), c('niter', 1000))
# layout.par\$cell.pointpointrad
# layout.par\$cell.pointcellrad

#sna::gplot.layout 


#p$data
stuff <- ggplot_build(p)
newplot <- stuff$data[[1]]

#psandspears <- p$data



newplot <- newplot[!is.na(newplot$x),]
newplot$type <- NA
newplot$type[grep('Otu', newplot$from)] <- '16S'
newplot$type[grep('CD3-', newplot$from)] <- 'CD3-'
newplot$type[grep('CD3\\+', newplot$from)] <- 'CD3+'
newplot$type[grep('ate', newplot$from)] <- 'VFA'
newplot$type[grep('\\.rate', newplot$from)] <- 'Enzyme'

newplot$type[which(is.na(newplot$type))] <- 'mRNA'


CD25 <- newplot[grep('CD25\\+', newplot$from),]
FoxP3 <- newplot[grep('FoxP3\\+', newplot$from),]
CD4 <- newplot[grep('CD4\\+', newplot$from),]
CD8 <- newplot[grep('CD8a\\+', newplot$from),]

newplot$label <- gsub('CD3./(CD4./CD8a.)/FoxP3./CD25.', '\\1', newplot$label)
newplot$label <- gsub('(.*):Otu[0-9]+', '\\1', newplot$label)
labbes <- unique(newplot[,c(2,4,5,8,30)])

# labbes$label <- gsub('(.*)\\([0-9]+\\)', "italic('\\1')", labbes$label)
# labbes$label[grep('CD4', labbes$label)] <- gsub('(.*)', "bold('\\1')", labbes$label[grep('CD4', labbes$label)])
# labbes$label[grep('ate', labbes$label)] <- gsub('(.*)', "plain('\\1')", labbes$label[grep('ate', labbes$label)])
# labbes$label[grep('^[A-Z]', labbes$label)] <- gsub('(.*)', "plain('\\1')", labbes$label[grep('^[A-Z]', labbes$label)]) 
# 
# labbes$label <- gsub('/', ' ', labbes$label)
# labbes$label <- gsub('([-\\+])', '\\1phantom(0)', labbes$label)
# labbes$label <- gsub('(.*)', 'expression(\\1)', labbes$label)

newplot$type

pl <- ggplot(newplot, aes(color=type)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color='black') +
  geom_point(data=CD25, aes(x=x, y=y, color='CD25+'), color='red', size=12, show.legend = FALSE)+
  geom_point(aes(x=x, y=y), size=8, show.legend = FALSE) +
  geom_rect(data=CD8, aes(xmin=x-.001,ymin=y-.03, xmax=x+.001, ymax=y+.03), fill='purple', color='black', show.legend = FALSE)+
  geom_rect(data=CD4, aes(xmin=x-.03,ymin=y-.001, xmax=x+.03, ymax=y+.001), fill='yellow', color='yellow', show.legend = FALSE)+
  geom_point(data=FoxP3, aes(x=x, y=y), color='blue', size=4)+
  geom_label_repel(data=labbes,aes(x=x, y=y, label=label, fill=type),point.padding = unit(.75, 'lines'), fontface="bold", color='black',
                   alpha=.75, show.legend=FALSE, size=3, parse = FALSE) +
  
  
  
  # this is all for the legend
  annotate('point', x=.982, y=1.08, color="#F8766D", size=8) + annotate('text', x=1.027, y=1.08 , label='16S', size=3, color='black') +
  annotate('point', x=.982, y=1.04, color="#619CFF", size=8) + annotate('text', x=1.027, y=1.04 , label='mRNA', size=3, color='black') +
  annotate('point', x=.982, y=1.00, color="#00BFC4", size=8) + annotate('text', x=1.027, y=1.00, label='other', size=3, color='black') +
  annotate('point', x=.982, y=.96, color="#00BA38", size=8) + annotate('text', x=1.027, y=.96 , label='CD3+', size=3, color='black') +
  annotate('point', x=.982, y=.92, color="#B79F00", size=8) + annotate('text', x=1.027, y=.92 , label='CD3-', size=3, color='black') +
  annotate('point', x=.982, y=.88, color="grey", size=8) + annotate('text', x=1.027, y=.88 , label='CD4+', size=3, color='black') +
  annotate('point', x=.982, y=.84, color="grey", size=8) + annotate('text', x=1.027, y=.84 , label='CD8+', size=3, color='black') +
  annotate('point', x=.982, y=.80, color="red", size=11) + 
  annotate('point', x=.982, y=.80, color="grey", size=8) + annotate('text', x=1.027, y=.80 , label='CD25+', size=3, color='black') +
  annotate('point', x=.982, y=.76, color="grey", size=8) + annotate('text', x=1.027, y=.76 , label='FoxP3+', size=3, color='black') +
  annotate('rect', xmin=.959, xmax=1.005, ymin=.879, ymax=.881, color="yellow", fill = 'yellow') +
  annotate('rect', xmin=.981, xmax=.983, ymin=.818, ymax=.862, color="black", fill = 'black') +
  annotate('point', x=.982, y=.76, color="blue", size=4) +
  
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        aspect.ratio = 1) 


pl

plotdata <- pl$data

proc.time() - ptm

scale_fill_brewer(palette = 'set1')
brewer.pal
#############


########### BUT.R ############

setwd('~/FS2/but_amp/from_sci/')

library(vegan)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(reshape2)

##### functions #####

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

##### end functions #####

# but taxonomy?

x <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
x <- unlist(strsplit(x, split = ' '))


blastn <- read.table('butreps_blastn.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)


blast <- read.table('butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)

colnames(blast) <- x
colnames(blastn) <- x


blast$tit1 <- gsub('.* \\[(.*)\\]', '\\1', blast$salltitles)
blast$tit1 <- paste(blast$tit1, blast$pident, sep = ' ')
blast$tit1
blast$salltitles


shared2 <- read.table('but2.shared', header = TRUE, stringsAsFactors = FALSE)
hist(rowSums(shared2[,-c(1,2,3)]), breaks = 100)

sort(rowSums(shared2[,-c(1,2,3)]))
rowSums(shared2[,-c(1,2,3)]) > 2000

blast$otu <- colnames(shared2)[-c(1,2,3)]



#shared2$Group == meta$group


#### this chunk creates a named vector that you can use to swap out the generic otuxxx labels for the blast results

butswip <- blast$otu
names(butswip) <- blast$tit1
butswap <- names(butswip)
names(butswap) <- butswip

butswap[blast$otu]
#



shared2[rowSums(shared2[,-c(1,2,3)]) < 2000,]$Group
rownames(shared2) <- shared2$Group

otu2 <- shared2[,-c(1,2,3)]



otu2 <- otu2[rowSums(otu2) > 2000,]

otu2.r <- rrarefy(otu2, 2992)
otu2.r <- otu2.r[,colSums(otu2.r) > 0]

meta <- read.table('~/FS2/16s/V4/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

met <- data.frame(rownames(otu2.r))
colnames(met) <- 'group'

met1 <- read.table('met1.txt', header = TRUE, stringsAsFactors = FALSE)
colnames(met1)[1] <- 'group'

met <- merge(met, met1, by = 'group', all = TRUE)
meta <- merge(meta, met, by = 'group', all = TRUE)

meta$tissue[is.na(meta$tissue)] <- 'feces'
meta$day.x[is.na(meta$day.x)] <- meta$day.y[is.na(meta$day.x)]
meta$pig_num[is.na(meta$pig_num)] <- meta$pigNumbers[is.na(meta$pig_num)]
meta$treatment.x[is.na(meta$treatment.x)] <- meta$treatment.y[is.na(meta$treatment.x)]
meta$experiment[is.na(meta$experiment)] <- 'FS2'
meta <- meta[,1:6]
meta <- meta[meta$experiment =='FS2',]

colnames(meta) <- c('group', 'experiment', 'tissue', 'day', 'pig_num', 'treatment')

meta$tissue[is.na(meta$day)][2] <- 'cec_cont_DNA'

meta$day[is.na(meta$day)] <- 21
meta$pig_num[is.na(meta$pig_num)] <- 75
meta$treatment[is.na(meta$treatment)] <- 'RPS'
meta <- meta[-57,]

meta$treatment <- gsub('RS', 'RPS', meta$treatment)
meta$day <- gsub('D(*)', '\\1', meta$day)

meta <- meta[meta$group %in% rownames(otu2.r),]

otu2.r <- otu2.r[rownames(otu2.r) %in% meta$group,]

meta$tissue[meta$tissue == 'cec_cont_DNA'] <- 'cec_cont_RNA'

meta$design <- paste(meta$tissue, meta$day, meta$treatment, sep = '_')

meta$group == rownames(otu2.r)

meta <- meta[match(rownames(otu2.r), meta$group),]

meta$group == rownames(otu2.r)


# for correlation
but2 <- data.frame(otu2.r)
but2$group <- rownames(but2)

#write.table(but2, '~/FS2/correlation/but_corr/but_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)

##### Alpha Diversity #####


meta$invsimpson <- diversity(otu2.r, index = 'invsimpson')
meta$shannon <- diversity(otu2.r, index = 'shannon')

#filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
#  ggplot() + geom_boxplot(aes(day, shannon, fill = treatment)) +
#  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time')


### USE THIS FIG FOR FECES ALPHA DIV ###
filter(meta, day %in% c(0,12,15,19,21) & tissue == 'feces') %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time')

filter(meta, day == 21) %>%
  ggplot() + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of tissues')

###

## nifty little thing to do wilcoxon tests on alpha diversity at each day between the two treatments

fecal_alpha_wilcox_results <- filter(meta, tissue == 'feces') %>% group_by(day) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(day, Wilcox = wilco$p.value)



# wilcoxon tests for tissues
tissue_alpha_wilcox_results <- filter(meta, day ==21) %>% group_by(tissue) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(tissue, Wilcox = wilco$p.value)

###




####### Ordinations ############


### bray curtis calc ###
meta$design <- as.factor(meta$design )

FS2.bray2 <- vegdist(otu2.r, method = 'bray')

attributes(FS2.bray2)$Labels == meta$group # still good!

##### dispersion stuff #####

dispers <- betadisper(FS2.bray2, group = meta$design)
pdispers <- permutest(dispers, pairwise = TRUE)
pdispers$pairwise
dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group
metadisp <- merge(meta, dispersdf, by = 'group')

dispgroups <- summarize(group_by(metadisp, design), average_dist=mean(dispers.distances))

dispgroups <- unique(inner_join(dispgroups, meta))

metadisp %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) + 
  geom_boxplot() + scale_fill_brewer(palette = 'Dark2') + ylim(c(.15,.7)) + 
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")


# maybe this is useful?  polish a little bit
#doesnt work #
fecaldisptime <-  dispgroups %>% filter(tissue == 'feces' & day < 22)
fecaldisptime$day <- as.numeric(fecaldisptime$day)
ggplot(fecaldisptime, aes(x=day, y=average_dist, color = treatment)) + 
  geom_vline(xintercept = 12)+ geom_point(size = 2.5) +  geom_line() + ylim(c(.375,.543))+
  ggtitle('Community Variability (Dispersion)',
          subtitle = "Vegan's betadisper(): how much variability is there in a group's community structure?") + 
  ylab("Average distance to the group median") + scale_color_brewer(palette = "Dark2")


fecaldisptime$average_dist

metadisp %>% filter(day == 21) %>% 
  ggplot(aes(x=tissue, y=dispers.distances, group=design, fill=treatment)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Distance to group centroid')



##### NMDS #####

FS2.mds <- metaMDS(FS2.bray2, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds$stress

FS2.nmds <-as.data.frame(FS2.mds$points)
FS2.nmds$group <- rownames(FS2.nmds)



FS2.metanmds <- merge(meta, FS2.nmds, by = 'group')
FS2.metanmds$design <- factor(FS2.metanmds$design)
FS2.metanmds$treatment <- factor(FS2.metanmds$treatment)
FS2.metanmds$tissue <- factor(FS2.metanmds$tissue)


FS2.metanmds$treatmentXday <- paste(FS2.metanmds$treatment, FS2.metanmds$day, sep = ' day ')
FS2.metanmds$treatmentXday <- factor(FS2.metanmds$treatmentXday)

###

FS2.nmds$group == FS2.metanmds$group
#FS2.metanmds$group <- as.character(FS2.metanmds$group)
FS2.metanmds <- FS2.metanmds[match(FS2.nmds$group,FS2.metanmds$group),] # ok I think this works as long as $group isnt a factor...
FS2.nmds$group == FS2.metanmds$group

# this is where the elipse weirdness happens, sometimes the .mds and .metanmds have different ordering....


ord <- ordiellipse(FS2.mds,FS2.metanmds$design, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean <- aggregate(FS2.metanmds[,10:11], list(group=FS2.metanmds$design), mean)



df_ell <- data.frame()
for (d in levels(FS2.metanmds$design)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(FS2.metanmds[FS2.metanmds$design == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}

df_ell$treatment <- NA
df_ell$treatment[grep('RPS', df_ell$group)] <- 'RPS'
df_ell$treatment[grep('control', df_ell$group)] <- 'control'
df_ell$tissue <- gsub('(.*)_[0-9]+_.*', '\\1', df_ell$group)
df_ell$day <-gsub('(.*)_([0-9]+)_.*', '\\2', df_ell$group)


colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')
FS2.metanmds <- merge(FS2.metanmds, NMDS.mean , by='design')

df_ell$treatmentXday <- paste(df_ell$treatment, df_ell$day, sep = ' ')
FS2.metanmds$day <- factor(FS2.metanmds$day)


### D0 vs D21 ordination
p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day %in% c(0,21)), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.7, alpha = .75) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, color=treatment), size = .3) + 
  geom_path(data = subset(df_ell, day %in% c(0, 21) & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment), size=1) +
  ggtitle('Feces community similarity', subtitle = 'Day 0 vs Day 21\nPERMANOVA: control: F=10.8, p=0.0001; RPS: F=11.5, p=0.0001') + 
  scale_color_brewer(palette = 'Dark2') + labs(caption = 'Ordination stress = 0.221')

p


###D12

p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 12), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==12  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 12') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p

# d15
p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 15), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==15  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 15') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p

# d19
p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 19), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==19  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 19') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p

# d21
p <- ggplot(data=subset(FS2.metanmds, tissue == 'feces' & day == 21), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==21  & tissue == 'feces'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Feces Day 21') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p + geom_text((aes(label=pig_num)))

################ ###############

p <- ggplot(data=subset(FS2.metanmds, tissue == 'cec_cont_RNA' & day == 21), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==21  & tissue == 'cec_cont_RNA'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Cec_cont_RNA Day 21') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p + geom_text((aes(label=pig_num)))


### dont know if I need these daily ordinations above...  #####


groupxday <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxday <- filter(groupxday, day %in% c(0:22) & tissue == 'feces')


FS2.metanmds %>% filter(tissue == 'feces' & day %in% c(0:21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_path(aes(group=treatment, color=treatment), size = 2, alpha=.7) +
  geom_point(aes(color=treatment), size=7) +
  geom_text(data=groupxday, aes(x=groupX, y=groupY, label=day))+
  scale_color_brewer(palette = 'Dark2')  + ggtitle('Fecal community structure over time', subtitle = 'Treatment group centroids') + ylab('NMDS2') + xlab('NMDS1')

# the above one is nice for showing progression over time maybe

##### #####

groupxtiss <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxtiss <- filter(groupxtiss, day == 21)

FS2.metanmds %>% filter(day %in% c(21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_point(aes(color=treatment), size=8) +
  geom_path(data = subset(df_ell, day == 21), aes(x=NMDS1, y=NMDS2, group=group, color=treatment), size=1) +
  geom_text(data=groupxtiss, aes(x=groupX, y=groupY, label=tissue))+
  scale_color_brewer(palette = 'Dark2') + labs(x='NMDS 1', y= 'NMDS 2', caption = 'Ordination stress: 0.16') +
  ggtitle('Bray-curtis based NMDS ordination of tissues at day 21') 



# I really like this one.  I think I will keep it.
# do I need stats here?



##### PERMANOVA with Adonis #####

x <- otu2.r
factors <- meta$design


PW.Adonis <- pairwise.adonis(x,factors,sim.method="bray",p.adjust.m = "bonferroni")

# write.table(PW.Adonis,"Adonis-Results.csv",sep=",")

PW.Adonis$pairs


fecVSfec <- grep('feces_.* vs feces_.*' ,PW.Adonis$pairs)
colVScol <- grep('colon_.* vs colon_.*' ,PW.Adonis$pairs)
cecVScec <- grep('cecum_.* vs cecum_.*' ,PW.Adonis$pairs)
rVSr <- grep('RNA_.* vs RNA_.*' ,PW.Adonis$pairs)
ilVSil <- grep('ileum_.* vs ileum_.*' ,PW.Adonis$pairs)

good <- PW.Adonis[c(fecVSfec, colVScol, cecVScec, rVSr, ilVSil),]

# FS2_feces_0_control vs FS2_feces_21_control 10.7868015 0.30141843 0.0001000     0.0435
# FS2_feces_21_control vs FS2_feces_0_RPS 11.5118295 0.31529040 0.0001000     0.0435


fecsametime <- grep("feces_([0-9]+)_control vs feces_\\1_RPS", good$pairs)
ilsametime <- grep("ileum_([0-9]+)_control vs ileum_\\1_RPS", good$pairs)
colsametime <- grep("colon_([0-9]+)_control vs colon_\\1_RPS", good$pairs)
cecsametime <- grep("cecum_([0-9]+)_control vs cecum_\\1_RPS", good$pairs)

good2 <- good[c(fecsametime, ilsametime, colsametime, cecsametime),]
#write.table(good2,"Adonis-Results.csv",sep=",") 

pwadon <- read.table("Adonis-Results.csv",sep=",")
just_poop <- pwadon[grep('feces_.* vs feces_.*', pwadon$pairs),]
# poop_time <- just_poop[grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", just_poop$pairs),]
poop_time <- just_poop
poop_time$day <- c(0,21,12,15,19)

poop_time$p.value <- round(poop_time$p.value, 4)
poop_time$day
filter(poop_time, day < 22) %>%
  ggplot(aes(x=day, y=F.Model)) +
  geom_line()+geom_point(size=2)+ geom_vline(xintercept = 12, color='purple') + geom_text(aes(y = F.Model + .2, label = paste('p =', p.value))) +
  ggtitle('Dissimilarity of fecal microbiota over time', subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint, how different are the two diets at each timepoint? ') + 
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%')

# this is nice...maybe add dispersion info on this too?

###  relative to D0 for each treatment  ###

reltoD0 <- grep('feces_0_control vs feces_[0-9]+_control', good$pairs)
reltoD0 <- c(reltoD0, grep('feces_0_RPS vs feces_[0-9]+_RPS', good$pairs))

R2D0 <- good[reltoD0,]
R2D0 <- add_row(.data = R2D0, pairs = 'feces_0_control vs feces_0_control', F.Model = 0, p.value = 1, R2 = 0, p.adjusted = 1)
R2D0 <- add_row(.data = R2D0, pairs = 'feces_0_RPS vs feces_0_RPS', F.Model = 0, p.value = 1, R2 = 0, p.adjusted = 1)
R2D0$treatment <- ifelse(grepl('control', R2D0$pairs), yes = 'control', no = 'RPS')
R2D0$day <- gsub('.* vs feces_([0-9]+)_.*', '\\1', R2D0$pairs)
R2D0$day <- as.numeric(R2D0$day)



filter(R2D0, day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=F.Model, color = treatment)) + geom_point(size=2.5) + geom_line() + geom_vline(xintercept = 12)+
  ggtitle('Dissimilarity of fecal microbiota compared to Day 0 over time',
          subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint,\nhow different are the two diets from their original community structures at each timepoint? ') + 
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%') + scale_color_brewer(palette = 'Dark2')

####### Writing stuff for correlations #######

meta$sample <- paste(meta$tissue, meta$day, meta$pig_num, sep = '_')


but.feces <- otu2.r[meta$tissue == 'feces' & meta$day %in% c(0,21),]
but.cec_cont_RNA <- otu2.r[meta$tissue == 'cec_cont_RNA',]
but.cecum <-  otu2.r[meta$tissue == 'cecum',]
but.colon <-  otu2.r[meta$tissue == 'colon',]
but.ileum <-  otu2.r[meta$tissue == 'ileum',]

otu <- as.data.frame(otu2.r)
otu$Group <- rownames(otu)
#write.table(meta, '~/FS2/correlation/but_corr/but_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(otu, '~/FS2/correlation/but_corr/but.all.shared', col.names = TRUE, row.names = TRUE, quote = FALSE)

##### Deseq2 stuff need to switch this over to butstuff #####



otu <- import_mothur(mothur_shared_file = 'but2.shared')
# taxo <- import_mothur(mothur_constaxonomy_file = 'V4.final.taxonomy')
#meta <- read.table(file = '', sep = '\t', header = TRUE)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]



FS2 <- phyloseq(otu)
FS2 <- merge_phyloseq(FS2, phy_meta)                       # combines the metadata with this phyloseq object
#colnames(tax_table(FS2)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS2 <- subset_samples(FS2, experiment == 'FS2')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.
FS2 <- prune_taxa(taxa_sums(FS2) > 10, FS2)        # removes OTUs that occur less than 10 times globally

#FS2.genus <- tax_glom(FS2, taxrank = "Genus")

######## D0 #########

#D0 vs D21

FS2.feces <- subset_samples(FS2, tissue == 'feces')











############ ###############
FS2.D0 <- subset_samples(FS2, day == 0)

sample_sums(FS2.D0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)

rowSums(FS2.D0@otu_table)

FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)

FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")


res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0 = res.D0[which(res.D0$padj < .05), ]
sigtab.D0 <- as.data.frame(sigtab.D0)
# No sigdiff genera at D0
#sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS2.D0)[rownames(sigtab.D0), ], "matrix"))
format(sigtab.D0$padj, scientific = TRUE)
sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "RPS", "Control")


deseq.D0 <- ggplot(sigtab.D0, aes(x=reorder(rownames(sigtab.D0), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0), y=log2FoldChange+.6, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
deseq.D0



#####   D12 #####



FS2.D12 <- subset_samples(FS2, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")


res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D12 = res.D12[which(res.D12$padj < .05), ]
sigtab.D12 <- as.data.frame(sigtab.D12)


format(sigtab.D12$padj, scientific = TRUE)
sigtab.D12$newp <- format(round(sigtab.D12$padj, digits = 3), scientific = TRUE)
sigtab.D12$Treatment <- ifelse(sigtab.D12$log2FoldChange >=0, "RPS", "Control")

sigtab.D12$name <- butswap[rownames(sigtab.D12)]

deseq.D12 <- ggplot(sigtab.D12, aes(x=reorder(rownames(sigtab.D12), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D12), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D12 feces')+ coord_flip()
deseq.D12

########### D15 ###########

#####   D15 #####



FS2.D15 <- subset_samples(FS2, day == 15)

sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)

rowSums(FS2.D15@otu_table)

FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")


res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D15 = res.D15[which(res.D15$padj < .05), ]
sigtab.D15 <- as.data.frame(sigtab.D15)


format(sigtab.D15$padj, scientific = TRUE)
sigtab.D15$newp <- format(round(sigtab.D15$padj, digits = 3), scientific = TRUE)
sigtab.D15$Treatment <- ifelse(sigtab.D15$log2FoldChange >=0, "RPS", "Control")

sigtab.D15$name <- butswap[rownames(sigtab.D15)]

deseq.D15 <- ggplot(sigtab.D15, aes(x=reorder(rownames(sigtab.D15), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D15), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D15 feces')+ coord_flip()
deseq.D15

#####   D19 #####



FS2.D19 <- subset_samples(FS2, day == 19)

sample_sums(FS2.D19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

rowSums(FS2.D19@otu_table)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")


res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D19 = res.D19[which(res.D19$padj < .05), ]
sigtab.D19 <- as.data.frame(sigtab.D19)


format(sigtab.D19$padj, scientific = TRUE)
sigtab.D19$newp <- format(round(sigtab.D19$padj, digits = 3), scientific = TRUE)
sigtab.D19$Treatment <- ifelse(sigtab.D19$log2FoldChange >=0, "RPS", "Control")

sigtab.D19$name <- butswap[rownames(sigtab.D19)]

deseq.D19 <- ggplot(sigtab.D19, aes(x=reorder(rownames(sigtab.D19), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D19), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D19 feces')+ coord_flip()
deseq.D19
#####   D21 #####



FS2.D21 <- subset_samples(FS2, tissue == 'feces' & day != 0)

#FS2.D21 <- subset_samples(FS2, tissue == 'feces' & day == 21)

sample_sums(FS2.D21)
FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)

rowSums(FS2.D21@otu_table)

FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")


res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21 = res.D21[which(res.D21$padj < .05), ]
sigtab.D21 <- as.data.frame(sigtab.D21)


format(sigtab.D21$padj, scientific = TRUE)
sigtab.D21$newp <- format(round(sigtab.D21$padj, digits = 3), scientific = TRUE)
sigtab.D21$Treatment <- ifelse(sigtab.D21$log2FoldChange >=0, "RPS", "Control")

sigtab.D21$name <- butswap[rownames(sigtab.D21)]

deseq.D21 <- ggplot(sigtab.D21, aes(x=reorder(rownames(sigtab.D21), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21 feces')+ coord_flip()
deseq.D21

#####   D21.ileum #####



FS2.D21.ileum <- subset_samples(FS2, day == 21 & tissue =='ileum')

sample_sums(FS2.D21.ileum)
FS2.D21.ileum <- prune_taxa(taxa_sums(FS2.D21.ileum) > 1, FS2.D21.ileum)

rowSums(FS2.D21.ileum@otu_table)

FS2.D21.ileum.De <- phyloseq_to_deseq2(FS2.D21.ileum, ~ design)

FS2.D21.ileum.De <- DESeq(FS2.D21.ileum.De, test = "Wald", fitType = "parametric")


res.D21.ileum = results(FS2.D21.ileum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.ileum = res.D21.ileum[which(res.D21.ileum$padj < .05), ]
sigtab.D21.ileum <- as.data.frame(sigtab.D21.ileum)


format(sigtab.D21.ileum$padj, scientific = TRUE)
sigtab.D21.ileum$newp <- format(round(sigtab.D21.ileum$padj, digits = 3), scientific = TRUE)
sigtab.D21.ileum$Treatment <- ifelse(sigtab.D21.ileum$log2FoldChange >=0, "RPS", "Control")

sigtab.D21.ileum$name <- butswap[rownames(sigtab.D21.ileum)]

deseq.D21.ileum <- ggplot(sigtab.D21.ileum, aes(x=reorder(rownames(sigtab.D21.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.ileum), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.ileum feces')+ coord_flip()
deseq.D21.ileum

#####   D21.cecum #####



FS2.D21.cecum <- subset_samples(FS2, day == 21 & tissue =='cecum')

sample_sums(FS2.D21.cecum)
FS2.D21.cecum <- prune_taxa(taxa_sums(FS2.D21.cecum) > 1, FS2.D21.cecum)

rowSums(FS2.D21.cecum@otu_table)

FS2.D21.cecum.De <- phyloseq_to_deseq2(FS2.D21.cecum, ~ design)

FS2.D21.cecum.De <- DESeq(FS2.D21.cecum.De, test = "Wald", fitType = "parametric")


res.D21.cecum = results(FS2.D21.cecum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.cecum = res.D21.cecum[which(res.D21.cecum$padj < .1), ]
sigtab.D21.cecum <- as.data.frame(sigtab.D21.cecum)


format(sigtab.D21.cecum$padj, scientific = TRUE)
sigtab.D21.cecum$newp <- format(round(sigtab.D21.cecum$padj, digits = 3), scientific = TRUE)
sigtab.D21.cecum$Treatment <- ifelse(sigtab.D21.cecum$log2FoldChange >=0, "RPS", "Control")

sigtab.D21.cecum$name <- butswap[rownames(sigtab.D21.cecum)]

deseq.D21.cecum <- ggplot(sigtab.D21.cecum, aes(x=reorder(rownames(sigtab.D21.cecum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.cecum), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.cecum feces')+ coord_flip()
deseq.D21.cecum





#####   D21.colon #####



FS2.D21.colon <- subset_samples(FS2, day == 21 & tissue =='colon')

sample_sums(FS2.D21.colon)
FS2.D21.colon <- prune_taxa(taxa_sums(FS2.D21.colon) > 1, FS2.D21.colon)

rowSums(FS2.D21.colon@otu_table)

FS2.D21.colon.De <- phyloseq_to_deseq2(FS2.D21.colon, ~ design)

FS2.D21.colon.De <- DESeq(FS2.D21.colon.De, test = "Wald", fitType = "parametric")


res.D21.colon = results(FS2.D21.colon.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.colon = res.D21.colon[which(res.D21.colon$padj < .1), ]
sigtab.D21.colon <- as.data.frame(sigtab.D21.colon)


format(sigtab.D21.colon$padj, scientific = TRUE)
sigtab.D21.colon$newp <- format(round(sigtab.D21.colon$padj, digits = 3), scientific = TRUE)
sigtab.D21.colon$Treatment <- ifelse(sigtab.D21.colon$log2FoldChange >=0, "RPS", "Control")

sigtab.D21.colon$name <- butswap[rownames(sigtab.D21.colon)]

deseq.D21.colon <- ggplot(sigtab.D21.colon, aes(x=reorder(rownames(sigtab.D21.colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.colon), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.colon feces')+ coord_flip()
deseq.D21.colon

#############  cec cont rna  ###########


FS2.D21.cec_cont_RNA <- subset_samples(FS2, day == 21 & tissue =='cec_cont_RNA')

sample_sums(FS2.D21.cec_cont_RNA)
FS2.D21.cec_cont_RNA <- prune_taxa(taxa_sums(FS2.D21.cec_cont_RNA) > 1, FS2.D21.cec_cont_RNA)

rowSums(FS2.D21.cec_cont_RNA@otu_table)

FS2.D21.cec_cont_RNA.De <- phyloseq_to_deseq2(FS2.D21.cec_cont_RNA, ~ design)

FS2.D21.cec_cont_RNA.De <- DESeq(FS2.D21.cec_cont_RNA.De, test = "Wald", fitType = "parametric")


res.D21.cec_cont_RNA = results(FS2.D21.cec_cont_RNA.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.cec_cont_RNA = res.D21.cec_cont_RNA[which(res.D21.cec_cont_RNA$padj < .05), ]
sigtab.D21.cec_cont_RNA <- as.data.frame(sigtab.D21.cec_cont_RNA)


format(sigtab.D21.cec_cont_RNA$padj, scientific = TRUE)
sigtab.D21.cec_cont_RNA$newp <- format(round(sigtab.D21.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.D21.cec_cont_RNA$Treatment <- ifelse(sigtab.D21.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

sigtab.D21.cec_cont_RNA$name <- butswap[rownames(sigtab.D21.cec_cont_RNA)]

deseq.D21.cec_cont_RNA <- ggplot(sigtab.D21.cec_cont_RNA, aes(x=reorder(rownames(sigtab.D21.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.cec_cont_RNA), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.cec_cont_RNA feces')+ coord_flip()
deseq.D21.cec_cont_RNA

############ ALl together now  ##########

FS2.all <- prune_taxa(taxa_sums(FS2) > 1, FS2)

FS2.all <- phyloseq_to_deseq2(FS2.all, ~design)



FS2.all.De <- DESeq(FS2.all, test = "Wald", fitType = "parametric")


res.D21.cec_cont_RNA = results(FS2.D21.cec_cont_RNA.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D21.cec_cont_RNA = res.D21.cec_cont_RNA[which(res.D21.cec_cont_RNA$padj < .05), ]
sigtab.D21.cec_cont_RNA <- as.data.frame(sigtab.D21.cec_cont_RNA)


format(sigtab.D21.cec_cont_RNA$padj, scientific = TRUE)
sigtab.D21.cec_cont_RNA$newp <- format(round(sigtab.D21.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.D21.cec_cont_RNA$Treatment <- ifelse(sigtab.D21.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

sigtab.D21.cec_cont_RNA$name <- butswap[rownames(sigtab.D21.cec_cont_RNA)]

deseq.D21.cec_cont_RNA <- ggplot(sigtab.D21.cec_cont_RNA, aes(x=reorder(rownames(sigtab.D21.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D21.cec_cont_RNA), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.cec_cont_RNA feces')+ coord_flip()
deseq.D21.cec_cont_RNA
## ileum stuff ##

feces <- as.data.frame(otu2.r[meta$tissue == 'feces' & meta$day == 21 & meta$pig_num %in% c(67:80),])
feces.meta <- meta[meta$tissue == 'feces'& meta$day == 21 & meta$pig_num %in% c(67:80),]


ileum <- as.data.frame(otu2.r[meta$tissue == 'ileum',])
ileum.meta <- meta[meta$tissue == 'ileum',]



ilgtfec <- ileum > feces
bigilOTUs <- colSums(ilgtfec) > 11
bigilOTUs <- ileum[,bigilOTUs]
bigilOTUs$treatment <- c(rep('control', 7), rep('RPS', 7))

boxplot(bigilOTUs$Otu014~bigilOTUs$treatment)
boxplot(bigilOTUs$Otu035~bigilOTUs$treatment)
boxplot(bigilOTUs$Otu040~bigilOTUs$treatment)
boxplot(bigilOTUs$Otu106~bigilOTUs$treatment)
boxplot(bigilOTUs$Otu122~bigilOTUs$treatment)

wilcox.test(bigilOTUs$Otu014~bigilOTUs$treatment)
wilcox.test(bigilOTUs$Otu035~bigilOTUs$treatment)
wilcox.test(bigilOTUs$Otu040~bigilOTUs$treatment)
wilcox.test(bigilOTUs$Otu106~bigilOTUs$treatment)
wilcox.test(bigilOTUs$Otu122~bigilOTUs$treatment)

library(reshape2)


######################### allflow  ###################

setwd('~/FS2/Immune/')


meta <- read.table('~/FS2/16s/V4/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

CD3neg <- read.table(file = 'CD3-.txt.csv', header = TRUE, sep = '\t', as.is = TRUE, check.names = FALSE)

CD3pos <- read.table(file = 'CD3+.txt.csv', header = TRUE, sep = '\t', as.is = TRUE, check.names = FALSE)

all.flow <- merge(CD3neg, CD3pos, by = 'Sample')

all.flow <- all.flow[-(grep('-', all.flow$Sample)),] # gets rid of duplicate sample #80

colnames(all.flow)

all.flow$Sample <- gsub(' ', '_', all.flow$Sample)
all.flow$Sample <- gsub('Specimen_([0-9]+_[A-Za-z]+)_.*', '\\1', all.flow$Sample)

all.flow$tissue <- gsub('[0-9]+_([A-Za-z]+)', '\\1', all.flow$Sample)
all.flow$pig_num <- gsub('([0-9]+)_[A-Za-z]+', '\\1', all.flow$Sample)

row.names(all.flow) <- all.flow$Sample




all.flow$treatment[all.flow$pig_num %in% meta$pig_num[meta$treatment == 'control']] <- 'control'
all.flow$treatment[all.flow$pig_num %in% meta$pig_num[meta$treatment == 'RPS']] <- 'RPS'


#write.table(all.flow, file = 'FS2_all_flow.txt', sep = '\t', row.names = FALSE)




####################################### From Tcell Ecol Notebook ######################################



library(vegan)
library(tidyverse)
library(plotly)
library(ggrepel)

meta <- read.table('~/FS2/16s/V4/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)
flow_counts <- read.table('FS2_all_flow.txt', header = TRUE, check.names = FALSE, sep = '\t', as.is = TRUE, comment.char = '')
#meta <- meta[meta$day == 21 & meta$tissue == 'cecum',]
meta <- meta[]
#rownames(flow_counts) <- flow_counts$Pig

#flow_counts <- flow_counts[-15,-c(2,3)]
#colnames(flow_counts)[1] <- 'pig_num'
#rownames(flow_counts) <- flow_counts$pig_num

rowSums(flow_counts[,-c(1,34,35,36)]) # this tells us the total number of CD3+ cells each sample has. 
#Note how some pigs had many more events collected than others.  Pig numbers are above the number of events 

# this feature of the data is analogous to uneven sequencing depth we commonly encounter across our samples


# We can rarefy this data to ensure that each pig is evenly sampled
#flow_counts[,-c(1,34,35,36)] <- rrarefy(flow_counts[,-c(1,34,35,36)], min(rowSums(flow_counts[,-c(1,34,35,36)])))

rowSums(flow_counts[,-c(1,34,35,36)])  # this tells us the total number of CD3+ cells each sample has after rarefying




flow_counts[,-c(1,34,35,36)] <- flow_counts[,-c(1,34,35,36)] / rowSums(flow_counts[,-c(1,34,35,36)]) # converts to relative abundance


### alpha div  ####
meta <- flow_counts[,c(1,34,35,36)]
H <- diversity(flow_counts[,-c(1,34,35,36)])
J <- H/log(specnumber(flow_counts[,-c(1,34,35,36)]))
meta$evenness <- J
meta$Shannon <- H

ggplot(meta) + geom_boxplot(aes(tissue, evenness, fill = treatment))
ggplot(meta) + geom_boxplot(aes(tissue, Shannon, fill = treatment))

#### ####
barflow <- flow_counts
#barflow$treatment <- c(rep( c("control", "RPS"), each=7))
rownames(barflow) <- barflow$Sample
barflow <- barflow[,-1]
barflow <- barflow %>% gather(attribute, value, -c(pig_num, treatment, tissue))

barflow$pig_num <- factor(barflow$pig_num)
#barflow$type[grep("CTL", barflow$attribute)] <- 'CTL'
#barflow$type[grep("TH", barflow$attribute)] <- 'TH'
#barflow$type[grep("DN", barflow$attribute)] <- 'DN'
#barflow$type[grep("DP", barflow$attribute)] <- 'DP'

#barflow$FoxP3[grep("FoxP3+", barflow$attribute)] <- TRUE
#barflow$FoxP3[grep("FoxP3-", barflow$attribute)] <- FALSE
#barflow$CD25[grep("CD25+", barflow$attribute)] <- TRUE
#barflow$CD25[grep("CD25-", barflow$attribute)] <- FALSE

#barflow$attribute[order()]
#CTL <- grep('CTL', barflow$attribute)
#TH <- grep('TH', barflow$attribute)
#DP <- grep('DP', barflow$attribute)
#DN <- grep('DN', barflow$attribute)

#barflow$activation <- gsub('[A-Z]+_(FoxP3[+-]_CD25[+-])', '\\1', barflow$attribute)
library(ggmosaic)
ggplot(data=subset(x = barflow, tissue == 'cec'))+geom_mosaic(aes(x=product(pig_num), fill=attribute, weight=value), offset = 0.005) + facet_grid(~treatment) +
  labs(x="treatment", y= "Proportion of CD3+ cells")+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#ggplotly()
ggplot(data=subset(barflow, tissue == 'cec'), aes(x=pig_num, y=value, fill=attribute)) + geom_bar(stat = 'identity')+ labs(y = 'Proportion of live single cells')

ggplot(data=barflow, aes(x=pig_num, y=value, fill=attribute)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Proportion of live single cells') +
  facet_grid(~tissue) + geom_segment(x=7.5, y=0, xend=7.5, yend=1)

barflow[barflow$tissue =='PBMC'& barflow$treat == 'RPS',]$value <- barflow[barflow$tissue =='PBMC'& barflow$treat == 'RPS',]$value * 7/6 #B/C missing 1 obs in RPS treat
barflow$value <- barflow$value/7
ggplot(data=barflow, aes(x=treatment, y=value, fill=attribute)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Proportion of live single cells') +
  facet_grid(~tissue)


#gplotly()
#ggplot(data=barflow, aes(x=pig_num, y=value, fill=attribute)) + geom_bar(stat = 'identity')+ labs(y = 'Proportion of live single cells') + facet_grid(~tissue+treatment)





#meta$shannon <- diversity(flow_counts, index = 'shannon')
#meta$invsimp <- diversity(flow_counts, index = 'invsimpson')

#invsim.test <- wilcox.test(diversity(flow_counts, index = 'invsimpson')~meta$treatment)
#shann.test <- wilcox.test(diversity(flow_counts, index = 'shannon')~meta$treatment)


rownames(flow_counts) <- flow_counts$Sample
flow.bray <- vegdist(flow_counts[,-c(1,34,35,36)], method = 'bray')
flow.mds <- metaMDS(flow.bray, k = 2,trymax = 1000, autotransform = FALSE)
stressplot(flow.mds)
nmds.stress <- flow.mds$stress
nmds.stress 
flow.mds$sratmx

meta <- flow_counts[,c(1,34,35,36)]

flow.nmds <-as.data.frame(flow.mds$points)
flow.nmds$Sample <- rownames(flow.nmds)

flow.meta.nmds <- merge(meta, flow.nmds, by = 'Sample')
flow.meta.nmds$treatment <- factor(flow.meta.nmds$treatment)
flow.meta.nmds$tissueXtreat <- paste(flow.meta.nmds$tissue, flow.meta.nmds$treatment, sep = '_')

ord <- ordiellipse(flow.mds,flow.meta.nmds$tissueXtreat, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.flow <- aggregate(flow.meta.nmds[,5:6], list(group=flow.meta.nmds$tissueXtreat), mean)

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


flow.meta.nmds$tissueXtreat <- factor(flow.meta.nmds$tissueXtreat)

df_ell <- data.frame()
for (d in levels(flow.meta.nmds$tissueXtreat)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(flow.meta.nmds[flow.meta.nmds$tissueXtreat == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}



flow.meta.nmds$centroidX <- NA
flow.meta.nmds$centroidY <- NA


flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'cec_control',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'cec_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'ln_control',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'ln_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'PBMC_control',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'PBMC_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'cec_RPS',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'cec_RPS']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'ln_RPS',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'ln_RPS']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'PBMC_RPS',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'PBMC_RPS']


flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'cec_control',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'cec_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'ln_control',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'ln_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'PBMC_control',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'PBMC_control']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'cec_RPS',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'cec_RPS']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'ln_RPS',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'ln_RPS']
flow.meta.nmds[flow.meta.nmds$tissueXtreat == 'PBMC_RPS',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'PBMC_RPS']




#adonis.feces <- adonis(flow_counts~meta$treatment, permutations = 999999)
#adon.pval <- adonis.feces$aov.tab$`Pr(>F)`[1]
#adonis.feces$aov.tab

ord.fit <- envfit(flow.mds, flow_counts, permutations = 99999)
spp.scrs <- as.data.frame(scores(ord.fit, display = 'vectors'))

spp.scrs <- spp.scrs[which(ord.fit$vectors$pvals < 0.01),]
spp.scrs <- spp.scrs/7
spp.scrs <- cbind(spp.scrs, cell_type = rownames(spp.scrs))
colnames(spp.scrs) <- c('MDS1', 'MDS2', 'cell_type')


p.floword <- ggplot(flow.meta.nmds, aes(x=MDS1, y = MDS2)) +
  #geom_text_repel(data=spp.scrs, aes(MDS1, MDS2, label = cell_type), size=3, alpha=.7)+
  #geom_segment(data = spp.scrs, aes(x=0, xend=spp.scrs$MDS1, y=0, yend=spp.scrs$MDS2), alpha=.5)+
  geom_point(data=flow.meta.nmds, aes(color=tissueXtreat), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Cecal T-cell community similarity: NMDS ordination using Bray-Curtis distances', subtitle = 'Treatments separated within Tissues') + 
  geom_segment(data = flow.meta.nmds,
               aes(x=flow.meta.nmds$MDS1,
                   xend=flow.meta.nmds$centroidX,
                   y=flow.meta.nmds$MDS2,
                   yend=flow.meta.nmds$centroidY,
                   color=tissueXtreat)) + scale_color_brewer(palette="Dark2")

p.floword

############################## luimping treatments ############################



flow.meta.nmds$tissue <- factor(flow.meta.nmds$tissue)



ord <- ordiellipse(flow.mds,flow.meta.nmds$tissue, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.flow <- aggregate(flow.meta.nmds[,5:6], list(group=flow.meta.nmds$tissue), mean)

df_ell <- data.frame()

for (d in levels(flow.meta.nmds$tissue)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(flow.meta.nmds[flow.meta.nmds$tissue == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}



flow.meta.nmds$centroidX <- NA
flow.meta.nmds$centroidY <- NA


flow.meta.nmds[flow.meta.nmds$tissue == 'cec',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'cec']
flow.meta.nmds[flow.meta.nmds$tissue == 'ln',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'ln']
flow.meta.nmds[flow.meta.nmds$tissue == 'PBMC',]$centroidX <- NMDS.mean.flow$MDS1[NMDS.mean.flow$group == 'PBMC']



flow.meta.nmds[flow.meta.nmds$tissue == 'cec',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'cec']
flow.meta.nmds[flow.meta.nmds$tissue == 'ln',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'ln']
flow.meta.nmds[flow.meta.nmds$tissue == 'PBMC',]$centroidY <- NMDS.mean.flow$MDS2[NMDS.mean.flow$group == 'PBMC']



#adonis.feces <- adonis(flow_counts~meta$treatment, permutations = 999999)
#adon.pval <- adonis.feces$aov.tab$`Pr(>F)`[1]
#adonis.feces$aov.tab

ord.fit <- envfit(flow.mds, flow_counts, permutations = 99999)
spp.scrs <- as.data.frame(scores(ord.fit, display = 'vectors'))

spp.scrs <- spp.scrs[which(ord.fit$vectors$pvals < 0.01),]
spp.scrs <- spp.scrs/7
spp.scrs <- cbind(spp.scrs, cell_type = rownames(spp.scrs))
colnames(spp.scrs) <- c('MDS1', 'MDS2', 'cell_type')





p.floword <- ggplot(flow.meta.nmds, aes(x=MDS1, y = MDS2)) +
  #geom_text_repel(data=spp.scrs, aes(MDS1, MDS2, label = cell_type), size=3, alpha=.7)+
  #geom_segment(data = spp.scrs, aes(x=0, xend=spp.scrs$MDS1, y=0, yend=spp.scrs$MDS2), alpha=.5)+
  geom_point(data=flow.meta.nmds, aes(color=tissue), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('T-cell community similarity: NMDS ordination using Bray-Curtis distances', subtitle = 'Treatments lumped by tissue') + 
  #geom_text(aes(label=pig_num))+
  geom_segment(data = flow.meta.nmds,
               aes(x=flow.meta.nmds$MDS1,
                   xend=flow.meta.nmds$centroidX,
                   y=flow.meta.nmds$MDS2,
                   yend=flow.meta.nmds$centroidY,
                   color=tissue)) + scale_color_brewer(palette='Set1')

p.floword
#############  Just cecal  ##############

flow_counts <- flow_counts[flow_counts$tissue == 'cec',-c(1,34,35,36)]
flow.bray <- vegdist(flow_counts, method = 'bray')
flow.mds <- metaMDS(flow.bray, k = 2,trymax = 1000, autotransform = FALSE)
stressplot(flow.mds)
nmds.stress <- flow.mds$stress
flow.nmds <-as.data.frame(flow.mds$points)
flow.nmds$pig_num <- rownames(flow.nmds)
flow.nmds$pig_num <- gsub('([0-9]+)_[a-z]+', '\\1', flow.nmds$pig_num)

meta2 <- meta[meta$tissue == 'cec',]

flow.meta.nmds <- merge(meta2, flow.nmds, by = 'pig_num')
flow.meta.nmds$treatment <- factor(flow.meta.nmds$treatment)


ord <- ordiellipse(flow.mds,flow.meta.nmds$treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')
NMDS.mean.flow <- aggregate(flow.meta.nmds[,5:6], list(group=flow.meta.nmds$treatment), mean)

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for (d in levels(flow.meta.nmds$treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(flow.meta.nmds[flow.meta.nmds$treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),group=d))
}



flow.meta.nmds$centroidX <- NA
flow.meta.nmds$centroidY <- NA


flow.meta.nmds[flow.meta.nmds$treatment == 'control',]$centroidX <- NMDS.mean.flow$MDS1[1]
flow.meta.nmds[flow.meta.nmds$treatment == 'RPS',]$centroidX <- NMDS.mean.flow$MDS1[2]

flow.meta.nmds[flow.meta.nmds$treatment == 'control',]$centroidY <- NMDS.mean.flow$MDS2[1]
flow.meta.nmds[flow.meta.nmds$treatment == 'RPS',]$centroidY <- NMDS.mean.flow$MDS2[2]


adonis.feces <- adonis(flow_counts~meta2$treatment, permutations = 9999)
adon.pval <- adonis.feces$aov.tab$`Pr(>F)`[1]
adonis.feces$aov.tab

bdisp.tcells <- betadisper(flow.bray, meta2$treatment)
bdisp.tcells

ord.fit <- envfit(flow.mds, flow_counts, permutations = 99999)
spp.scrs <- as.data.frame(scores(ord.fit, display = 'vectors'))
spp.scrs <- spp.scrs[which(ord.fit$vectors$pvals < 0.01),]
spp.scrs <- spp.scrs/7
spp.scrs <- cbind(spp.scrs, cell_type = rownames(spp.scrs))
colnames(spp.scrs) <- c('MDS1', 'MDS2', 'cell_type')


p.floword <- ggplot(flow.meta.nmds, aes(x=MDS1, y = MDS2)) +
  geom_segment(data = spp.scrs, aes(x=0, xend=spp.scrs$MDS1, y=0, yend=spp.scrs$MDS2), alpha=.5)+
  geom_point(data=flow.meta.nmds, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Cecal T-cell community similarity: NMDS ordination using Bray-Curtis distances', subtitle = paste('PERMANOVA p-value = ', adon.pval)) + 
  geom_segment(data = flow.meta.nmds,
               aes(x=flow.meta.nmds$MDS1,
                   xend=flow.meta.nmds$centroidX,
                   y=flow.meta.nmds$MDS2,
                   yend=flow.meta.nmds$centroidY,
                   color=treatment)) + 
  geom_text_repel(data=spp.scrs, aes(MDS1, MDS2, label = cell_type), size=4, alpha=.7)+
  scale_color_brewer(palette="Dark2")

p.floword

####################################  boxplots  ########################################




fin <- data.frame()

# this for loop does wilcoxon tests for each cell type and saves the results into a new column that I can use to label??
barflow$attribute <- factor(barflow$attribute)

barflow$tissue <- factor(barflow$tissue)


barflowtoo <- barflow

barflowtoo$tisXatt <- paste(barflowtoo$tissue, barflowtoo$attribute, sep = '_')


barflowtoo$tisXatt <- factor(barflowtoo$tisXatt)







#for (tiss in levels(barflow$tissue)){
#  temp1 <- barflow[which(barflow$tissue == tiss),]
#  for (lev in levels(temp1$attribute)){
#    temp <- temp1[which(temp1$attribute == lev),]
#    twilc <- wilcox.test(temp$value ~ temp$treatment)
#    temp$pvalue <- twilc$p.value
#    fin <- rbind(fin,temp)
#  }
#}
#




fin <- data.frame()

for (lev in levels(barflowtoo$tisXatt)){
  temp <- barflowtoo[which(barflowtoo$tisXatt == lev),]
  twilc <- wilcox.test(temp$value ~ temp$treatment)
  temp$pvalue <- twilc$p.value
  fin <- rbind(fin,temp)
}


fin$pvalue2 <- round(fin$pvalue, 4)
fin$pvalue2 <- factor(paste('wilcox pvalue = ', fin$pvalue2, sep = ''))


fin$value <- fin$value * 100

p <- ggplot(subset(fin, tissue == 'cec'), aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+
  facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette="Dark2")+
  labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))


p


