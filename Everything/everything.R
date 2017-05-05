setwd('~/FS2/Everything/')

ptm <- proc.time()

set.seed(7777777)

library(lubridate)
library(vegan)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(reshape2)
library(ccrepe)
library(geomnet)
library(Hmisc)
library(ggrepel) # just put this up here...

system('ls')

################# FUNCTIIONS ###############


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

ccrepe_to_ggnet <- function(ccrepe.obj, pcut=0.05, spearcut=0.6){
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
      #print(x)
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

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


pairwise.adonis <- function(x,factors, sim.method, p.adjust.m){
  #this function taken from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
  
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

########## 16S  #########

# reading in the data #

tax <- read.table('V4.final.taxonomy', header = TRUE, stringsAsFactors = FALSE)
shared <- read.table('V4.final.shared', header=TRUE, stringsAsFactors = FALSE)
rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

# NTC stuff #

ntc <- shared[rownames(shared) == 'NTC',]
rowSums(ntc) # 71 sequences in my NTC
ntc[,which(ntc >1)]
# this could be due to index read errors or low level contamination
# end NTC #

shared <- shared[rownames(shared) != 'NTC',] # removes NTC from otu table

# building metadata from sample names, there were a couple of other experiments on this
# MiSeq run, the RPS experiment is FS2 #

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

# writing out meta and shared for correlation, these have an extra 'group' column

meta_for_corr <- filter(meta, experiment == 'FS2')
meta_for_corr$sample <- paste(meta_for_corr$tissue, meta_for_corr$day, meta_for_corr$pig_num, sep = '_')
#write.table(meta_for_corr, '~/FS2/correlation/but_corr/16S_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
shared2 <- shared[meta$experiment == 'FS2',]
shared2$group <- rownames(shared2)
#write.table(shared2, '~/FS2/correlation/but_corr/16S_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(design, file = 'V4.design', quote = FALSE, sep = '\t', row.names = FALSE)
#write.table(FS2.accnos, file = 'FS2.accnos', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
#write.table(meta, file = 'V4.metadata.txt', quote = FALSE, sep = '\t', row.names = FALSE)
# I wrote these out so I can load them later or in other scripts so I dont have to run this chunk constantly

# alpha diversity calcs #

FS2 <- shared[meta$experiment == "FS2",]
FS2.meta <- meta[meta$experiment == 'FS2',]
FS2 <- FS2[FS2.meta$tissue != 'cec_cont_DNA',]
# my cec_cont DNA samples failed miserably.  I think it has to do with residual RNAlater salts inhibiting the PCR reactions
FS2.meta <- FS2.meta[FS2.meta$tissue != 'cec_cont_DNA',]  

#write.table(FS2, '~/FS2/correlation/but_corr/FS216S.shared', quote = FALSE, sep = '\t', row.names = TRUE)

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

# USE THIS FIG FOR FECES ALPHA DIV #
fectime.alpha.16S <- filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time') +
  annotate('segment', x=4.7, xend=5.2, y=30, yend=30, size=.3) + annotate('text', x=4.9, y=31, label='* p=0.029', size=4)

#fectime.alpha.16S

# nifty little thing to do wilcoxon tests on alpha diversity at each day between the two treatments

fecal_alpha_wilcox_results <- FecalTime.meta %>% group_by(day) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(day, Wilcox = wilco$p.value)
#

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
p.tissue.alpha.16S <- ggplot(Tissue21.meta) + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') +
  ggtitle('Alpha diversity of tissues') + annotate('text', x=4, y=36.5, label= '* p = 0.022 *')+
  annotate('segment', x=3.8, xend=4.2, y=35, yend=35)

#p.tissue.alpha.16S

FS2.meta$design <- as.factor(FS2.meta$design) # do I need this?

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

min(rowSums(cec_cont_RNA.shared)) # 21805 seqs per sample in cec cont RNA
cec_cont_RNA.shared <- as.data.frame(rrarefy(cec_cont_RNA.shared ,min(rowSums(cec_cont_RNA.shared))))

#boxplot(diversity(cec_cont_RNA.shared, index = 'invsimpson')~cec_cont_RNA.meta$treatment)
#wilcox.test((diversity(cec_cont_RNA.shared, index = 'invsimpson')~cec_cont_RNA.meta$treatment))

# START ORDS #
# these ordinations use the OTU tables generated above, they only contain one type of tissue at d21 and are rarefied to the min number of reads in that tissue category

# Feces Ordination #

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


# Statistical Tests #

adonis.feces <- adonis(feces.shared~feces.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.feces, group = FS2.metanmds.feces$treatment)
permutest(BDISP, permutations = 9999)

#

p.feces <- ggplot(FS2.metanmds.feces, aes(x=MDS1, y = MDS2)) +
  geom_point(data=FS2.metanmds.feces, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Fecal community similarity (Day 21): NMDS ordination using Bray-Curtis distances',
          subtitle = 'Group similarity (PERMANOVA): p = 1e-05, F = 5.66, R2 = 0.18\nGroup dispersion (PERMDISP2): p = 0.0015, F = 13.12') + 
  geom_segment(data = FS2.metanmds.feces,
               aes(x=FS2.metanmds.feces$MDS1,
                   xend=FS2.metanmds.feces$centroidX,
                   y=FS2.metanmds.feces$MDS2,
                   yend=FS2.metanmds.feces$centroidY,
                   color=treatment)) + scale_color_brewer(palette="Dark2") + 
  labs(caption = 'Ordination stress = 0.14')

p.feces

# ileum ordination #

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

# Statistical Tests #

BDISP <- betadisper(FS2.bray.ileum, group = FS2.metanmds.ileum$treatment)
permutest(BDISP, permutations = 9999)
adonis(ileum.shared~ileum.meta$treatment, permutations = 9999)

# plot #

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


p.ileum

# cecum ordination #

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

# Statistical Tests #

adonis(cecum.shared~cecum.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.cecum, group = FS2.metanmds.cecum$treatment)
permutest(BDISP, permutations = 9999)

# plot #

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
p.cecum

# colon ord #

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

# Statistical tests #

adonis(colon.shared~colon.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.colon, group = FS2.metanmds.colon$treatment)
permutest(BDISP, permutations = 9999)

# plot #

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

p.colon

# cec_cont_RNA ordination #

FS2.bray.cec_cont_RNA <- vegdist(cec_cont_RNA.shared, method = 'bray')
FS2.mds.cec_cont_RNA <- metaMDS(FS2.bray.cec_cont_RNA, k = 2,trymax = 1000, autotransform = FALSE)
FS2.mds.cec_cont_RNA$stress

FS2.nmds.cec_cont_RNA <-as.data.frame(FS2.mds.cec_cont_RNA$points)
FS2.nmds.cec_cont_RNA$group <- rownames(FS2.nmds.cec_cont_RNA)

FS2.metanmds.cec_cont_RNA <- merge(meta, FS2.nmds.cec_cont_RNA, by = 'group')
FS2.metanmds.cec_cont_RNA$design <- factor(FS2.metanmds.cec_cont_RNA$design)
FS2.metanmds.cec_cont_RNA$treatment <- factor(FS2.metanmds.cec_cont_RNA$treatment)


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


# Statistical Tests #

adonis(cec_cont_RNA.shared~cec_cont_RNA.meta$treatment, permutations = 9999)
BDISP <- betadisper(FS2.bray.cec_cont_RNA, group = FS2.metanmds.cec_cont_RNA$treatment)
permutest(BDISP, permutations = 9999)

# plot #

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

p.cec_cont_RNA

# Deseq2 differential abundance #

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

# D0 #

FS2.D0 <- subset_samples(FS2.genus, day == 0)

sample_sums(FS2.D0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)

rowSums(FS2.D0@otu_table)

FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)

FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")


res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0 = res.D0[which(res.D0$padj < .05), ]


# No sigdiff genera at D0

# D12 #

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
# deseq.D12

# 4 genera enriched in RPS pigs, none in control, not a good figure....

# D15 #


FS2.D15 <- subset_samples(FS2.genus, day == 15)
sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)


FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")

# NO DIFF ABUND GENERA AT D15

# D19 #

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


deseq.D19 <- ggplot(sigtab.D19, aes(x=reorder(rownames(sigtab.D19), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D19), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: feces', subtitle = 'Day 19')+ coord_flip()
deseq.D19

# not a bad fig, but maybe better for supplement, and use ordinations to show timepoints prior to D21

# D21 #

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
deseq.D21

# lumping all post D0 timepoints #

FS2.all.feces <- subset_samples(FS2.genus, day %in% c(12,15,19,21) & tissue == 'feces')



FS2.all.feces <- prune_taxa(taxa_sums(FS2.all.feces) > 1, FS2.all.feces)
FS2.all.feces.De <- phyloseq_to_deseq2(FS2.all.feces, ~ treatment)

FS2.all.feces.De <- DESeq(FS2.all.feces.De, test = "Wald", fitType = "parametric")


res.feces.all = results(FS2.all.feces.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')

sigtab.feces.all = res.feces.all[which(res.feces.all$padj < .05), ]
sigtab.feces.all = cbind(as(sigtab.feces.all, "data.frame"), as(tax_table(FS2.all.feces)[rownames(sigtab.feces.all), ], "matrix"))
format(sigtab.feces.all$padj, scientific = TRUE)
sigtab.feces.all$newp <- format(round(sigtab.feces.all$padj, digits = 3), scientific = TRUE)
sigtab.feces.all$Treatment <- ifelse(sigtab.feces.all$log2FoldChange >=0, "RPS", "Control")


deseq.feces.all <- ggplot(sigtab.feces.all, aes(x=reorder(rownames(sigtab.feces.all), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.feces.all), y=0, label = paste(Family, Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(), 
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Differentially abundant genera: feces', subtitle = 'Day 12, 15, 19, 21 lumped')+ coord_flip()
deseq.feces.all

# This is all D21 tissue stuff #

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

#

sigtab.ileum = res.ileum[which(res.ileum$padj < 0.05), ]
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
deseq.ileum

sigtab.cecum = res.cecum[which(res.cecum$padj < 0.05), ]
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
deseq.cecum


sigtab.colon = res.colon[which(res.colon$padj < 0.05), ]
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
deseq.colon


sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < 0.05), ]
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
deseq.cec_cont_RNA

# Now doing everything over just at OTU level this time   #

# this is for supplementary table stuff and also for limiting network visualizations to only differentially abundant features #

# MIGHT NOT NEED TO READ THESE IN AGAIN #
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

FS2.D0 <- subset_samples(FS2, day == 0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)
FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)
FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")

res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0 = res.D0[which(res.D0$padj < 0.1), ]

sigtab.D0 = cbind(as(sigtab.D0, "data.frame"), as(tax_table(FS2.D0)[rownames(sigtab.D0), ], "matrix"))
format(sigtab.D0$padj, scientific = TRUE)
sigtab.D0$newp <- format(round(sigtab.D0$padj, digits = 3), scientific = TRUE)
sigtab.D0$Treatment <- ifelse(sigtab.D0$log2FoldChange >=0, "RPS", "Control")
sigtab.D0$tissue <- 'D0_feces'
# write.table(sigtab.D0, file = './OTUs/D0_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D12 #

FS2.D12 <- subset_samples(FS2, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")

res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D12 = res.D12[which(res.D12$padj < 0.1), ]
sigtab.D12 = cbind(as(sigtab.D12, "data.frame"), as(tax_table(FS2.D12)[rownames(sigtab.D12), ], "matrix"))
format(sigtab.D12$padj, scientific = TRUE)
sigtab.D12$newp <- format(round(sigtab.D12$padj, digits = 3), scientific = TRUE)
sigtab.D12$Treatment <- ifelse(sigtab.D12$log2FoldChange >=0, "RPS", "Control")
sigtab.D12$tissue <- 'D12_feces'
sigtab.D12$otu <- rownames(sigtab.D12)

#write.table(sigtab.D12, file = './OTUs/D12_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D15 #

FS2.D15 <- subset_samples(FS2, day == 15)
sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)

FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")

res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D15 = res.D15[which(res.D15$padj < 0.1), ]
sigtab.D15 = cbind(as(sigtab.D15, "data.frame"), as(tax_table(FS2.D15)[rownames(sigtab.D15), ], "matrix"))
format(sigtab.D15$padj, scientific = TRUE)
sigtab.D15$newp <- format(round(sigtab.D15$padj, digits = 3), scientific = TRUE)
sigtab.D15$Treatment <- ifelse(sigtab.D15$log2FoldChange >=0, "RPS", "Control")
sigtab.D15$tissue <- 'D15_feces'
sigtab.D15$otu <- rownames(sigtab.D15)

#write.table(sigtab.D15, file = './OTUs/D15_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D19 #

FS2.D19 <- subset_samples(FS2, day == 19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")

res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D19 = res.D19[which(res.D19$padj < 0.1), ]
sigtab.D19 = cbind(as(sigtab.D19, "data.frame"), as(tax_table(FS2.D19)[rownames(sigtab.D19), ], "matrix"))
format(sigtab.D19$padj, scientific = TRUE)
sigtab.D19$newp <- format(round(sigtab.D19$padj, digits = 3), scientific = TRUE)
sigtab.D19$Treatment <- ifelse(sigtab.D19$log2FoldChange >=0, "RPS", "Control")
sigtab.D19$tissue <- 'D19_feces'
sigtab.D19$otu <- rownames(sigtab.D19)

#write.table(sigtab.D19, file = './OTUs/D19_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D21 #

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
sigtab.D21$tissue <- 'D21_feces'
sigtab.D21$otu <- rownames(sigtab.D21)

#write.table(sigtab.D21, file = './OTUs/D21_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# Day 12,15,19,21 lumped #

FS2.D21 <- subset_samples(FS2, day %in% c(12,15,19,21) & tissue == 'feces')

FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)
FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")

res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.lumped_feces = res.D21[which(res.D21$padj < 0.05), ]
sigtab.lumped_feces = cbind(as(sigtab.lumped_feces, "data.frame"), as(tax_table(FS2.D21)[rownames(sigtab.lumped_feces), ], "matrix"))
format(sigtab.lumped_feces$padj, scientific = TRUE)
sigtab.lumped_feces$newp <- format(round(sigtab.lumped_feces$padj, digits = 3), scientific = TRUE)
sigtab.lumped_feces$Treatment <- ifelse(sigtab.lumped_feces$log2FoldChange >=0, "RPS", "Control")
sigtab.lumped_feces$tissue <- 'lumped_feces'
sigtab.lumped_feces$otu <- rownames(sigtab.lumped_feces)

#write.table(sigtab.lumped_feces, file = './OTUs/feceslumped_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# This is all D21 tissue stuff #

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

sigtab.ileum = res.ileum[which(res.ileum$padj < 0.1), ]
sigtab.ileum = cbind(as(sigtab.ileum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.ileum), ], "matrix"))
format(sigtab.ileum$padj, scientific = TRUE)
sigtab.ileum$newp <- format(round(sigtab.ileum$padj, digits = 3), scientific = TRUE)
sigtab.ileum$Treatment <- ifelse(sigtab.ileum$log2FoldChange >=0, "RPS", "Control")
sigtab.ileum$tissue <- 'ileum'
sigtab.ileum$otu <- rownames(sigtab.ileum)

sigtab.cecum = res.cecum[which(res.cecum$padj < 0.1), ]
sigtab.cecum = cbind(as(sigtab.cecum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cecum), ], "matrix"))
format(sigtab.cecum$padj, scientific = TRUE)
sigtab.cecum$newp <- format(round(sigtab.cecum$padj, digits = 3), scientific = TRUE)
sigtab.cecum$Treatment <- ifelse(sigtab.cecum$log2FoldChange >=0, "RPS", "Control")
sigtab.cecum$tissue <- 'cecum'
sigtab.cecum$otu <- rownames(sigtab.cecum)

sigtab.colon = res.colon[which(res.colon$padj < 0.1), ]
sigtab.colon = cbind(as(sigtab.colon, "data.frame"), as(tax_table(FS2)[rownames(sigtab.colon), ], "matrix"))
format(sigtab.colon$padj, scientific = TRUE)
sigtab.colon$newp <- format(round(sigtab.colon$padj, digits = 3), scientific = TRUE)
sigtab.colon$Treatment <- ifelse(sigtab.colon$log2FoldChange >=0, "RPS", "Control")
sigtab.colon$tissue <- 'colon'
sigtab.colon$otu <- rownames(sigtab.colon)

sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < 0.1), ]
sigtab.cec_cont_RNA = cbind(as(sigtab.cec_cont_RNA, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cec_cont_RNA), ], "matrix"))
format(sigtab.cec_cont_RNA$padj, scientific = TRUE)
sigtab.cec_cont_RNA$newp <- format(round(sigtab.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.cec_cont_RNA$Treatment <- ifelse(sigtab.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")
sigtab.cec_cont_RNA$tissue <- 'cec_cont_RNA'
sigtab.cec_cont_RNA$otu <- rownames(sigtab.cec_cont_RNA)

# write.table(sigtab.cec_cont_RNA, file = './OTUs/cec_cont_RNA_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(sigtab.colon, file = './OTUs/colon_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(sigtab.cecum, file = './OTUs/cecum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(sigtab.ileum, file = './OTUs/ileum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

dif_ab_16S <- rbind(sigtab.D12,
              sigtab.D15,
              sigtab.D19,
              sigtab.D21,
              sigtab.lumped_feces,
              sigtab.cec_cont_RNA,
              sigtab.colon,
              sigtab.cecum,
              sigtab.ileum)

write.table(dif_ab_16S, 'dif_ab_16s.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# end OTU level stuff  #
# Below here everything is rarrefied to 4200 sequences #

FS2.4200 <- read.table("V4.FS2.4200.shared", header = TRUE)
rownames(FS2.4200) <- FS2.4200$Group
meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)

meta <- meta[meta$group %in% rownames(FS2.4200),]
meta <- meta[meta$design != "FS2_cec_cont_DNA_21_RPS",]

FS2.4200 <- FS2.4200[rownames(FS2.4200) %in% meta$group,]

FS2.4200 <- FS2.4200[,-c(1,2,3)]

# OTU 87 STUFF #

tax <- extract_mothur_tax('V4.final.taxonomy')

timepo <- FS2.4200[,-1]/4200

timepo$group <- rownames(timepo)
timepo <- merge(meta, timepo, by = 'group')
timepo$day <- as.numeric(timepo$day)


timepo2 <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(tissue == 'feces', day <  22)
timepo2$otu2 <- factor(timepo2$otu)

tiss <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(day <  22)


p.otu87.tiss <- tiss %>% filter(otu == 'Otu00087' & day ==21) %>% 
  ggplot(aes(x=tissue, y=value))  +
  geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  ggtitle('Relative proportion of OTU87 in each tissue a Day 21') +
  scale_fill_brewer(palette = 'Dark2')


p.otu87.fec <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(85)]) %>% 
  ggplot(aes(x=day, y=value))  +
  geom_vline(xintercept = 12, colour='black')+
  geom_boxplot(aes(x=day, y=value, group=design, fill=treatment)) + ggtitle('Relative proportion of OTU87 in feces over time') + scale_fill_brewer(palette = 'Dark2')

# ORDINATION OVER TIME #

FS2.bray2 <- vegdist(FS2.4200, method = 'bray')

# dispersion stuff, this has its own metathing so hopefully later steps still work #

dispers <- betadisper(FS2.bray2, group = meta$design)
pdispers <- permutest(dispers, pairwise = TRUE)
pdispers$pairwise
dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group
metadisp <- merge(meta, dispersdf, by = 'group')

dispgroups <- summarise(group_by(metadisp, design), average_dist=mean(dispers.distances))

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

# I think I should include this  figure in the supplement

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

FS2.nmds$group == FS2.metanmds$group
FS2.metanmds$group <- as.character(FS2.metanmds$group)
FS2.metanmds <- FS2.metanmds[match(FS2.nmds$group,FS2.metanmds$group),] # ok I think this works as long as $group isnt a factor...
FS2.nmds$group == FS2.metanmds$group

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

groupxday <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxday <- filter(groupxday, day %in% c(0:22) & tissue == 'feces')


FS2.metanmds %>% filter(tissue == 'feces' & day %in% c(0:21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_path(aes(group=treatment, color=treatment), size = 2, alpha=.7) +
  geom_point(aes(color=treatment), size=7) +
  geom_text(data=groupxday, aes(x=groupX, y=groupY, label=day))+
  scale_color_brewer(palette = 'Dark2')  + ggtitle('Fecal community structure over time', subtitle = 'Treatment group centroids') + ylab('NMDS2') + xlab('NMDS1')

# the above one is nice for showing progression over time maybe

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

# phylum level grouping #

meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)

otu <- import_mothur(mothur_shared_file = 'V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'V4.final.taxonomy')


phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]

FS2 <- phyloseq(otu, taxo)    
FS2 <- merge_phyloseq(FS2, phy_meta)  
colnames(tax_table(FS2)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.

FS2 <- subset_samples(FS2, experiment =='FS2')
FS2.phylum <- tax_glom(FS2, 'Phylum')

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

p.weaning.phyla <- filter(FS2.phylum.relabund, tissue == 'feces' & day %in% c(0,21) & variable %in% goodphy) %>%
  ggplot(aes(x=day, y=value, fill=treatment))+ geom_boxplot(aes(group=design)) +
  scale_fill_brewer(palette = 'Dark2') + expand_limits(y=0) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  ggtitle('Weaning effects on major phyla in fecal microbiota') + 
  ylab('proportion of total community')

# boom!  use it next to ordination showing D0 to D21 shift

D21phy <- levels(FS2.phylum.relabund$variable)[c(1,2,3,4,5,7,8,10,12)]

p.fecesphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'feces' & day == 21 & variable %in% D21phy) %>% 
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Feces') + 
  ylab('proportion of total community')


p.cecphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'cecum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal mucosa') +
  ylab('proportion of total community')

p.colonphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'colon' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Colonic mucosa') +
  ylab('proportion of total community')

p.ileumphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'ileum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Ileal mucosa') + 
  ylab('proportion of total community')

p.cec_cont_RNAphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'cec_cont_RNA' & day == 21 & variable %in% D21phy) %>% 
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal contents (RNA)') + 
  ylab('proportion of total community')

FS2.phylum.relabund$tissue <- factor(FS2.phylum.relabund$tissue, levels = c('ileum', 'cecum', 'cec_cont_RNA','cec_cont_DNA', 'colon', 'feces'))

# These are good, but maybe its too much... can I concentrate this at all?

p.allphyla.D21 <- filter(FS2.phylum.relabund, tissue != 'cec_cont_DNA', day == 21 & variable %in% D21phy) %>% 
  ggplot(aes(x=tissue, y=value, fill=treatment, group=design)) +
  geom_boxplot()+ scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal contents (RNA)') + 
  ylab('proportion of total community')
# hmm... perhaps....

########### but amplicon  ############

# but taxonomy

x <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
x <- unlist(strsplit(x, split = ' '))

blastn <- read.table('butreps_blastn.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)

blast <- read.table('butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)

colnames(blast) <- x
colnames(blastn) <- x

blast$tit1 <- gsub('.* \\[(.*)\\]', '\\1', blast$salltitles)
blast$tit1 <- paste(blast$tit1, blast$pident, sep = ' ')  # change this if need to alter names of but otus
blast$tit1
blast$salltitles

shared2 <- read.table('but2.shared', header = TRUE, stringsAsFactors = FALSE)
hist(rowSums(shared2[,-c(1,2,3)]), breaks = 100)

sort(rowSums(shared2[,-c(1,2,3)]))
rowSums(shared2[,-c(1,2,3)]) > 2000

blast$otu <- colnames(shared2)[-c(1,2,3)]

# this chunk creates a named vector that you can use to swap out the generic otuxxx labels for the blast results

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

meta <- read.table('V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)


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

meta$tissue[meta$tissue == 'cec_cont_DNA'] <- 'cec_cont_RNA'  # mislabeled? check this out

meta$design <- paste(meta$tissue, meta$day, meta$treatment, sep = '_')

meta$group == rownames(otu2.r)

meta <- meta[match(rownames(otu2.r), meta$group),]

meta$group == rownames(otu2.r)

# for correlation
but2 <- data.frame(otu2.r)
but2$group <- rownames(but2)

#write.table(but2, '~/FS2/correlation/but_corr/but_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Alpha Diversity #

# igraph diversity is conflicting with vegan diversity here...#

meta$invsimpson <- diversity(otu2.r, index = 'invsimpson')
meta$shannon <- diversity(otu2.r, index = 'shannon')

# Shannon

# filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
#  ggplot() + geom_boxplot(aes(day, shannon, fill = treatment)) +
#  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time')

# USE THIS FIG FOR FECES ALPHA DIV #
p.butfeces.alphatime <- filter(meta, day %in% c(0,12,15,19,21) & tissue == 'feces') %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of feces over time')

p.buttissue.alpha <- filter(meta, day == 21) %>%
  ggplot() + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Alpha diversity of tissues')

## nifty little thing to do wilcoxon tests on alpha diversity at each day between the two treatments

but_fecal_alpha_wilcox_results <- filter(meta, tissue == 'feces') %>% group_by(day) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(day, Wilcox = wilco$p.value)

# wilcoxon tests for tissues
but_tissue_alpha_wilcox_results <- filter(meta, day ==21) %>% group_by(tissue) %>%
  do(wilco = wilcox.test(invsimpson~treatment, data=., paired=FALSE)) %>%
  summarise(tissue, Wilcox = wilco$p.value)

# Ordinations #

# bray curtis calc #
meta$design <- as.factor(meta$design )

FS2.bray2 <- vegdist(otu2.r, method = 'bray')

attributes(FS2.bray2)$Labels == meta$group # still good!

# dispersion stuff #

dispers <- betadisper(FS2.bray2, group = meta$design)
pdispers <- permutest(dispers, pairwise = TRUE)
pdispers$pairwise
dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group
metadisp <- merge(meta, dispersdf, by = 'group')

dispgroups <- summarise(group_by(metadisp, design), average_dist=mean(dispers.distances))

dispgroups <- unique(inner_join(dispgroups, meta))

p.butfecal.dispersion.time <- metadisp %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) + 
  geom_boxplot() + scale_fill_brewer(palette = 'Dark2') + ylim(c(.15,.7)) + 
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")

# maybe this is useful?  polish a little bit

fecaldisptime <-  dispgroups %>% filter(tissue == 'feces' & day < 22)
fecaldisptime$day <- as.numeric(fecaldisptime$day)
p.butfecal.dispersionAVG <- ggplot(fecaldisptime, aes(x=day, y=average_dist, color = treatment)) + 
  geom_vline(xintercept = 12)+ geom_point(size = 2.5) +  geom_line() + ylim(c(.375,.543))+
  ggtitle('Community Variability (Dispersion)',
          subtitle = "Vegan's betadisper(): how much variability is there in a group's community structure?") + 
  ylab("Average distance to the group median") + scale_color_brewer(palette = "Dark2")

fecaldisptime$average_dist

p.but.tissue.dispersion <- metadisp %>% filter(day == 21) %>% 
  ggplot(aes(x=tissue, y=dispers.distances, group=design, fill=treatment)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Dark2') + ggtitle('Distance to group centroid')

# but NMDS Ordination #

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

FS2.nmds$group == FS2.metanmds$group
#FS2.metanmds$group <- as.character(FS2.metanmds$group)
FS2.metanmds <- FS2.metanmds[match(FS2.nmds$group,FS2.metanmds$group),] 
FS2.nmds$group == FS2.metanmds$group


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

# D12

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

#

p <- ggplot(data=subset(FS2.metanmds, tissue == 'cec_cont_RNA' & day == 21), aes(x=MDS1, y = MDS2)) +
  geom_point(aes(color=treatment, shape=day), size=2.1) +
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, linetype=day, color=treatment)) + 
  geom_path(data = subset(df_ell, day ==21  & tissue == 'cec_cont_RNA'), aes(x=NMDS1, y=NMDS2, group=treatmentXday, color=treatment, linetype=day), size=1) +
  ggtitle('Cec_cont_RNA Day 21') + 
  geom_point(aes(x=groupX, y=groupY, color=treatment, shape=day),size=5) + scale_color_brewer(palette = 'Dark2')

p + geom_text((aes(label=pig_num)))

# dont know if I need these daily ordinations above...#

groupxday <- unique(FS2.metanmds[,c(4,5,12:14)]) # this little dataframe is needed for annotations that dont look bad.
groupxday <- filter(groupxday, day %in% c(0:22) & tissue == 'feces')

FS2.metanmds %>% filter(tissue == 'feces' & day %in% c(0:21)) %>% ggplot(aes(x=groupX, y=groupY)) +
  geom_path(aes(group=treatment, color=treatment), size = 2, alpha=.7) +
  geom_point(aes(color=treatment), size=7) +
  geom_text(data=groupxday, aes(x=groupX, y=groupY, label=day))+
  scale_color_brewer(palette = 'Dark2')  + ggtitle('Fecal community structure over time', subtitle = 'Treatment group centroids') + ylab('NMDS2') + xlab('NMDS1')

# the above one is nice for showing progression over time maybe

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

# PERMANOVA with Adonis #

PW.Adonis <- pairwise.adonis(otu2.r,meta$design, sim.method="bray", p.adjust.m = "bonferroni")

# write.table(PW.Adonis,"Adonis-Results.csv",sep=",")

fecVSfec <- grep('feces_.* vs feces_.*' ,PW.Adonis$pairs)
colVScol <- grep('colon_.* vs colon_.*' ,PW.Adonis$pairs)
cecVScec <- grep('cecum_.* vs cecum_.*' ,PW.Adonis$pairs)
rVSr <- grep('RNA_.* vs RNA_.*' ,PW.Adonis$pairs)
ilVSil <- grep('ileum_.* vs ileum_.*' ,PW.Adonis$pairs)

good <- PW.Adonis[c(fecVSfec, colVScol, cecVScec, rVSr, ilVSil),]

fecsametime <- grep("feces_([0-9]+)_control vs feces_\\1_RPS", good$pairs)
ilsametime <- grep("ileum_([0-9]+)_control vs ileum_\\1_RPS", good$pairs)
colsametime <- grep("colon_([0-9]+)_control vs colon_\\1_RPS", good$pairs)
cecsametime <- grep("cecum_([0-9]+)_control vs cecum_\\1_RPS", good$pairs)

good2 <- good[c(fecsametime, ilsametime, colsametime, cecsametime),]
write.table(good2,"but.Adonis-Results.csv",sep=",") 

pwadon <- read.table("but.Adonis-Results.csv",sep=",")
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

# Writing stuff for correlations #

meta$sample <- paste(meta$tissue, meta$day, meta$pig_num, sep = '_')

but.feces <- otu2.r[meta$tissue == 'feces' & meta$day %in% c(0,21),]
but.cec_cont_RNA <- otu2.r[meta$tissue == 'cec_cont_RNA',]
but.cecum <-  otu2.r[meta$tissue == 'cecum',]
but.colon <-  otu2.r[meta$tissue == 'colon',]
but.ileum <-  otu2.r[meta$tissue == 'ileum',]
rownames(but.cecum)
otu <- as.data.frame(otu2.r)
otu$Group <- rownames(otu)

#write.table(meta, '~/FS2/correlation/but_corr/but_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(otu, '~/FS2/correlation/but_corr/but.all.shared', col.names = TRUE, row.names = TRUE, quote = FALSE)

# Deseq2 stuff #

otu <- import_mothur(mothur_shared_file = 'but2.shared')
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
phy_meta <- phy_meta[,-1]

FS2 <- phyloseq(otu)
FS2 <- merge_phyloseq(FS2, phy_meta)                       # combines the metadata with this phyloseq object
FS2 <- subset_samples(FS2, experiment == 'FS2')
FS2 <- prune_samples(sample_sums(FS2) > 700, FS2)  # This removes samples that have fewer than 700 sequences associated with them.
FS2 <- prune_taxa(taxa_sums(FS2) > 10, FS2)        # removes OTUs that occur less than 10 times globally

# D0 #

FS2.D0 <- subset_samples(FS2, day == 0)

sample_sums(FS2.D0)
FS2.D0 <- prune_taxa(taxa_sums(FS2.D0) > 1, FS2.D0)

rowSums(FS2.D0@otu_table)

FS2.D0.De <- phyloseq_to_deseq2(FS2.D0, ~ design)

FS2.D0.De <- DESeq(FS2.D0.De, test = "Wald", fitType = "parametric")

res.D0 = results(FS2.D0.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D0 = res.D0[which(res.D0$padj < .05), ]
sigtab.but.D0 <- as.data.frame(sigtab.but.D0)
format(sigtab.but.D0$padj, scientific = TRUE)
sigtab.but.D0$newp <- format(round(sigtab.but.D0$padj, digits = 3), scientific = TRUE)
sigtab.but.D0$Treatment <- ifelse(sigtab.but.D0$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D0$name <- butswap[rownames(sigtab.but.D0)]
sigtab.but.D0$tissue <- 'D0_feces'
sigtab.but.D0$otu <- rownames(sigtab.but.D0)

deseq.D0 <- ggplot(sigtab.but.D0, aes(x=reorder(rownames(sigtab.but.D0), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D0), y=0, label = name), size=3)+ labs(x="otu")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
deseq.D0

#  D12 #

FS2.D12 <- subset_samples(FS2, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")

res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D12 = res.D12[which(res.D12$padj < .05), ]
sigtab.but.D12 <- as.data.frame(sigtab.but.D12)

format(sigtab.but.D12$padj, scientific = TRUE)
sigtab.but.D12$newp <- format(round(sigtab.but.D12$padj, digits = 3), scientific = TRUE)
sigtab.but.D12$Treatment <- ifelse(sigtab.but.D12$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D12$name <- butswap[rownames(sigtab.but.D12)]
sigtab.but.D12$tissue <- 'D12_feces'
sigtab.but.D12$otu <- rownames(sigtab.but.D12)

deseq.D12 <- ggplot(sigtab.but.D12, aes(x=reorder(rownames(sigtab.but.D12), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D12), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D12 feces')+ coord_flip()
deseq.D12

# D15 #

FS2.D15 <- subset_samples(FS2, day == 15)

sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)

rowSums(FS2.D15@otu_table)

FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")

res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D15 = res.D15[which(res.D15$padj < .05), ]
sigtab.but.D15 <- as.data.frame(sigtab.but.D15)

format(sigtab.but.D15$padj, scientific = TRUE)
sigtab.but.D15$newp <- format(round(sigtab.but.D15$padj, digits = 3), scientific = TRUE)
sigtab.but.D15$Treatment <- ifelse(sigtab.but.D15$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D15$name <- butswap[rownames(sigtab.but.D15)]
sigtab.but.D15$tissue <- 'D15_feces'
sigtab.but.D15$otu <- rownames(sigtab.but.D15)

deseq.D15 <- ggplot(sigtab.but.D15, aes(x=reorder(rownames(sigtab.but.D15), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D15), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D15 feces')+ coord_flip()
deseq.D15

#  D19 #

FS2.D19 <- subset_samples(FS2, day == 19)

sample_sums(FS2.D19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

rowSums(FS2.D19@otu_table)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")

res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D19 = res.D19[which(res.D19$padj < .05), ]
sigtab.but.D19 <- as.data.frame(sigtab.but.D19)

format(sigtab.but.D19$padj, scientific = TRUE)
sigtab.but.D19$newp <- format(round(sigtab.but.D19$padj, digits = 3), scientific = TRUE)
sigtab.but.D19$Treatment <- ifelse(sigtab.but.D19$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D19$name <- butswap[rownames(sigtab.but.D19)]
sigtab.but.D19$tissue <- 'D19_feces'
sigtab.but.D19$otu <- rownames(sigtab.but.D19)

deseq.D19 <- ggplot(sigtab.but.D19, aes(x=reorder(rownames(sigtab.but.D19), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D19), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D19 feces')+ coord_flip()
deseq.D19

#  D21 #

FS2.D21 <- subset_samples(FS2, tissue == 'feces' & day == 21)

sample_sums(FS2.D21)
FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)

rowSums(FS2.D21@otu_table)

FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")

res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D21 = res.D21[which(res.D21$padj < .05), ]
sigtab.but.D21 <- as.data.frame(sigtab.but.D21)

format(sigtab.but.D21$padj, scientific = TRUE)
sigtab.but.D21$newp <- format(round(sigtab.but.D21$padj, digits = 3), scientific = TRUE)
sigtab.but.D21$Treatment <- ifelse(sigtab.but.D21$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D21$name <- butswap[rownames(sigtab.but.D21)]
sigtab.but.D21$tissue <- 'D21_feces'
sigtab.but.D21$otu <- rownames(sigtab.but.D21)

deseq.D21 <- ggplot(sigtab.but.D21, aes(x=reorder(rownames(sigtab.but.D21), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21 feces')+ coord_flip()
deseq.D21

# D12 15 19 21 LUMPED #

FS2.butlump <- subset_samples(FS2, tissue == 'feces' & day != 0)

sample_sums(FS2.butlump)
FS2.butlump <- prune_taxa(taxa_sums(FS2.butlump) > 1, FS2.butlump)

rowSums(FS2.butlump@otu_table)

FS2.butlump.De <- phyloseq_to_deseq2(FS2.butlump, ~ treatment)

FS2.butlump.De <- DESeq(FS2.butlump.De, test = "Wald", fitType = "parametric")

res.butlump = results(FS2.butlump.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.butlump = res.butlump[which(res.butlump$padj < .05), ]
sigtab.but.butlump <- as.data.frame(sigtab.but.butlump)

format(sigtab.but.butlump$padj, scientific = TRUE)
sigtab.but.butlump$newp <- format(round(sigtab.but.butlump$padj, digits = 3), scientific = TRUE)
sigtab.but.butlump$Treatment <- ifelse(sigtab.but.butlump$log2FoldChange >=0, "RPS", "Control")

sigtab.but.butlump$name <- butswap[rownames(sigtab.but.butlump)]
sigtab.but.butlump$tissue <- 'lumped_feces'
sigtab.but.butlump$otu <- rownames(sigtab.but.butlump)

deseq.butlump <- ggplot(sigtab.but.butlump, aes(x=reorder(rownames(sigtab.but.butlump), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.butlump), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: butlump feces')+ coord_flip()
deseq.butlump

# D21.ileum #

FS2.D21.ileum <- subset_samples(FS2, day == 21 & tissue =='ileum')

sample_sums(FS2.D21.ileum)
FS2.D21.ileum <- prune_taxa(taxa_sums(FS2.D21.ileum) > 1, FS2.D21.ileum)

rowSums(FS2.D21.ileum@otu_table)

FS2.D21.ileum.De <- phyloseq_to_deseq2(FS2.D21.ileum, ~ design)

FS2.D21.ileum.De <- DESeq(FS2.D21.ileum.De, test = "Wald", fitType = "parametric")

res.D21.ileum = results(FS2.D21.ileum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D21.ileum = res.D21.ileum[which(res.D21.ileum$padj < .05), ]
sigtab.but.D21.ileum <- as.data.frame(sigtab.but.D21.ileum)

format(sigtab.but.D21.ileum$padj, scientific = TRUE)
sigtab.but.D21.ileum$newp <- format(round(sigtab.but.D21.ileum$padj, digits = 3), scientific = TRUE)
sigtab.but.D21.ileum$Treatment <- ifelse(sigtab.but.D21.ileum$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D21.ileum$name <- butswap[rownames(sigtab.but.D21.ileum)]
sigtab.but.D21.ileum$tissue <- 'ileum'
sigtab.but.D21.ileum$otu <- rownames(sigtab.but.D21.ileum)

deseq.D21.ileum <- ggplot(sigtab.but.D21.ileum, aes(x=reorder(rownames(sigtab.but.D21.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21.ileum), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.ileum feces')+ coord_flip()
deseq.D21.ileum

#  D21.cecum #

FS2.D21.cecum <- subset_samples(FS2, day == 21 & tissue =='cecum')

sample_sums(FS2.D21.cecum)
FS2.D21.cecum <- prune_taxa(taxa_sums(FS2.D21.cecum) > 1, FS2.D21.cecum)

rowSums(FS2.D21.cecum@otu_table)

FS2.D21.cecum.De <- phyloseq_to_deseq2(FS2.D21.cecum, ~ design)

FS2.D21.cecum.De <- DESeq(FS2.D21.cecum.De, test = "Wald", fitType = "parametric")

res.D21.cecum = results(FS2.D21.cecum.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D21.cecum = res.D21.cecum[which(res.D21.cecum$padj < .1), ]
sigtab.but.D21.cecum <- as.data.frame(sigtab.but.D21.cecum)


format(sigtab.but.D21.cecum$padj, scientific = TRUE)
sigtab.but.D21.cecum$newp <- format(round(sigtab.but.D21.cecum$padj, digits = 3), scientific = TRUE)
sigtab.but.D21.cecum$Treatment <- ifelse(sigtab.but.D21.cecum$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D21.cecum$name <- butswap[rownames(sigtab.but.D21.cecum)]
sigtab.but.D21.cecum$tissue <- 'cecum'
sigtab.but.D21.cecum$otu <- rownames(sigtab.but.D21.cecum)


deseq.D21.cecum <- ggplot(sigtab.but.D21.cecum, aes(x=reorder(rownames(sigtab.but.D21.cecum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21.cecum), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.cecum feces')+ coord_flip()
deseq.D21.cecum

# D21.colon #

FS2.D21.colon <- subset_samples(FS2, day == 21 & tissue =='colon')

sample_sums(FS2.D21.colon)
FS2.D21.colon <- prune_taxa(taxa_sums(FS2.D21.colon) > 1, FS2.D21.colon)

rowSums(FS2.D21.colon@otu_table)

FS2.D21.colon.De <- phyloseq_to_deseq2(FS2.D21.colon, ~ design)

FS2.D21.colon.De <- DESeq(FS2.D21.colon.De, test = "Wald", fitType = "parametric")

res.D21.colon = results(FS2.D21.colon.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D21.colon = res.D21.colon[which(res.D21.colon$padj < .1), ]
sigtab.but.D21.colon <- as.data.frame(sigtab.but.D21.colon)

format(sigtab.but.D21.colon$padj, scientific = TRUE)
sigtab.but.D21.colon$newp <- format(round(sigtab.but.D21.colon$padj, digits = 3), scientific = TRUE)
sigtab.but.D21.colon$Treatment <- ifelse(sigtab.but.D21.colon$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D21.colon$name <- butswap[rownames(sigtab.but.D21.colon)]
sigtab.but.D21.colon$tissue <- 'colon'
sigtab.but.D21.colon$otu <- rownames(sigtab.but.D21.colon)

deseq.D21.colon <- ggplot(sigtab.but.D21.colon, aes(x=reorder(rownames(sigtab.but.D21.colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21.colon), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.colon feces')+ coord_flip()
deseq.D21.colon

#  cec cont rna  #

FS2.D21.cec_cont_RNA <- subset_samples(FS2, day == 21 & tissue =='cec_cont_RNA')

sample_sums(FS2.D21.cec_cont_RNA)
FS2.D21.cec_cont_RNA <- prune_taxa(taxa_sums(FS2.D21.cec_cont_RNA) > 1, FS2.D21.cec_cont_RNA)

rowSums(FS2.D21.cec_cont_RNA@otu_table)

FS2.D21.cec_cont_RNA.De <- phyloseq_to_deseq2(FS2.D21.cec_cont_RNA, ~ design)

FS2.D21.cec_cont_RNA.De <- DESeq(FS2.D21.cec_cont_RNA.De, test = "Wald", fitType = "parametric")

res.D21.cec_cont_RNA = results(FS2.D21.cec_cont_RNA.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.but.D21.cec_cont_RNA = res.D21.cec_cont_RNA[which(res.D21.cec_cont_RNA$padj < .05), ]
sigtab.but.D21.cec_cont_RNA <- as.data.frame(sigtab.but.D21.cec_cont_RNA)

format(sigtab.but.D21.cec_cont_RNA$padj, scientific = TRUE)
sigtab.but.D21.cec_cont_RNA$newp <- format(round(sigtab.but.D21.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.but.D21.cec_cont_RNA$Treatment <- ifelse(sigtab.but.D21.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

sigtab.but.D21.cec_cont_RNA$name <- butswap[rownames(sigtab.but.D21.cec_cont_RNA)]
sigtab.but.D21.cec_cont_RNA$tissue <- 'cec_cont_RNA'
sigtab.but.D21.cec_cont_RNA$otu <- rownames(sigtab.but.D21.cec_cont_RNA)

deseq.D21.cec_cont_RNA <- ggplot(sigtab.but.D21.cec_cont_RNA, aes(x=reorder(rownames(sigtab.but.D21.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21.cec_cont_RNA), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette="Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.cec_cont_RNA feces')+ coord_flip()
deseq.D21.cec_cont_RNA

dif_ab_but <- rbind(sigtab.but.D12,
              sigtab.but.D15,
              sigtab.but.D19,
              sigtab.but.D21,
              sigtab.but.butlump,
              sigtab.but.D21.ileum,
              sigtab.but.D21.cecum,
              sigtab.but.D21.colon,
              sigtab.but.D21.cec_cont_RNA)

write.table(dif_ab_but, 'dif_ab_but.txt', row.names = TRUE, col.names = TRUE, quote = FALSE, sep='\t')

######################### Flow cytometry analysis  ###################

meta <- read.table('V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

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

#

meta <- read.table('V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

flow_counts <- read.table('FS2_all_flow.txt', header = TRUE, check.names = FALSE, sep = '\t', as.is = TRUE, comment.char = '')

rowSums(flow_counts[,-c(1,34,35,36)]) # this tells us the total number of live single cells sorted for each sample. 

rowSums(flow_counts[,-c(1,34,35,36)]) 

flow_counts[,-c(1,34,35,36)] <- flow_counts[,-c(1,34,35,36)] / rowSums(flow_counts[,-c(1,34,35,36)]) # converts to relative abundance

# alpha div  #
meta <- flow_counts[,c(1,34,35,36)]

flow_counts %>% group_by(tissue)

H <- diversity(flow_counts[,-c(1,34,35,36)])
J <- H/log(specnumber(flow_counts[,-c(1,34,35,36)]))
meta$evenness <- J
meta$Shannon <- H

ggplot(meta) + geom_boxplot(aes(tissue, evenness, fill = treatment))
ggplot(meta) + geom_boxplot(aes(tissue, Shannon, fill = treatment))

#

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
#library(ggmosaic)
#ggplot(data=subset(x = barflow, tissue == 'cec'))+geom_mosaic(aes(x=product(pig_num), fill=attribute, weight=value), offset = 0.005) + facet_grid(~treatment) +
#  labs(x="treatment", y= "Proportion of CD3+ cells")+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#ggplotly()
#ggplot(data=subset(barflow, tissue == 'cec'), aes(x=pig_num, y=value, fill=attribute)) + geom_bar(stat = 'identity')+ labs(y = 'Proportion of live single cells')

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

ord.fit <- envfit(flow.mds, flow_counts, permutations = 9999)
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

# luimping treatments #

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

#adonis.flow <- adonis(flow_counts~meta$treatment, permutations = 9999)
#adon.pval <- adonis.flow$aov.tab$`Pr(>F)`[1]
#adonis.flow$aov.tab

ord.fit <- envfit(flow.mds, flow_counts, permutations = 9999)
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
#  Just cecal  #

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

ord.fit <- envfit(flow.mds, flow_counts, permutations = 9999)
spp.scrs <- as.data.frame(scores(ord.fit, display = 'vectors'))
spp.scrs <- spp.scrs[which(ord.fit$vectors$pvals < 0.01),]
spp.scrs <- spp.scrs/7
spp.scrs <- cbind(spp.scrs, cell_type = rownames(spp.scrs))
colnames(spp.scrs) <- c('MDS1', 'MDS2', 'cell_type')


p.floword <- ggplot(flow.meta.nmds, aes(x=MDS1, y = MDS2)) +
  #geom_segment(data = spp.scrs, aes(x=0, xend=spp.scrs$MDS1, y=0, yend=spp.scrs$MDS2), alpha=.5)+
  geom_point(data=flow.meta.nmds, aes(color=treatment), size=2.5) +
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  ggtitle('Cecal T-cell community similarity: NMDS ordination using Bray-Curtis distances', subtitle = paste('PERMANOVA p-value = ', adon.pval)) + 
  geom_segment(data = flow.meta.nmds,
               aes(x=flow.meta.nmds$MDS1,
                   xend=flow.meta.nmds$centroidX,
                   y=flow.meta.nmds$MDS2,
                   yend=flow.meta.nmds$centroidY,
                   color=treatment)) + 
  #geom_text_repel(data=spp.scrs, aes(MDS1, MDS2, label = cell_type), size=4, alpha=.7)+
  scale_color_brewer(palette="Dark2")

p.floword

#  boxplots  #

fin <- data.frame()
# this for loop does wilcoxon tests for each cell type and saves the results into a new column that I can use to label??
barflow$attribute <- factor(barflow$attribute)

barflow$tissue <- factor(barflow$tissue)


barflowtoo <- barflow

barflowtoo$tisXatt <- paste(barflowtoo$tissue, barflowtoo$attribute, sep = '_')


barflowtoo$tisXatt <- factor(barflowtoo$tisXatt)

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

# just the significantly different populations in all tissues
p.flow.allsig <- filter(fin, pvalue < 0.05) %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+
  facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette="Dark2")+
  labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))


p.flow.allsig


# just the significantly different populations only in the cecum
p.flow.cecsig <- filter(fin, pvalue < 0.1 & tissue == 'cec') %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+
  facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette="Dark2")+
  labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))


p.flow.cecsig

################# Nitrophenyl linked sugars assay ####################

# D0 stuff #
plateinfoD0 <- read.table('D0_plate_layout.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plateinfoD0 <- plateinfoD0[,-2]
plateinfoD067 <- read.table('D0_6_7_plate_layout.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plateinfoD067 <- plateinfoD067[,-2]

D0_plate2_3 <- read.table('D0_2_3.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
D0_plate4_5 <- read.table('D0_4_5.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
D0_plate6_7 <- read.table('D0_6_7.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
D0platesulf <- read.table('sulf_0_21.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)


# have to do 102 on the filter because well is a factor and it is filtering based on factor level, level 102 corresponds to 101 because blank is a level
D0plate2_3.melt <- melt(data = D0_plate2_3, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
D0plate2_3.melt$substrate <- ifelse(as.numeric(D0plate2_3.melt$well) <102, 'N-acetyl-B-D-glucosamine', 'B-D-glucopyranoside')
D0plate4_5.melt <- melt(data = D0_plate4_5, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
D0plate4_5.melt$substrate <- ifelse(as.numeric(D0plate4_5.melt$well) <102, 'B-D-xylopyranoside', 'a-D-galactopyranoside')
D0plate6_7.melt <- melt(data = D0_plate6_7, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
D0plate6_7.melt$substrate <- ifelse(as.numeric(D0plate6_7.melt$well) <102, 'B-D-galactopyranoside', 'a-L-fucopyranoside')
D0platesulf.melt <- melt(data = D0platesulf, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
D0platesulf.melt$substrate <- '4-nitrophenyl sulfate'
D0platesulf.melt <- filter(D0platesulf.melt, as.numeric(well) < 102)

D0plate6_7.melt <- merge(D0plate6_7.melt, plateinfoD067, by = 'well', all = TRUE )
D0 <- rbind(D0plate2_3.melt, D0plate4_5.melt, D0platesulf.melt)
D0 <- merge(D0, plateinfoD0, by = 'well', all = TRUE)
D0 <- rbind(D0, D0plate6_7.melt)


# D21 stuff #
plateinfo <- read.table('D21_plate_layout.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plateinfo <- plateinfo[,-2]
plate2_4 <- read.table('D21_2_4.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plate3_5 <- read.table('d21_3_5.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plate6_7 <- read.table('d21_6_7.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
platesulf <- read.table('sulf_0_21.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)


# have to do 102 on the filter because well is a factor and it is filtering based on factor level, level 102 corresponds to 101 because blank is a level
plate2_4.melt <- melt(data = plate2_4, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
plate2_4.melt$substrate <- ifelse(as.numeric(plate2_4.melt$well) <102, 'N-acetyl-B-D-glucosamine', 'B-D-xylopyranoside')
plate3_5.melt <- melt(data = plate3_5, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
plate3_5.melt$substrate <- ifelse(as.numeric(plate3_5.melt$well) <102, 'B-D-glucopyranoside', 'a-D-galactopyranoside')
plate6_7.melt <- melt(data = plate6_7, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
plate6_7.melt$substrate <- ifelse(as.numeric(plate6_7.melt$well) <102, 'B-D-galactopyranoside', 'a-L-fucopyranoside')
platesulf.melt <- melt(data = platesulf, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
platesulf.melt$substrate <- '4-nitrophenyl sulfate'
platesulf.melt <- filter(platesulf.melt, as.numeric(well) > 101)

master <- rbind(plate2_4.melt, plate3_5.melt, plate6_7.melt, platesulf.melt)

master <- merge(plateinfo, master, by = 'well', all = TRUE)

master <- rbind(D0, master)

master$seconds <- period_to_seconds(hms(master$Time))

control <- c(67,68,69,70,71,72,73,81,82,83,84,85,86,87)

master$treatment <- ifelse(master$pig %in% control, 'control', 'RPS')

master <- master[order(master$seconds),]

master$treatment[master$pig == 'NTC'] <- 'NC'
master$day <- factor(master$day)
master <- master[master$pig != '',]

master$group <- paste(master$substrate, master$day, master$tissue, master$pig)

filter(master, seconds >4000 & seconds<10000 & pig != 'NTC') %>%
  ggplot(aes(x=seconds, y=abs405, group=group, color=treatment, linetype=day)) +
  geom_path()+ facet_wrap(~substrate, scales = 'free')

# calc rates with lm and check on r2 to see how linear my kinetic data really is.  

master.results <- filter(master, seconds >4000 & seconds<10000) %>%
  group_by(group) %>%
  do(data.frame(setNames(as.list(coef(lm(abs405~seconds, data=.))[2]),
                         c('rate'))))

master.results2 <- filter(master, seconds >4000 & seconds<10000) %>%
  group_by(group) %>%
  do(data.frame(setNames(as.list(summary(lm(abs405~seconds, data=.))$r.squared), c('r.squared'))))


master.results <- merge(master.results, master.results2, by = 'group', all = TRUE)

master.results$substrate <- gsub('(.*e) ([0-9]+) (.*[A-Za-z]+) [0-9A-Z]+', '\\1', master.results$group)
master.results$day <- gsub('(.*e) ([0-9]+) (.*[A-Za-z]+) [0-9A-Z]+', '\\2', master.results$group)
master.results$tissue <- gsub('(.*e) ([0-9]+) (.*[A-Za-z]+) [0-9A-Z]+', '\\3', master.results$group)
master.results$pig <- gsub('(.*e) ([0-9]+) (.*[A-Za-z]+) ([0-9A-Z]+)', '\\4', master.results$group)
master.results$treatment <- ifelse(master.results$pig %in% control, 'control', 'RPS')
master.results$treatment[master.results$pig == 'NTC'] <- 'NTC'
master.results$design <- paste(master.results$substrate, master.results$day, master.results$tissue, master.results$treatment, sep = ' ')
master.results$group <- paste(master.results$day, master.results$tissue, master.results$treatment)
master.results$group2 <- paste(master.results$day, master.results$tissue, master.results$pig)




comm_mat <- select(master.results, c(10,2,4)) %>% spread(key = substrate, value = rate, drop = TRUE)
comm_mat <- comm_mat[-grep('NTC', comm_mat$group2),]
comm_mat$day <- gsub('([0-9]+) ([A-Za-z]+_?[A-Za-z]+) ([0-9]+)', '\\1', comm_mat$group2)
comm_mat$tissue <- gsub('([0-9]+) ([A-Za-z]+_?[A-Za-z]+) ([0-9]+)', '\\2', comm_mat$group2)
comm_mat$pig <- gsub('([0-9]+) ([A-Za-z]+_?[A-Za-z]+) ([0-9]+)', '\\3', comm_mat$group2)
comm_mat$treatment <- ifelse(comm_mat$pig %in% control, 'control', 'RPS')
comm_mat$`4-nitrophenyl sulfate`[comm_mat$`4-nitrophenyl sulfate` < 0] <- 0 # one of the pigs had a negative rate....
comm_mat$muc_deg <- rowSums(comm_mat[,c(2,4,5,8)])

#write.table(comm_mat, file='glycoside_hydrolase_matrix.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)


filter(comm_mat, tissue == 'feces') %>%
  ggplot(aes(x=day, y=muc_deg, fill = treatment))+
  geom_boxplot()

master.results$substrate <- factor(master.results$substrate,
                                   levels = c("4-nitrophenyl sulfate",
                                              "a-L-fucopyranoside",
                                              "N-acetyl-B-D-glucosamine",
                                              "a-D-galactopyranoside",
                                              "B-D-galactopyranoside",
                                              "B-D-xylopyranoside",
                                              "B-D-glucopyranoside"))

#master.results$day <- as.numeric(master.results$day)
filter(master.results, pig != 'NTC' & tissue == 'feces' & r.squared > 0.7 & substrate != 'B-D-galactopyranoside') %>%
  ggplot(aes(x=day, y=rate, group=design)) +
  geom_boxplot(aes(fill=treatment)) + 
  scale_fill_brewer(palette = 'Dark2') + facet_wrap(~substrate, scales = 'free') + expand_limits(y=0)


filter(master.results, pig != 'NTC' & tissue != 'feces' & r.squared > 0.7 & substrate != 'B-D-galactopyranoside') %>%
  ggplot(aes(x=day, y=rate, group=design)) +
  geom_boxplot(aes(fill=treatment)) + 
  scale_fill_brewer(palette = 'Dark2') + facet_wrap(~substrate, scales = 'free') + expand_limits(y=0)

master.results.filter <- filter(master.results, pig != 'NTC' & tissue == 'feces' & r.squared > 0.7)
master.results.filter2 <- filter(master.results, pig != 'NTC' & tissue == 'cec_cont' & r.squared > 0.7)


tests <- data.frame()
for (x in levels(factor(master.results.filter$substrate))){
  print(x)
  temp <- filter(master.results.filter, substrate == x)
  temp <- pairwise.wilcox.test(temp$rate, temp$design)$p.value
  temp <- melt(temp)
  tests <- rbind(tests, temp) 
}

tests <- tests[!is.na(tests$value),]
tests <- tests[tests$value < 0.05,]

tests.cec <- data.frame()
for (x in levels(factor(master.results.filter2$substrate))){
  print(x)
  temp <- filter(master.results.filter2, substrate == x)
  temp <- pairwise.wilcox.test(temp$rate, temp$design)$p.value
  temp <- melt(temp)
  tests.cec <- rbind(tests.cec, temp) 
}

tests.cec <- tests.cec[!is.na(tests.cec$value),]
tests.cec <- tests.cec[tests.cec$value < 0.1,]

# hmm <- master.results %>% filter(substrate == 'N-acetyl-B-D-glucosamine' & tissue == 'cec_cont')
# master.results %>% filter(substrate == 'N-acetyl-B-D-glucosamine'& tissue != 'cec_cont' & treatment != 'NTC') %>% ggplot(aes(x=day, y=rate, group=design)) +
#   geom_boxplot(aes(fill=treatment)) + geom_jitter()+
#   facet_wrap(~substrate, scales = 'free')

#####################  VFAs  ####################

vfa.cec <- read.table('cec_redo_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfa.fec <- read.table('fecalplasma_R.txt', header = TRUE, stringsAsFactors = FALSE)

vfas <- rbind(vfa.cec, vfa.fec)
vfas[,2:16] <- vfas[,2:16]*3
vfas$BCFA <- vfas$isobutyrate + vfas$isovalerate
vfas$total <- rowSums(vfas[,2:16])
colnames(vfas)[8] <- 'lactate'


vfas$sample <- paste(vfas$tissue, 21, vfas$pig_num, sep = '_')
vfas$treatment <- ifelse(vfas$pig_num %in% c(67,68,69,70,71,72,73,81,82,83,84,85,86,87), 'control', 'RPS')
vfas$design <- paste(vfas$tissue, vfas$treatment, sep = '_')
colnames(vfas)
vfas.melt <- melt(data = vfas, id.vars = c(1,17, 20:22)) 


#vfas.melt  %>% ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
#vfas.melt %>% filter(tissue == 'portal') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
#vfas.melt %>% filter(tissue == 'cecum') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')

p.vfas.final <- vfas.melt %>% filter(variable %in% c('propionate', 'butyrate', 'lactate','BCFA', 'succinate', 'total') & tissue != 'portal') %>%
  ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  geom_boxplot() + facet_wrap(~variable, scales = 'free', nrow = 3, ncol = 2) + scale_fill_brewer(palette = 'Dark2') + ylab('Concentration (mM)') + ggtitle('VFA concentrations at day 21')

p.vfas.final

########### NEED TO ADD VFA STATISTICAL TESTS HERE ###############

############### BUTNET ###############

# but taxonomy #

xx <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
xx <- unlist(strsplit(xx, split = ' '))

blast <- read.table('butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
shared2 <- read.table('but2.shared', header = TRUE, stringsAsFactors = FALSE)

colnames(blast) <- xx

blast$tit1 <- gsub('.* \\[(.*)\\]', '\\1', blast$salltitles)
blast$tit1

blast$otu <- colnames(shared2)[-c(1,2,3)]
blast$gen <- gsub('(\\[?[A-Za-z]+\\]? [A-Za-z]+).*', '\\1', blast$tit1)
blast$gen <- gsub('uncultured ', '', blast$gen)
blast$gen <- gsub(' sp', '', blast$gen)
blast$gen <- gsub(' bacterium', '', blast$gen)
blast$gen <- gsub('\\[', '', blast$gen)
blast$gen <- gsub('\\]', '', blast$gen)

blast$lab <- paste(blast$gen, blast$otu, sep = ' ')

#### this chunk creates a named vector that you can use to swap out the generic otuxxx labels for the blast results

butswip <- blast$otu
names(butswip) <- blast$lab
butswap <- names(butswip)
names(butswap) <- butswip

butswap[blast$otu]


# Read in data #

tax <- extract_mothur_tax('V4.final.taxonomy')
tax <- tax[,-c(2,3)]
swap <- otu_tax_labels(tax)

meta <- read.table('16S_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE) 

otu <- read.table('16S_shared_forcorr.txt', header = TRUE) %>% filter(group %in% meta$group)

rownames(otu) <- meta$sample

# differential abundance data for but and 16S otus, previously calculated by DeSeq2
dif_ab_16S <- read.table('dif_ab_16s.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
dif_ab_but <- read.table('dif_ab_but.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)


dif_ab_16S$node <- swap[dif_ab_16S$otu]
dif_ab_but$node <- butswap[dif_ab_but$otu]

enrich <- rbind(dif_ab_16S[,c(13,14,15,16,17)], dif_ab_but[,c(7,8,10,11,12)])

# but data #

but <- read.table('but_shared_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
butmeta <- read.table('but_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
but$group == butmeta$group
butmeta$sample <- paste(butmeta$tissue, butmeta$day, butmeta$pig_num, sep = '_')
rownames(but) <- butmeta$sample

# getting everything in order #
# Only feces used in this correlation

meta <- meta %>% filter(tissue == 'feces')

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
rowSums(but)

# sqrt transformation, comment out to revert to relative abundance.  I think the sqrt transform
# is appropriate here, It is common throughout ecology and it serves to correct the well established
# amplification bias seen in amplicon based sequencing projects
# Some sequences are amplified preferentially so their abundance is inflated, other sequences are amplified
# poorly so their abundance is underestimated.  A sqrt transformation helps correct this by reducing the  relative abundance
# of highly abundant features and increasing the relative abundance of less abundant features, yet the ranks of abundances are unaffected.
# It probably doesnt make a difference.

but <- sqrt(but)
otu <- sqrt(otu)

otu <- otu/rowSums(otu)
but <- but/rowSums(but)

rowSums(but)
rowSums(otu)

otu <- otu[,((colSums(otu)/length(otu[,1]))*100) >0.001] # removes 16S otus with less than 0.01% average abundance
but <- but[,((colSums(but)/length(otu[,1]))*100) >0.001] # removes but otus with less than 0.01% average abundance

colnames(otu) <- swap[colnames(otu)]
colnames(but) <- butswap[colnames(but)]

rownames(but) == rownames(otu)

# fecals only #
but_16S.all <- ccrepe(but, otu, verbose = FALSE, min.subj = 7)
but.all <- ccrepe(but, verbose = FALSE, min.subj = 7)
otu.all <-  ccrepe(otu, verbose = FALSE, min.subj = 7)


but_16S.all.sigs <- ccrepe_to_ggnet(but_16S.all)
but.all.sigs <- ccrepe_to_ggnet(but.all)
otu.all.sigs <- ccrepe_to_ggnet(otu.all, pcut = 0.01)

nodes <- rbind(gather_nodes(otu, '16S'), 
               gather_nodes(but, 'but'))


# this limits the butnet to only features that are diffabund

but_16S.all.sigs <- but_16S.all.sigs[but_16S.all.sigs$from %in% enrich$node & but_16S.all.sigs$to %in% enrich$node,] # one of the nodes in a but-otu connection must be diff abund between the two groups

but.all.sigs <- but.all.sigs[but.all.sigs$from %in% enrich$node & but.all.sigs$to %in% enrich$node,] # but-but connections must both be diff abund
otu.all.sigs <- otu.all.sigs[otu.all.sigs$from %in% enrich$node & otu.all.sigs$to %in% enrich$node,] # otu-otu connections must both be diff abund

all <- rbind(but.all.sigs, but_16S.all.sigs, otu.all.sigs)

all <- fortify(as.edgedf(all), nodes)

filtered4 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 4)

p.butnet1 <- ggplot(filtered4, aes(from_id = from_id, to_id = to_id, label=from, color=type)) + 
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
           aes(color = type, label = from_id),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
           labelgeom = 'text') +
  theme_net()

#p.butnet1


enrich$from2 <- swap[enrich$otu]
net <- as.data.frame(ggplot_build(p.butnet1)$data)

colnames(dif_ab_16S[,c(13,14,15,16,17)])
colnames(dif_ab_but[,c(7,8,10,11,12)])

colnames(enrich)[5] <- 'from'
enrich <- enrich[!(enrich$tissue %in% c('ileum', 'cecum', 'cec_cont_RNA', 'colon')),]

net2 <- merge(enrich, net, by = 'from', all = TRUE)

net2 <- net2[!is.na(net2$x),]
net2$type <- NA
net2$type[grep('Otu[0-9][0-9][0-9]', net2$from)] <- 'but'
net2$type[grep('Otu[0-9][0-9][0-9][0-9][0-9]', net2$from)] <- '16S'

net2$label <- butswap[net2$label]
net2$label[is.na(net2$label)] <- swap[net2$from[is.na(net2$label)] ]

colnames(net2)
net3 <- unique(net2[,c(1,3,9,10)])
net4 <- unique(net2[,c(9,10,11,12)])

net4 <- unique(na.omit(net4[net4 != net4[,c(3,4,1,2)],]))

net3$type <- NA
net3$type[grep('Otu.....', net3$from)] <- '16S'
net3$type[grep('Otu...$', net3$from)] <- 'but'

net3$from <- gsub('(.*)\\(.*\\).*','\\1',net3$from)
net3$from <- gsub('(.*) Otu...', '\\1', net3$from)

p.butnet.final <- ggplot(net2, aes(x=x, y=y)) + 
  
  geom_point(data=net3, aes(color = Treatment), alpha=.5, size=5, show.legend = TRUE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=8, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=10, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=12, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=14, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=16, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=18, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=20, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=22, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=24, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=26, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=28, show.legend = FALSE) +
  geom_point(data=net3, aes(color = Treatment), alpha=.1, size=30, show.legend = FALSE) +
  
  geom_segment(data=net4, aes(x=x, y=y, xend=xend, yend=yend), alpha=0.5) +
  
  geom_point(aes(shape = type, fill = type), size=5) + 
  geom_label_repel(data=net3, aes(label=from, fill=type), size=3, alpha=.7, show.legend = FALSE) +
  scale_color_brewer(palette = 'Dark2') + scale_shape_manual(values=c(21,24)) +
  scale_fill_brewer(palette = 'Accent') + theme(panel.background = element_blank(),
                                                axis.title = element_blank(),
                                                axis.text = element_blank(),
                                                axis.ticks = element_blank())
p.butnet.final


############################# Correlation Network ####################

# read in data #

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

# uncomment for only cd3+ cells
# flowc <- read.table("Flow_abs_counts.txt", header = TRUE, comment.char = '', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
# rownames(flowc) <- flowc$Pig
# flowc <- flowc[-15,-c(1,2,3)]
# flowr <- flowc/rowSums(flowc)
# rowSums(flowr)  #checking rows sum to 1
# 
# colnames(flowr)

# meta <- read.table('~/FS2/16s/V4/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

meta <- read.table('V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)
meta <- meta[meta$day == 21 & meta$tissue == 'cecum',]

# misc #
##################### NEED TO CHANGE THIS SECTION HERE #################
misc <- read.table('miscforcorr.txt', header = TRUE, stringsAsFactors = FALSE, sep = ',', as.is = TRUE, check.names = FALSE)

rownames(misc) <- misc$pig_num

misc2 <- misc
misc2$treatment <- c(rep('control', 7), rep('RPS', 7))

boxplot(misc2$`ng IgA/mg dry contents`~misc2$treatment)
wilcox.test(misc2$`ng IgA/mg dry contents`~misc2$treatment)

############# Need to make IgA figure still ###########

tax <- extract_mothur_tax('V4.final.taxonomy')
swap <- otu_tax_labels(tax)

# VFAS #
vfa.cec <- read.table('cec_redo_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfa.fec <- read.table('fecalplasma_R.txt', header = TRUE, stringsAsFactors = FALSE) 

vfas <- rbind(vfa.cec, vfa.fec)
vfas[,2:16] <- vfas[,2:16]*3  # correct for the dillution in our GC protocol
vfas$BCFA <- vfas$isobutyrate + vfas$isovalerate
vfas$total <- rowSums(vfas[,2:16])
colnames(vfas)[8] <- 'lactate'

vfas$sample <- paste(vfas$tissue, 21, vfas$pig_num, sep = '_')
vfas$treatment <- ifelse(vfas$pig_num %in% c(67,68,69,70,71,72,73,81,82,83,84,85,86,87), 'control', 'RPS')
vfas$design <- paste(vfas$tissue, vfas$treatment, sep = '_')
colnames(vfas)

cecbact <- read.table('cecum.shared', header = TRUE)
rownames(cecbact) <- meta$pig_num
cecbact <- cecbact[,-c(1,2,3)]
rowSums(cecbact)  #31441 seqs per sample
cecbact <- cecbact/rowSums(cecbact)  # converts to relative abundance
cecbact <- cecbact[,colSums(cecbact)>.0001] #removes OTUs with less than 0.01% abundance across all cecal tissue samples
colnames(cecbact) <- swap[colnames(cecbact)]

# qpcr deltaCTs

qPCR <- read.table('Cec_tissue_forcorr.txt', header = TRUE)
qPCR <- aggregate(x=qPCR, data=qPCR, by=list(qPCR$pig_num), FUN=mean)
rownames(qPCR) <- qPCR$pig_num
qPCR <- qPCR[,-c(1,2)]

qPCR$ACTb/qPCR$DEFb1
deltaCT <- qPCR - qPCR$ACTb

deltaCT.m <- as.matrix(deltaCT)
deltaCT.m <- 1/deltaCT.m[,-1]

vfas.m <- as.matrix(filter(vfas, tissue == 'cecum')[,-c(1,7,12,16,17,19,20:22)])

flowrm <- as.matrix(flowr)
cecbactm <- as.matrix(cecbact)
cecbactcontm <- as.matrix(cecbactcont)

deltaCT.m[,8]

misc2$mucexpr <- deltaCT.m[,8]

miscm <- as.matrix(misc)
miscm <- miscm[,c(2,3,4,5)]

# Correlation calculations #

# THIS IS THE ORIGINAL SWITCH BACK IF WEIRD
TxB <- ccrepe(x = cecbact, y = flowr, min.subj = 7, verbose = FALSE)
bac <- ccrepe(x = cecbact, min.subj = 7, verbose = FALSE, compare.within.x = TRUE)
im <- ccrepe(x = flowr, min.subj = 7, verbose = FALSE)

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

# Gathering node data #

nodes <- rbind(gather_nodes(flowrm, 'T-cell'), 
               gather_nodes(cecbact, '16s'), 
               gather_nodes(vfas.m, 'VFA'),
               gather_nodes(deltaCT.m, 'mRNA'),
               gather_nodes(misc))

nodes$type[grep('CD3-', nodes$node)] <- 'CD3neg'
nodes$type[grep('CD3\\+', nodes$node)] <- 'CD3pos'
# Convert to geomnet format #

TxB.sigs <- ccrepe_to_ggnet(TxB, spearcut = 0.6, pcut = 0.05)
im.sigs <- ccrepe_to_ggnet(im, spearcut = 0.6, pcut = 0.05)

vf.sig <- rcorr_to_ggnet(vfaVSflow, pcut = 0.015)
qPCR_flow.sig <- rcorr_to_ggnet(qPCRvsFlow, pcut = 0.05)
misc_qPCR.sig <- rcorr_to_ggnet(miscVSqPCR, pcut = 0.05)
misc_flow.sig <- rcorr_to_ggnet(miscVSflow, pcut = 0.05)
vfa_misc.sig <- rcorr_to_ggnet(vfaVSmisc, pcut = 0.05)
vfa_qPCR.sig <- rcorr_to_ggnet(vfaVSqPCR, pcut = 0.05)

# these below need to be pruned #

vfbac.sig <- rcorr_to_ggnet(vfaVSbact, pcut = 0.001)
qPCR_bact.sig <- rcorr_to_ggnet(qPCRvsBac, pcut = 0.001)
misc_bac.sig <- rcorr_to_ggnet(miscVSbac, pcut = 0.05)
bac.sigs <- ccrepe_to_ggnet(bac, pcut = 0.05, spearcut = 0.6)


# removing stuff from bac.sigs and vf.sigs #

# removes tcell-tcell correlations calculated by rcorr()

vf.sig <- vf.sig[!(grepl('CD25', vf.sig$from) & grepl('CD25', vf.sig$to)),] 
misc_flow.sig <- misc_flow.sig[!(grepl('CD25', misc_flow.sig$from) & grepl('CD25', misc_flow.sig$to)),] 
misc_bac.sig<- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),] 

vfbac.sig <- vfbac.sig[!(grepl('Otu', vfbac.sig$from) & grepl('Otu', vfbac.sig$to)),] 

qPCR_bact.sig <- qPCR_bact.sig[!(grepl('Otu', qPCR_bact.sig$from) & grepl('Otu', qPCR_bact.sig$to)),] 
misc_bac.sig <- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),] 

# this makes it so that bac to bac correlations are limited to those OTUs which also correlate with other features

bac.sigs <- rbind(bac.sigs[bac.sigs$from %in% TxB.sigs$to,],
                  bac.sigs[bac.sigs$to %in% TxB.sigs$to,],
                  bac.sigs[bac.sigs$from %in% vfbac.sig$to,],
                  bac.sigs[bac.sigs$to %in% vfbac.sig$from,])

qPCR_flow.sig <- qPCR_flow.sig[!(grepl('CD25', qPCR_flow.sig$from) & grepl('CD25', qPCR_flow.sig$to)),]

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

filtered <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 5)

# Plotting 

set.seed(12812)

p <- ggplot(filtered, aes(from_id = from_id, to_id = to_id, label=from_id, color=type)) + 
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
           aes(color = type),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = FALSE, fontsize=3.5, singletons = FALSE,labelcolour="black",
           labelgeom = 'label') +
  theme_net()
p

stuff <- ggplot_build(p)
newplot <- stuff$data[[1]]

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

#

proc.time() - ptm

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
