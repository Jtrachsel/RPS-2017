setwd('~/Documents/FS2/')

ptm <- proc.time()

set.seed(357)

library(vegan)
library(tidyverse)
library(ggrepel)
library(phyloseq)
library(DESeq2)
library(reshape2)
library(ccrepe)
library(geomnet)
library(Hmisc)
library(ggrepel)
library(funfuns)
system('ls')

######### 16S rRNA gene amplicons  #########
# reading in the data #

tax <- read.table('./data/V4.final.taxonomy', header = TRUE, stringsAsFactors = FALSE)
shared <- read.table('./data/V4.final.shared', header=TRUE, stringsAsFactors = FALSE)
rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

# NTC stuff #

ntc <- shared[rownames(shared) == 'NTC',]
rowSums(ntc) # 71 sequences in my NTC
ntc[,which(ntc >1)]
# this could be due to index read errors or low level contamination
# end NTC #

shared <- shared[rownames(shared) != 'NTC',] # removes NTC from otu table
rownames(shared)
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
write.table(meta_for_corr, './data/16S_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
shared2 <- shared[meta$experiment == 'FS2',]
shared2$group <- rownames(shared2)
write.table(shared2, './data/16S_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(design, file = 'V4.design', quote = FALSE, sep = '\t', row.names = FALSE)
#write.table(FS2.accnos, file = 'FS2.accnos', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
#write.table(meta, file = './output/V4.metadata.txt', quote = FALSE, sep = '\t', row.names = FALSE)
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
# SUPPLEMENT
p.fectime.alpha.16S <- filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ggtitle('Alpha diversity of feces over time') +
  annotate('segment', x=4.7, xend=5.2, y=30, yend=30, size=.3) + annotate('text', x=4.9, y=31, label='* p=0.029', size=4)

p.fectime.alpha.16S

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
# supplement
p.tissue.alpha.16S <- ggplot(Tissue21.meta) + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  ggtitle('Alpha diversity of tissues') + annotate('text', x=4, y=36.5, label= '* p = 0.022 *')+
  annotate('segment', x=3.8, xend=4.2, y=35, yend=35)

p.tissue.alpha.16S

# Deseq2 differential abundance #

otu <- import_mothur(mothur_shared_file = './data/V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './data/V4.final.taxonomy')

meta <- read.table(file = './output/V4.metadata.txt', sep = '\t', header = TRUE)
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
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
deseq.D12

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
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
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


sigtab.D21[reorder(rownames(sigtab.D21), sigtab.D21$log2FoldChange),]

sigtab.D21 <- sigtab.D21[order(sigtab.D21$log2FoldChange),]

respire.feces.d21 <- sort(c(grep('Camp', sigtab.D21$Genus),
  grep('Sutt', sigtab.D21$Genus),
  grep('Strep', sigtab.D21$Genus),
  #grep('Trep', sigtab.D21$Genus),
  grep('Helic', sigtab.D21$Genus),
  grep('Muci', sigtab.D21$Genus)))

sigtab.D21$Family <- as.character(sigtab.D21$Family)
sigtab.D21$Genus <- as.character(sigtab.D21$Genus)

resp.anno <- sigtab.D21[respire.feces.d21,]
resp.anno$x <- respire.feces.d21
resp.anno$lab <- paste(resp.anno$Family, resp.anno$Genus)

sigtab.D21$Family[respire.feces.d21] <- ''
sigtab.D21$Genus[respire.feces.d21] <- ''

p.deseq.D21.feces <- ggplot(sigtab.D21, aes(x=reorder(rownames(sigtab.D21), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(y=0, label = paste(Family, Genus, sep = ' ')),fontface='italic', size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + geom_text(data=resp.anno, aes(x=x, y=0, label=lab), fontface='bold.italic', size=4)+
  theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Feces')+ coord_flip()
p.deseq.D21.feces

# MAYBE THIS ONE FOR FIG 6

# 
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

sigtab.ileum <- sigtab.ileum[order(sigtab.ileum$log2FoldChange),]

respire.feces.d21 <- sort(c(grep('Camp', sigtab.ileum$Genus),
                            grep('Sutt', sigtab.ileum$Genus),
                            grep('Strep', sigtab.ileum$Genus),
                            #grep('Trep', sigtab.ileum$Genus),
                            grep('Helic', sigtab.ileum$Genus),
                            grep('Muci', sigtab.ileum$Genus)))

sigtab.ileum$Family <- as.character(sigtab.ileum$Family)
sigtab.ileum$Genus <- as.character(sigtab.ileum$Genus)


resp.anno <- sigtab.ileum[respire.feces.d21,]
resp.anno$x <- respire.feces.d21
resp.anno$lab <- paste(resp.anno$Family, resp.anno$Genus)

sigtab.ileum$Family[respire.feces.d21] <- ''
sigtab.ileum$Genus[respire.feces.d21] <- ''


p.deseq.ileum <- ggplot(sigtab.ileum, aes(x=reorder(rownames(sigtab.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) + geom_text(aes(x=rownames(sigtab.ileum), y=0, label = paste(Family, Genus, sep = ' ')), fontface='italic', size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + geom_text(data=resp.anno, aes(x=x, y=0, label=lab), fontface='bold.italic', size=4)+
  theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank()) +
  ggtitle('Ileal mucosa')+ coord_flip()

# THIS ONE for fig5 #
p.deseq.ileum

sigtab.cecum = res.cecum[which(res.cecum$padj < 0.05), ]
sigtab.cecum = cbind(as(sigtab.cecum, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cecum), ], "matrix"))
format(sigtab.cecum$padj, scientific = TRUE)
sigtab.cecum$newp <- format(round(sigtab.cecum$padj, digits = 3), scientific = TRUE)
sigtab.cecum$Treatment <- ifelse(sigtab.cecum$log2FoldChange >=0, "RPS", "Control")

sigtab.cecum <- sigtab.cecum[order(sigtab.cecum$log2FoldChange),]

respire.feces.d21 <- sort(c(grep('Camp', sigtab.cecum$Genus),
                            grep('Sutt', sigtab.cecum$Genus),
                            grep('Strep', sigtab.cecum$Genus),
                            #grep('Trep', sigtab.cecum$Genus),
                            grep('Helic', sigtab.cecum$Genus),
                            grep('Muci', sigtab.cecum$Genus)))

sigtab.cecum$Family <- as.character(sigtab.cecum$Family)
sigtab.cecum$Genus <- as.character(sigtab.cecum$Genus)


resp.anno <- sigtab.cecum[respire.feces.d21,]
resp.anno$x <- respire.feces.d21
resp.anno$lab <- paste(resp.anno$Family, resp.anno$Genus)

sigtab.cecum$Family[respire.feces.d21] <- ''
sigtab.cecum$Genus[respire.feces.d21] <- ''



p.deseq.cecum <- ggplot(sigtab.cecum, aes(x=reorder(rownames(sigtab.cecum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) + geom_text(aes(x=rownames(sigtab.cecum), y=0, label = paste(Family, Genus, sep = ' ')),fontface='italic', size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + geom_text(data=resp.anno, aes(x=x, y=0, label=lab), fontface='bold.italic', size=4)+
  theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Cecal mucosa')+ coord_flip()
############## THIS ONE for fig 5 ###############
p.deseq.cecum


sigtab.colon = res.colon[which(res.colon$padj < 0.05), ]
sigtab.colon = cbind(as(sigtab.colon, "data.frame"), as(tax_table(FS2.genus)[rownames(sigtab.colon), ], "matrix"))
format(sigtab.colon$padj, scientific = TRUE)
sigtab.colon$newp <- format(round(sigtab.colon$padj, digits = 3), scientific = TRUE)
sigtab.colon$Treatment <- ifelse(sigtab.colon$log2FoldChange >=0, "RPS", "Control")

sigtab.colon <- sigtab.colon[order(sigtab.colon$log2FoldChange),]

respire.feces.d21 <- sort(c(grep('Camp', sigtab.colon$Genus),
                            grep('Sutt', sigtab.colon$Genus),
                            grep('Strep', sigtab.colon$Genus),
                            #grep('Trep', sigtab.colon$Genus),
                            grep('Helic', sigtab.colon$Genus),
                            grep('Muci', sigtab.colon$Genus)))

sigtab.colon$Family <- as.character(sigtab.colon$Family)
sigtab.colon$Genus <- as.character(sigtab.colon$Genus)


resp.anno <- sigtab.colon[respire.feces.d21,]
resp.anno$x <- respire.feces.d21
resp.anno$lab <- paste(resp.anno$Family, resp.anno$Genus)

sigtab.colon$Family[respire.feces.d21] <- ''
sigtab.colon$Genus[respire.feces.d21] <- ''


p.deseq.colon <- ggplot(sigtab.colon, aes(x=reorder(rownames(sigtab.colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) + geom_text(aes(x=rownames(sigtab.colon), y=0, label = paste(Family, Genus, sep = ' ')), fontface='italic', size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + geom_text(data=resp.anno, aes(x=x, y=0, label=lab), fontface='bold.italic', size=4)+
  theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_blank(),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_blank(),
                                             axis.ticks.y = element_blank())+ ggtitle('Colonic mucosa')+ coord_flip()
#THIS ONE FOR FIG 5
p.deseq.colon


sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < 0.05), ]
sigtab.cec_cont_RNA = cbind(as(sigtab.cec_cont_RNA, "data.frame"), as(tax_table(FS2)[rownames(sigtab.cec_cont_RNA), ], "matrix"))
format(sigtab.cec_cont_RNA$padj, scientific = TRUE)
sigtab.cec_cont_RNA$newp <- format(round(sigtab.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
sigtab.cec_cont_RNA$Treatment <- ifelse(sigtab.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")

p.deseq.cec_cont_RNA <- ggplot(sigtab.cec_cont_RNA, aes(x=reorder(rownames(sigtab.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) + geom_text(aes(x=rownames(sigtab.cec_cont_RNA), y=0, label = Genus), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: cecal contents (RNA)')+ coord_flip()
p.deseq.cec_cont_RNA
# not good for fig, but still interesting

# Now doing everything over just at OTU level this time   #

# this is for supplementary table stuff and also for limiting network visualizations to only differentially abundant features #

# Reading in unsubsampled OTU table and taxonomy #
otu <- import_mothur(mothur_shared_file = './data/V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './data/V4.final.taxonomy')
meta <- read.table(file = './data/V4.metadata.txt', sep = '\t', header = TRUE)
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
otu.sigtab.D0 = res.D0[which(res.D0$padj < 0.1), ]

otu.sigtab.D0 = cbind(as(otu.sigtab.D0, "data.frame"), as(tax_table(FS2.D0)[rownames(otu.sigtab.D0), ], "matrix"))
format(otu.sigtab.D0$padj, scientific = TRUE)
otu.sigtab.D0$newp <- format(round(otu.sigtab.D0$padj, digits = 3), scientific = TRUE)
otu.sigtab.D0$Treatment <- ifelse(otu.sigtab.D0$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.D0$tissue <- 'D0_feces'
# write.table(otu.sigtab.D0, file = './OTUs/D0_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D12 #

FS2.D12 <- subset_samples(FS2, day == 12)

sample_sums(FS2.D12)
FS2.D12 <- prune_taxa(taxa_sums(FS2.D12) > 1, FS2.D12)

rowSums(FS2.D12@otu_table)

FS2.D12.De <- phyloseq_to_deseq2(FS2.D12, ~ design)

FS2.D12.De <- DESeq(FS2.D12.De, test = "Wald", fitType = "parametric")

res.D12 = results(FS2.D12.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
otu.sigtab.D12 = res.D12[which(res.D12$padj < 0.1), ]
otu.sigtab.D12 = cbind(as(otu.sigtab.D12, "data.frame"), as(tax_table(FS2.D12)[rownames(otu.sigtab.D12), ], "matrix"))
format(otu.sigtab.D12$padj, scientific = TRUE)
otu.sigtab.D12$newp <- format(round(otu.sigtab.D12$padj, digits = 3), scientific = TRUE)
otu.sigtab.D12$Treatment <- ifelse(otu.sigtab.D12$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.D12$tissue <- 'D12_feces'
otu.sigtab.D12$otu <- rownames(otu.sigtab.D12)

#write.table(otu.sigtab.D12, file = './OTUs/D12_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D15 #

FS2.D15 <- subset_samples(FS2, day == 15)
sample_sums(FS2.D15)
FS2.D15 <- prune_taxa(taxa_sums(FS2.D15) > 1, FS2.D15)

FS2.D15.De <- phyloseq_to_deseq2(FS2.D15, ~ design)

FS2.D15.De <- DESeq(FS2.D15.De, test = "Wald", fitType = "parametric")

res.D15 = results(FS2.D15.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
otu.sigtab.D15 = res.D15[which(res.D15$padj < 0.1), ]
otu.sigtab.D15 = cbind(as(otu.sigtab.D15, "data.frame"), as(tax_table(FS2.D15)[rownames(otu.sigtab.D15), ], "matrix"))
format(otu.sigtab.D15$padj, scientific = TRUE)
otu.sigtab.D15$newp <- format(round(otu.sigtab.D15$padj, digits = 3), scientific = TRUE)
otu.sigtab.D15$Treatment <- ifelse(otu.sigtab.D15$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.D15$tissue <- 'D15_feces'
otu.sigtab.D15$otu <- rownames(otu.sigtab.D15)

#write.table(otu.sigtab.D15, file = './OTUs/D15_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D19 #

FS2.D19 <- subset_samples(FS2, day == 19)
FS2.D19 <- prune_taxa(taxa_sums(FS2.D19) > 1, FS2.D19)

FS2.D19.De <- phyloseq_to_deseq2(FS2.D19, ~ design)

FS2.D19.De <- DESeq(FS2.D19.De, test = "Wald", fitType = "parametric")

res.D19 = results(FS2.D19.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
otu.sigtab.D19 = res.D19[which(res.D19$padj < 0.1), ]
otu.sigtab.D19 = cbind(as(otu.sigtab.D19, "data.frame"), as(tax_table(FS2.D19)[rownames(otu.sigtab.D19), ], "matrix"))
format(otu.sigtab.D19$padj, scientific = TRUE)
otu.sigtab.D19$newp <- format(round(otu.sigtab.D19$padj, digits = 3), scientific = TRUE)
otu.sigtab.D19$Treatment <- ifelse(otu.sigtab.D19$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.D19$tissue <- 'D19_feces'
otu.sigtab.D19$otu <- rownames(otu.sigtab.D19)

#write.table(otu.sigtab.D19, file = './OTUs/D19_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# D21 #

FS2.D21 <- subset_samples(FS2, day == 21 & tissue == 'feces')

FS2.D21 <- prune_taxa(taxa_sums(FS2.D21) > 1, FS2.D21)
FS2.D21.De <- phyloseq_to_deseq2(FS2.D21, ~ treatment)

FS2.D21.De <- DESeq(FS2.D21.De, test = "Wald", fitType = "parametric")

res.D21 = results(FS2.D21.De, cooksCutoff = FALSE, pAdjustMethod = 'BH')
otu.sigtab.D21 = res.D21[which(res.D21$padj < 0.05), ]
otu.sigtab.D21 = cbind(as(otu.sigtab.D21, "data.frame"), as(tax_table(FS2.D21)[rownames(otu.sigtab.D21), ], "matrix"))
format(otu.sigtab.D21$padj, scientific = TRUE)
otu.sigtab.D21$newp <- format(round(otu.sigtab.D21$padj, digits = 3), scientific = TRUE)
otu.sigtab.D21$Treatment <- ifelse(otu.sigtab.D21$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.D21$tissue <- 'D21_feces'
otu.sigtab.D21$otu <- rownames(otu.sigtab.D21)

#write.table(otu.sigtab.D21, file = './OTUs/D21_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

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

otu.sigtab.ileum = res.ileum[which(res.ileum$padj < 0.1), ]
otu.sigtab.ileum = cbind(as(otu.sigtab.ileum, "data.frame"), as(tax_table(FS2)[rownames(otu.sigtab.ileum), ], "matrix"))
format(otu.sigtab.ileum$padj, scientific = TRUE)
otu.sigtab.ileum$newp <- format(round(otu.sigtab.ileum$padj, digits = 3), scientific = TRUE)
otu.sigtab.ileum$Treatment <- ifelse(otu.sigtab.ileum$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.ileum$tissue <- 'ileum'
otu.sigtab.ileum$otu <- rownames(otu.sigtab.ileum)

otu.sigtab.cecum = res.cecum[which(res.cecum$padj < 0.1), ]
otu.sigtab.cecum = cbind(as(otu.sigtab.cecum, "data.frame"), as(tax_table(FS2)[rownames(otu.sigtab.cecum), ], "matrix"))
format(otu.sigtab.cecum$padj, scientific = TRUE)
otu.sigtab.cecum$newp <- format(round(otu.sigtab.cecum$padj, digits = 3), scientific = TRUE)
otu.sigtab.cecum$Treatment <- ifelse(otu.sigtab.cecum$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.cecum$tissue <- 'cecum'
otu.sigtab.cecum$otu <- rownames(otu.sigtab.cecum)

otu.sigtab.colon = res.colon[which(res.colon$padj < 0.1), ]
otu.sigtab.colon = cbind(as(otu.sigtab.colon, "data.frame"), as(tax_table(FS2)[rownames(otu.sigtab.colon), ], "matrix"))
format(otu.sigtab.colon$padj, scientific = TRUE)
otu.sigtab.colon$newp <- format(round(otu.sigtab.colon$padj, digits = 3), scientific = TRUE)
otu.sigtab.colon$Treatment <- ifelse(otu.sigtab.colon$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.colon$tissue <- 'colon'
otu.sigtab.colon$otu <- rownames(otu.sigtab.colon)

otu.sigtab.cec_cont_RNA = res.cec_cont_RNA[which(res.cec_cont_RNA$padj < 0.1), ]
otu.sigtab.cec_cont_RNA = cbind(as(otu.sigtab.cec_cont_RNA, "data.frame"), as(tax_table(FS2)[rownames(otu.sigtab.cec_cont_RNA), ], "matrix"))
format(otu.sigtab.cec_cont_RNA$padj, scientific = TRUE)
otu.sigtab.cec_cont_RNA$newp <- format(round(otu.sigtab.cec_cont_RNA$padj, digits = 3), scientific = TRUE)
otu.sigtab.cec_cont_RNA$Treatment <- ifelse(otu.sigtab.cec_cont_RNA$log2FoldChange >=0, "RPS", "Control")
otu.sigtab.cec_cont_RNA$tissue <- 'cec_cont_RNA'
otu.sigtab.cec_cont_RNA$otu <- rownames(otu.sigtab.cec_cont_RNA)

# write.table(otu.sigtab.cec_cont_RNA, file = './OTUs/cec_cont_RNA_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(otu.sigtab.colon, file = './OTUs/colon_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(otu.sigtab.cecum, file = './OTUs/cecum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')
# write.table(otu.sigtab.ileum, file = './OTUs/ileum_OTUs.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

dif_ab_16S <- rbind(otu.sigtab.D12,
                    otu.sigtab.D15,
                    otu.sigtab.D19,
                    otu.sigtab.D21,
                    otu.sigtab.cec_cont_RNA,
                    otu.sigtab.colon,
                    otu.sigtab.cecum,
                    otu.sigtab.ileum)#,
#otu.sigtab.wean)
############ COMMENTED HERE #################
#dif_ab_16S$name <- swap[dif_ab_16S$otu]
#dif_ab_16S <- dif_ab_16S[,c(1:14,17,15,16)]
write.table(dif_ab_16S, './data/dif_ab_16s.txt', row.names = TRUE, quote = FALSE, col.names = TRUE, sep = '\t')

# end OTU level stuff  #
# Below here everything is rarrefied to 4200 sequences #

FS2.4200 <- read.table("./data/V4.FS2.4200.shared", header = TRUE)
rownames(FS2.4200) <- FS2.4200$Group
meta <- read.table(file = './data/V4.metadata.txt', sep = '\t', header = TRUE)

meta <- meta[meta$group %in% rownames(FS2.4200),]
meta <- meta[meta$design != "FS2_cec_cont_DNA_21_RPS",]
meta <- meta[meta$day <22,]

FS2.4200 <- FS2.4200[rownames(FS2.4200) %in% meta$group,]

FS2.4200 <- FS2.4200[,-c(1,2,3)]


meta$group == rownames(FS2.4200)
# OTU 87 STUFF #

tax <- extract_mothur_tax('./data/V4.final.taxonomy')

timepo <- FS2.4200[,-1]/4200

timepo$group <- rownames(timepo)
timepo <- merge(meta, timepo, by = 'group')
timepo$day <- as.numeric(timepo$day)


timepo2 <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(tissue == 'feces', day <  22)
timepo2$otu2 <- factor(timepo2$otu)

tiss <- timepo %>% gather(otu, value, -(group:treatXday)) %>% filter(day <  22)


p.otu87.tiss <- tiss %>% filter(otu == 'Otu00087' & day ==21) %>%
  ggplot(aes(x=tissue, y=value))  +
  geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment), show.legend = TRUE) +
  #ggtitle('Relative proportion of OTU87 in each tissue a Day 21') +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('proportion of total community')

p.otu87.fec <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(85)]) %>%
  ggplot(aes(x=day, y=value))  +
  #geom_vline(xintercept = 12, colour='black')+
  geom_boxplot(aes(x=day, y=value, group=design, fill=treatment), show.legend = FALSE) +# ggtitle('Relative proportion of OTU87 in feces over time')
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+ ylab('proportion of total community')


###
# temp <- timepo2 %>% filter(otu %in% levels(timepo2$otu2)[c(85)])
# 
# temp <- tiss %>% filter(otu == 'Otu00087' & day ==21 & tissue == 'feces')
# # ORDINATION OVER TIME #

FS2.bray2 <- vegdist(FS2.4200, method = 'bray')

# dispersion stuff, this has its own metathing so hopefully later steps still work #

dispers <- betadisper(FS2.bray2, group = meta$design)
pdispers <- permutest(dispers, pairwise = TRUE, permutations = 1000)

dispersion_supp_16S <- as.data.frame(pdispers$pairwise$permuted)
dispersion_supp_16S$comp <- rownames(dispersion_supp_16S)

dispersion_supp_16S <- dispersion_supp_16S[grep('FS2_(.*)_[0-9]+_[A-Za-z]+-FS2_\\1_[0-9]+_[A-Za-z]+', dispersion_supp_16S$comp),]

colnames(dispersion_supp_16S)[1] <- '16S p value'
dispersion_supp_16S$comp <- gsub('FS2_', '', dispersion_supp_16S$comp)
dispersion_supp_16S$comp <- gsub('-', ' vs ', dispersion_supp_16S$comp)

# write.table(dispersion_supp_16S, file = './output/16S_dispersion_stats.txt', quote = FALSE, row.names = FALSE, sep = '\t')


dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group
metadisp <- merge(meta, dispersdf, by = 'group')

dispgroups <- summarise(group_by(metadisp, design), average_dist=mean(dispers.distances), sd_dist=sd(dispers.distances), se_dist=sd_dist/sqrt(n()))

dispgroups <- unique(inner_join(dispgroups, meta[,-c(1,2,5,7,9)]))

dispgroups %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) +
  geom_boxplot() + scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylim(c(.15,.7)) +
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")





metadisp %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) +
  geom_boxplot() + scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylim(c(.15,.7)) +
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")


fecaldisptime <-  dispgroups %>% filter(tissue == 'feces' & day < 22)

p.16s.disp.time <- ggplot(fecaldisptime, aes(x=day, y=average_dist, color = treatment)) +
  geom_vline(xintercept = 12, color='purple', size=.75, alpha=.7)+ geom_point(size = 2.5) +  geom_line(show.legend = FALSE) + ylim(c(.2,.5))+ xlim(c(0,22))+
  #ggtitle('Community Variability (Dispersion)',
  #       subtitle = "Vegan's betadisper(): how much variability is there in a group's community structure?") +
  ylab("Avg dist to group mean") + scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + 
  geom_errorbar(aes(ymin=average_dist-se_dist, ymax=average_dist+se_dist), width=0.2, size=1)+
  annotate(geom = 'text', x=21, y=.45, label='* p=0.001 *')
p.16s.disp.time

#############FIG 2 16S DISPER #############
# I think I should include this  figure in the supplement

# these ordinations contain all samples and are rarefied to 4200 seqs per sample
rownames(FS2.4200) == meta$group
rownames(meta) <- meta$group
rownames(meta) == rownames(FS2.4200)



FS2.ord <- NMDS_ellipse(metadata = meta,
                        OTU_table = FS2.4200,
                        grouping_set = 'design',
                        MDS_trymax = 1000, 
                        autotransform = FALSE)

FS2.ord[[3]]$species
stressplot(FS2.ord[[3]])




FS2.metanmds <- FS2.ord[[1]]
df_ell <- FS2.ord[[2]]

df_ell$treatment <- NA
df_ell$treatment[grep('RPS', df_ell$group)] <- 'RPS'
df_ell$treatment[grep('control', df_ell$group)] <- 'control'
df_ell$tissue <- gsub('FS2_(.*)_[0-9]+_.*', '\\1', df_ell$group)
df_ell$day <-gsub('FS2_(.*)_([0-9]+)_.*', '\\2', df_ell$group)


df_ell$treatmentXday <- paste(df_ell$treatment, df_ell$day, sep = ' ')
FS2.metanmds$day <- factor(FS2.metanmds$day)

########## 16S TIS ORD #######

groupxtiss <- FS2.metanmds %>% filter(day == 21) %>% select(tissue, day, treatment, centroidX, centroidY) %>% unique() 

p.tissord <- groupxtiss %>% ggplot(aes(x=centroidX, y=centroidY)) +
  geom_point(aes(color=treatment), size=8, show.legend = FALSE) +
  geom_path(data = subset(df_ell, day == 21), aes(x=NMDS1, y=NMDS2, group=group, color=treatment), size=1, show.legend = FALSE) +
  geom_text(data=groupxtiss, aes(x=centroidX, y=centroidY, label=tissue), show.legend = FALSE)+
  scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + labs(x='NMDS 1', y= 'NMDS 2')# + 
  #geom_point(data = filter(FS2.metanmds, day == 21), aes(x=MDS1, y=MDS2, color=treatment), alpha=.3)#+
#ggtitle('Bray-curtis based NMDS ordination of tissues at day 21')
p.tissord



# I really like this one.  I think I will keep it.
# do I need stats here?

#PERMANOVA with Adonis

PW.Adonis <- pairwise.adonis(FS2.4200,meta$design, sim.method="bray", p.adjust.m = "bonferroni")

fecVSfec <- grep('.*_feces_.* vs .*_feces_.*' ,PW.Adonis$pairs)
colVScol <- grep('.*_colon_.* vs .*_colon_.*' ,PW.Adonis$pairs)
cecVScec <- grep('.*_cecum_.* vs .*_cecum_.*' ,PW.Adonis$pairs)
rVSr <- grep('FS2_cec_cont_RNA_.* vs FS2_cec_cont_RNA_.*' ,PW.Adonis$pairs)
ilVSil <- grep('.*_ileum_.* vs .*_ileum_.*' ,PW.Adonis$pairs)

good <- PW.Adonis[c(fecVSfec, colVScol, cecVScec, rVSr, ilVSil),]

# FS2_feces_0_control vs FS2_feces_21_control 10.7868015 0.30141843 0.0001000     0.0435
# FS2_feces_21_control vs FS2_feces_0_RPS 11.5118295 0.31529040 0.0001000     0.0435


fecsametime <- grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", good$pairs)
ilsametime <- grep("FS2_ileum_([0-9]+)_control vs FS2_ileum_\\1_RPS", good$pairs)
colsametime <- grep("FS2_colon_([0-9]+)_control vs FS2_colon_\\1_RPS", good$pairs)
cecsametime <- grep("FS2_cecum_([0-9]+)_control vs FS2_cecum_\\1_RPS", good$pairs)
rsametime <- grep('RNA', good$pairs)
good2 <- good[c(fecsametime, ilsametime, colsametime, cecsametime,rsametime),]
write.table(good2,"./output/Adonis-Results.csv",sep=",")

just_poop <- good2[grep('feces', good2$pairs),]
poop_time <- just_poop[grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", just_poop$pairs),]
poop_time$day <- c(0,12,15,19,21)
poop_time$p.value <- round(poop_time$p.value, 4)

filter(poop_time, day <22)%>%
  ggplot(aes(x=day, y=F.Model)) +
  geom_line()+geom_point(size=2)+ geom_vline(xintercept = 12, color='purple') + geom_text(aes(y = F.Model + .2, label = paste('p =', p.value))) +
  ggtitle('Dissimilarity of fecal microbiota over time', subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint, how different are the two diets at each timepoint? ') +
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%')


#saving to combine with buts
adonfobuts <- filter(poop_time, day <22)
# this is nice...maybe add dispersion info on this too?

# phylum level grouping #

meta <- read.table(file = './data/V4.metadata.txt', sep = '\t', header = TRUE)

otu <- import_mothur(mothur_shared_file = './data/V4.final.shared')
taxo <- import_mothur(mothur_constaxonomy_file = './data/V4.final.taxonomy')


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

# should probably try a gather here instead
FS2.phylum.relabund <- melt(FS2.phylum.relabund, id.vars = c(1, 22:29) )

FS2.phylum.relabund$pig_num <- factor(FS2.phylum.relabund$pig_num)
FS2.phylum.relabund$day <- factor(FS2.phylum.relabund$day)

goodphy <- levels(FS2.phylum.relabund$variable)[c(1:8,13)]


#figure 1A
p.weaning.phyla <- filter(FS2.phylum.relabund, tissue == 'feces' & day %in% c(0,21) & variable %in% goodphy) %>%
  ggplot(aes(x=day, y=value, fill=treatment))+ geom_boxplot(aes(group=design), show.legend = FALSE) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + expand_limits(y=0) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  ggtitle('Weaning effects on major phyla in fecal microbiota') +
  ylab('proportion of total community')

# boom!  use it next to ordination showing D0 to D21 shift

D21phy <- levels(FS2.phylum.relabund$variable)[c(1,2,3,4,5,7,8,10,12)]

p.fecesphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'feces' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Feces') +
  ylab('proportion of total community')


p.cecphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'cecum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal mucosa') +
  ylab('proportion of total community')

p.colonphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'colon' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Colonic mucosa') +
  ylab('proportion of total community')

p.ileumphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'ileum' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Ileal mucosa') +
  ylab('proportion of total community')

p.cec_cont_RNAphyla.D21 <- filter(FS2.phylum.relabund, tissue == 'cec_cont_RNA' & day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal contents (RNA)') +
  ylab('proportion of total community')

FS2.phylum.relabund$tissue <- factor(FS2.phylum.relabund$tissue, levels = c('ileum', 'cecum', 'cec_cont_RNA','cec_cont_DNA', 'colon', 'feces'))

# These are good, but maybe its too much... can I concentrate this at all?

p.allphyla.D21 <- filter(FS2.phylum.relabund, tissue != 'cec_cont_DNA', day == 21 & variable %in% D21phy) %>%
  ggplot(aes(x=tissue, y=value, fill=treatment, group=design)) +
  geom_boxplot()+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  facet_wrap(~variable, shrink = TRUE, scales = 'free') +
  expand_limits(y=0) + ggtitle('Proportions of major phyla at day 21', subtitle = 'Cecal contents (RNA)') +
  ylab('proportion of total community')
# hmm... perhaps....

########### but amplicon  ############

# but taxonomy

x <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
x <- unlist(strsplit(x, split = ' '))

blastn <- read.table('./data/butreps_blastn.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)

blast <- read.table('./data/butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)

colnames(blast) <- x
colnames(blastn) <- x

blast$tit1 <- gsub('.* \\[(.*)\\]', '\\1', blast$salltitles)
blast$tit1 <- paste(blast$tit1, blast$pident, sep = ' ')  # change this if need to alter names of but otus

# read in mothur shared file #
shared2 <- read.table('./data/but2.shared', header = TRUE, stringsAsFactors = FALSE)
#hist(rowSums(shared2[,-c(1,2,3)]), breaks = 100)

#shared2[57,]
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
library(vegan)
otu2.r <- rrarefy(otu2, 2992)
otu3.rel <- otu2/rowSums(otu2)
otu2.r <- otu2.r[,colSums(otu2.r) > 0]

meta <- read.table('./data/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)


met <- data.frame(rownames(otu2.r))
colnames(met) <- 'group'

butmet <- read.table('./data/butmeta.txt', header = TRUE, stringsAsFactors = FALSE)
colnames(butmet)[1] <- 'group'

met <- merge(met, butmet, by = 'group', all = TRUE)
meta <- merge(meta, met, by = 'group', all = TRUE)

meta$tissue[is.na(meta$tissue)] <- 'feces'
meta$day.x[is.na(meta$day.x)] <- meta$day.y[is.na(meta$day.x)]
meta$pig_num[is.na(meta$pig_num)] <- meta$pigNumbers[is.na(meta$pig_num)]
meta$treatment.x[is.na(meta$treatment.x)] <- meta$treatment.y[is.na(meta$treatment.x)]
meta$experiment[is.na(meta$experiment)] <- 'FS2'
meta <- meta[,1:6]
meta <- meta[meta$experiment =='FS2',]

colnames(meta) <- c('group', 'experiment', 'tissue', 'day', 'pig_num', 'treatment')

#meta$tissue[is.na(meta$day)][2] <- 'cec_cont_DNA' # this is the old line...
meta$tissue[is.na(meta$day)] <- 'cec_cont_DNA'  # replaced it with this...

meta[is.na(meta$day),]
meta[is.na(meta$pig_num),]


meta$day[is.na(meta$day)] <- 21
meta$pig_num[is.na(meta$pig_num)] <- 75
meta$treatment[is.na(meta$treatment)] <- 'RPS'

meta[57,]
meta <- meta[-57,]  # why? this was to remove pig 97 from the data 
meta[57,]
meta$treatment <- gsub('RS', 'RPS', meta$treatment)
meta$day <- gsub('D(*)', '\\1', meta$day)

meta <- meta[meta$group %in% rownames(otu2.r),]

otu2.r <- otu2.r[rownames(otu2.r) %in% meta$group,]

meta$tissue[meta$tissue == 'cec_cont_DNA'] <- 'cec_cont_RNA'  # mislabeled? check this out

meta$design <- paste(meta$tissue, meta$day, meta$treatment, sep = '_')

meta$group == rownames(otu2.r)

meta <- meta[match(rownames(otu2.r), meta$group),]

meta$group == rownames(otu2.r)




#meta2 <- meta[meta$day !=0,]

#otu3.r <- otu2.r[rownames(otu2.r)%in%meta2$group,] # this is for NMDS ord, exlude D0

meta$group == rownames(otu2.r)
#meta2$group == rownames(otu3.r)

# for correlation
but2 <- data.frame(otu3.rel)
but2$group <- rownames(but2)
but2[57,]
but2 <- but2[-57,]
but2[57,]
write.table(but2, './data/but_shared_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Alpha Diversity #

# igraph diversity is conflicting with vegan diversity here...#

meta$invsimpson <- vegan::diversity(otu2.r, index = 'invsimpson')
meta$shannon <- vegan::diversity(otu2.r, index = 'shannon')

# Shannon

# filter(FecalTime.meta, day %in% c(0,12,15,19,21)) %>%
#  ggplot() + geom_boxplot(aes(day, shannon, fill = treatment)) +
#  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ggtitle('Alpha diversity of feces over time')

# USE THIS FIG FOR FECES ALPHA DIV #
p.butfeces.alphatime <- filter(meta, day %in% c(0,12,15,19,21) & tissue == 'feces') %>%
  ggplot() + geom_boxplot(aes(day, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ggtitle('Alpha diversity of feces over time')

p.buttissue.alpha <- filter(meta, day == 21) %>%
  ggplot() + geom_boxplot(aes(tissue, invsimpson, fill = treatment)) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ggtitle('Alpha diversity of tissues')

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
#FS2.bray3 <- vegdist(otu3.r, method = 'bray')

attributes(FS2.bray2)$Labels == meta$group # still good!

# dispersion stuff #

############ here for the dispersion stuff #################

dispers <- betadisper(FS2.bray2, group = meta$design)
dispers$distances

pdispers <- permutest(dispers, pairwise = TRUE)
anova(dispers)
pdispers$pairwise$observed

dispersdf <- data.frame(dispers$distances)
dispersdf$group <- rownames(dispersdf)
meta$group == dispersdf$group

dispersion_supp_but <- as.data.frame(pdispers$pairwise$permuted)
dispersion_supp_but$comp <- rownames(dispersion_supp_but)
dispersion_supp_but <- dispersion_supp_but[grep('(.*)_[0-9]+_[A-Za-z]+-\\1_[0-9]+_[A-Za-z]+', dispersion_supp_but$comp),]
colnames(dispersion_supp_but)[1] <- 'but p value'
#dispersion_supp_but$comp <- gsub('FS2_', '', dispersion_supp_but$comp)
dispersion_supp_but$comp <- gsub('-', ' vs ', dispersion_supp_but$comp)

#write.table(dispersion_supp_but, file = 'but_dispersion_stats.txt', quote = FALSE, row.names = FALSE, sep = '\t')

dispersion_supp_all <- merge(dispersion_supp_16S, dispersion_supp_but, by = 'comp')


dispersion_supp_all$comp

goodones <- c(grep('(.*)_([0-9]+)_(.*) vs .*_[0-9]+_\\3',dispersion_supp_all$comp),
  grep('(.*)_([0-9]+)_(.*) vs .*_\\2_.*',dispersion_supp_all$comp))

dispersion_supp_all <- dispersion_supp_all[goodones,]

write.table(dispersion_supp_all,'./output/dispersion_supp_all.txt' ,sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
but.metadisp <- merge(meta, dispersdf, by = 'group')

# this is a bad way to do this command.  SHould probably switch to pipes %>% %>% %>% %>% %>% %>% %>% %>% 

but.dispgroups <- summarise(group_by(.data = but.metadisp, design), average_dist=mean(dispers.distances))

but.dispgroups <- unique(inner_join(but.dispgroups, meta))
but.dispgroups$type <- 'but'


#############but dispersion fig 2 ###########
p.butfecal.dispersion.time <- but.metadisp %>% filter(tissue == 'feces' & day %in% c(0:22)) %>%
  ggplot(aes(x=day, y=dispers.distances, fill = treatment, group = design)) +
  geom_boxplot() + scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  ylab("distance to the group median") + ggtitle("Fecal beta diversity dispersion over time")

# maybe this is useful?  polish a little bit

but.fecaldisptime <-  but.dispgroups %>% filter(tissue == 'feces' & day < 22)
but.fecaldisptime$day <- as.numeric(but.fecaldisptime$day)

############ THIS FOR FIG2 BUT DISP ################

my_title <- expression(paste(italic("but"), " community dispersion"))
p.butfecal.dispersionAVG <- ggplot(but.fecaldisptime, aes(x=day, y=average_dist, color = treatment)) +
  geom_vline(xintercept = 12)+ geom_point(size = 2.5) +  geom_line() + ylim(c(.3,.6))+
  #ggtitle(my_title,
  #        subtitle = "Vegan's betadisper(): how much variability is there in a group's community structure?") +
  ylab("distance to the group mean") + scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS'))

#fecaldisptime$average_dist
#########Supplement #########
p.but.tissue.dispersion <- but.metadisp %>% filter(day == 21) %>%
  ggplot(aes(x=tissue, y=dispers.distances, group=design, fill=treatment)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ggtitle('Distance to group centroid')

# but NMDS Ordination #

rownames(meta) <- meta$group
rownames(otu2.r) == rownames(meta)


FS2.butord <- NMDS_ellipse(metadata = meta, OTU_table = otu2.r, grouping_set = 'design', 
                           autotransform = FALSE)

stressplot(FS2.butord[[3]])

FS2.but.metanmds <- FS2.butord[[1]]

df_ell <- FS2.butord[[2]]

df_ell$treatment <- NA
df_ell$treatment[grep('RPS', df_ell$group)] <- 'RPS'
df_ell$treatment[grep('control', df_ell$group)] <- 'control'
df_ell$tissue <- gsub('(.*)_[0-9]+_.*', '\\1', df_ell$group)
df_ell$day <-gsub('(.*)_([0-9]+)_.*', '\\2', df_ell$group)

# colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')
# FS2.but.metanmds <- merge(FS2.but.metanmds, NMDS.mean , by='design')

df_ell$treatmentXday <- paste(df_ell$treatment, df_ell$day, sep = ' ')
FS2.but.metanmds$day <- factor(FS2.but.metanmds$day)


colnames(FS2.but.metanmds)
groupxday <- unique(FS2.but.metanmds[,c(3,4,6,7,12,13)]) # this little dataframe is needed for annotations that dont look bad.
groupxday <- filter(groupxday, day %in% c(0:22) & tissue == 'feces')

p.but.time.ord <- FS2.but.metanmds %>% filter(tissue == 'feces' & day %in% c(0:21)) %>% ggplot(aes(x=centroidX, y=centroidY)) +
  geom_line(aes(group=treatment, color=treatment), size = 2, alpha=.7) +
  geom_point(aes(color=treatment), size=7) +
  geom_text(data=groupxday, aes(x=centroidX, y=centroidY, label=day))+
  scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS'))  + ggtitle('Fecal community structure over time', subtitle = 'Treatment group centroids') + ylab('NMDS2') + xlab('NMDS1')
p.but.time.ord


########### FIGURE 2 ###########
# the above one is nice for showing progression over time maybe

groupxtiss <- unique(FS2.but.metanmds[,c(3,4,6,7,12,13)]) # this little dataframe is needed for annotations that dont look bad.
groupxtiss <- filter(groupxtiss, day == 21)

p.but.tiss.ord <- FS2.but.metanmds %>% filter(day %in% c(21)) %>% ggplot(aes(x=centroidX, y=centroidY)) +
  geom_point(aes(color=treatment), size=8, show.legend = FALSE) +
  geom_path(data = subset(df_ell, day == 21), aes(x=NMDS1, y=NMDS2, group=group, color=treatment), size=1, show.legend = FALSE) +
  geom_text(data=groupxtiss, aes(x=centroidX, y=centroidY, label=tissue), show.legend = FALSE)+
  scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + labs(x='NMDS 1', y= 'NMDS 2') #+
#ggtitle('Bray-curtis based NMDS ordination of tissues at day 21')
p.but.tiss.ord

#########Figure 2 ###########
# I really like this one.  I think I will keep it.
# do I need stats here?

# PERMANOVA with Adonis #

PW.Adonis <- pairwise.adonis(otu2.r,meta$design, sim.method="bray", p.adjust.m = "bonferroni")

# write.table(PW.Adonis,"Adonis-Results.csv",sep=",")

fecVSfec <- grep('feces_.* vs feces_.*' ,PW.Adonis$pairs)
colVScol <- grep('colon_.* vs colon_.*' ,PW.Adonis$pairs)
cecVScec <- grep('cecum_.* vs cecum_.*' ,PW.Adonis$pairs)
rVSr <- grep('cec_cont_RNA_.* vs cec_cont_RNA_.*' ,PW.Adonis$pairs)
ilVSil <- grep('ileum_.* vs ileum_.*' ,PW.Adonis$pairs)

good <- PW.Adonis[c(fecVSfec, colVScol, cecVScec, rVSr, ilVSil),]

fecsametime <- grep("feces_([0-9]+)_control vs feces_\\1_RPS", good$pairs)
ilsametime <- grep("ileum_([0-9]+)_control vs ileum_\\1_RPS", good$pairs)
colsametime <- grep("colon_([0-9]+)_control vs colon_\\1_RPS", good$pairs)
cecsametime <- grep("cecum_([0-9]+)_control vs cecum_\\1_RPS", good$pairs)
rsametime <- grep('RNA', good$pairs)
pwadon.but <- good[c(fecsametime, ilsametime, colsametime, cecsametime,rsametime),]
write.table(pwadon.but,"./output//but.Adonis-Results.csv",sep=",")

#pwadon.but <- read.table("but.Adonis-Results.csv",sep=",")
just_poop.but <- pwadon.but[grep('feces_.* vs feces_.*', pwadon.but$pairs),]

# poop_time.but <- just_poop.but[grep("FS2_feces_([0-9]+)_control vs FS2_feces_\\1_RPS", just_poop.but$pairs),]
poop_time.but <- just_poop.but
poop_time.but$day <- c(0,21,12,15,19)

poop_time.but$p.value <- round(poop_time.but$p.value, 4)
poop_time.but$day
filter(poop_time.but, day < 22) %>%
  ggplot(aes(x=day, y=F.Model)) +
  geom_line()+geom_point(size=2)+ geom_vline(xintercept = 12, color='purple') + geom_text(aes(y = F.Model + .2, label = paste('p =', p.value))) +
  ggtitle('Dissimilarity of fecal microbiota over time', subtitle = 'PERMANOVA F statistic, control vs RPS at each timepoint, how different are the two diets at each timepoint? ') +
  labs(caption='Vertical line represents diet change: Lactose from 10% to 2.5%')

# this is nice...maybe add dispersion info on this too?
############## TRYING BOTH 16S AND BUT ON SAME GRAPH #################
butadon <- filter(poop_time.but, day < 22)
butadon$type <- 'but'
adonfobuts$type <- '16S'

butadon$p.adjusted <- p.adjust(butadon$p.value)
adonfobuts$p.adjusted <- p.adjust(adonfobuts$p.value)


alladon <- rbind(butadon, adonfobuts)
alladon$p.adjusted[alladon$p.adjusted > 0.05] <- NA

alladon$p.adjusted

alladon$lab <- ifelse(alladon$p.adjusted < 0.001, 'p < 0.001', 
                      ifelse(alladon$p.adjusted < 0.01, 'p < 0.01', 
                             ifelse(alladon$p.adjusted < 0.05, 'p < 0.05')))

#round(alladon$p.adjusted, 3)
#alladon$lab <- paste('p =', alladon$p.adjusted)
#alladon$lab[grep('NA', alladon$lab)] <- NA

p.alladon <- ggplot(alladon, aes(x=day, y=F.Model, group=type, color=type)) +
  geom_line()+geom_point(size=3)+ geom_vline(xintercept = 12, color='purple') +
  geom_text(aes(label=lab), nudge_y = .2, fontface = 'bold') +
  scale_color_brewer(palette = 'Set1') + ylab(label = 'PERMANOVA F stat.') 





##########################################################################
# Writing stuff for correlations #

meta$sample <- paste(meta$tissue, meta$day, meta$pig_num, sep = '_')

but.feces <- otu2.r[meta$tissue == 'feces' & meta$day %in% c(0,21),]
but.cec_cont_RNA <- otu2.r[meta$tissue == 'cec_cont_RNA',]
but.cecum <-  otu2.r[meta$tissue == 'cecum',]
but.colon <-  otu2.r[meta$tissue == 'colon',]
but.ileum <-  otu2.r[meta$tissue == 'ileum',]
rownames(but.cecum)

otu <- as.data.frame(otu3.rel)
otu$Group <- rownames(otu)

##################### TEST TEST ##################
#
#

#timepo$group <- rownames(timepo)
#rowSums(but2)
but3 <- otu3.rel
but3$group <- rownames(but3)
but3 <- merge(meta, but3, by = 'group')
but3$day <- as.numeric(but3$day)


but22 <- but3 %>% gather(otu, value, -(group:design)) %>% filter(tissue == 'feces', day <  22)
but22$otu2 <- factor(but22$otu)

tiss <- but3 %>% gather(otu, value, -(group:design)) %>% filter(day <  22)


p.but.otu67.tiss <- tiss %>% filter(otu == 'Otu067' & day ==21) %>%
  ggplot(aes(x=tissue, y=value))  +
  geom_boxplot(aes(x=tissue, y=value, group=design, fill=treatment), show.legend = TRUE) +
  #ggtitle('Relative proportion of OTU87 in each tissue a Day 21') +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('proportion of total community')

p.but.otu67.fec <- but22 %>% filter(otu == 'Otu067') %>%
  ggplot(aes(x=day, y=value))  +
  #geom_vline(xintercept = 12, colour='black')+
  geom_boxplot(aes(x=day, y=value, group=design, fill=treatment), show.legend = FALSE) +# ggtitle('Relative proportion of OTU87 in feces over time')
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+ ylab('proportion of total community')

#switchd <- but22 %>% filter(otu == 'Otu067' & day == 21)

############## end test


write.table(meta, './data/but_meta_forcorr.txt', col.names = TRUE, row.names = FALSE, quote = FALSE)
# write.table(otu, './data/but.all.shared', col.names = TRUE, row.names = TRUE, quote = FALSE)

# Deseq2 stuff #

otu <- import_mothur(mothur_shared_file = './data/but2.shared')
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

p.but.deseq.D0 <- ggplot(sigtab.but.D0, aes(x=reorder(rownames(sigtab.but.D0), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D0), y=0, label = name), size=3)+ labs(x="otu")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant genera: feces')+ coord_flip()
p.but.deseq.D0

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

p.but.deseq.D12 <- ggplot(sigtab.but.D12, aes(x=reorder(rownames(sigtab.but.D12), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D12), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D12 feces')+ coord_flip()
p.but.deseq.D12

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

p.but.deseq.D15 <- ggplot(sigtab.but.D15, aes(x=reorder(rownames(sigtab.but.D15), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D15), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D15 feces')+ coord_flip()
p.but.deseq.D15

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

p.but.deseq.D19 <- ggplot(sigtab.but.D19, aes(x=reorder(rownames(sigtab.but.D19), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D19), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D19 feces')+ coord_flip()
p.but.deseq.D19

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

# now a lot of garbage to make sure the italics are correct...

sigtab.but.D21$PID <- sub('.* ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21$name)
sigtab.but.D21$ital <- sub('(.*) ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21$name)
sigtab.but.D21$ital <- sub('uncultured ', '', sigtab.but.D21$ital)
sigtab.but.D21$ital <- sub('\\[', '', sigtab.but.D21$ital)
sigtab.but.D21$ital <- sub('\\]', '', sigtab.but.D21$ital)

sigtab.but.D21$ital2 <- sub('(.*) (sp\\. ?.*)', '\\1', sigtab.but.D21$ital)
sigtab.but.D21$nonital <- sub('(.*) (sp\\. ?.*)', '\\2', sigtab.but.D21$ital)
sigtab.but.D21$nonital[sigtab.but.D21$ital2 == sigtab.but.D21$nonital] <- ''


sigtab.but.D21$ital2[3] <- 'Firmicutes'
sigtab.but.D21$nonital[3] <- 'bacterium CAG:176'

sigtab.but.D21$ital2[6] <- 'Coprococcus catus'
sigtab.but.D21$nonital[6] <- 'GD/7'

sigtab.but.D21$ital2[9] <- 'Lachnospiraceae'
sigtab.but.D21$nonital[9] <- 'bacterium V9D3004'

sigtab.but.D21$ital2[13] <- 'Mogibacterium timidum'
sigtab.but.D21$nonital[13] <- 'ATCC 33093'

sigtab.but.D21$ital2[17] <- 'Peptostreptococcaceae'
sigtab.but.D21$nonital[17] <- 'bacterium oral taxon 113 str. W5053'

sigtab.but.D21$ital2[26] <- 'Firmicutes'
sigtab.but.D21$nonital[26] <- 'bacterium CAG:176'

sigtab.but.D21$ital2[28] <- 'Firmicutes'
sigtab.but.D21$nonital[28] <- 'bacterium CAG:176'

sigtab.but.D21$ital2[30] <- 'Lachnospiraceae'
sigtab.but.D21$nonital[30] <- 'bacterium P6B14'

sigtab.but.D21$ital2[31] <- 'Clostridiales'
sigtab.but.D21$nonital[31] <- 'bacterium CoAT_30-4c'

# end garbage #

# need ggscinames here... and the grid package

sigtab.but.D21$imp <- sigtab.but.D21$otu %in% c('Otu007','Otu034', 'Otu045','Otu067')

p.but.deseq.D21 <- ggplot(sigtab.but.D21, aes(x=reorder(rownames(sigtab.but.D21), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text_sciname(aes(x=rownames(sigtab.but.D21), y=0, sci=ital2, nonsci=nonital, important=imp), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21 feces')+ coord_flip()
p.but.deseq.D21

library(grid)
library(ggscinames)

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

# italics stuff...

sigtab.but.D21.ileum$PID <- sub('.* ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.ileum$name)
sigtab.but.D21.ileum$ital <- sub('(.*) ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.ileum$name)
sigtab.but.D21.ileum$ital <- sub('uncultured ', '', sigtab.but.D21.ileum$ital)
sigtab.but.D21.ileum$ital <- sub('\\[', '', sigtab.but.D21.ileum$ital)
sigtab.but.D21.ileum$ital <- sub('\\]', '', sigtab.but.D21.ileum$ital)

sigtab.but.D21.ileum$ital2 <- sub('(.*) (sp\\. ?.*)', '\\1', sigtab.but.D21.ileum$ital)
sigtab.but.D21.ileum$nonital <- sub('(.*) (sp\\. ?.*)', '\\2', sigtab.but.D21.ileum$ital)
sigtab.but.D21.ileum$nonital[sigtab.but.D21.ileum$ital2 == sigtab.but.D21.ileum$nonital] <- ''
# sigtab.but.D21.ileum$ital2[10] <- 'Clostridiales'
# sigtab.but.D21.ileum$nonital[10] <- 'bacterium CoAT_30-4c'

# end garbage #


sigtab.but.D21.ileum$imp <- sigtab.but.D21.ileum$otu %in% c('Otu007','Otu034', 'Otu045','Otu067')

sigtab.but.D21.ileum$PID <- paste('(', round(as.numeric(sigtab.but.D21.ileum$PID),0), ')', sep='')

p.but.deseq.D21.ileum <- ggplot(sigtab.but.D21.ileum, aes(x=reorder(rownames(sigtab.but.D21.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) +
  geom_text_sciname(aes(x=rownames(sigtab.but.D21.ileum), y=0, sci = ital2, nonsci=paste(nonital, PID, sep = ' '), important=imp), size=3)+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
        axis.title.x=element_text(size = 10), axis.title.y=element_blank())+ ggtitle('Ileal mucosa')+ coord_flip()

p.but.deseq.D21.ileum


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

# itallics stuff

sigtab.but.D21.cecum$PID <- sub('.* ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.cecum$name)
sigtab.but.D21.cecum$ital <- sub('(.*) ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.cecum$name)
sigtab.but.D21.cecum$ital <- sub('uncultured ', '', sigtab.but.D21.cecum$ital)
sigtab.but.D21.cecum$ital <- sub('\\[', '', sigtab.but.D21.cecum$ital)
sigtab.but.D21.cecum$ital <- sub('\\]', '', sigtab.but.D21.cecum$ital)

sigtab.but.D21.cecum$ital2 <- sub('(.*) (sp\\. .*)', '\\1', sigtab.but.D21.cecum$ital)
sigtab.but.D21.cecum$nonital <- sub('(.*) (sp\\. .*)', '\\2', sigtab.but.D21.cecum$ital)
sigtab.but.D21.cecum$nonital[sigtab.but.D21.cecum$ital2 == sigtab.but.D21.cecum$nonital] <- ''
sigtab.but.D21.cecum$ital2[10] <- 'Clostridiales'
sigtab.but.D21.cecum$nonital[10] <- 'bacterium CoAT_30-4c'

# end garbage #

sigtab.but.D21.cecum$imp <- sigtab.but.D21.cecum$otu %in% c('Otu007','Otu034', 'Otu045','Otu067')

sigtab.but.D21.cecum$PID <- paste('(', round(as.numeric(sigtab.but.D21.cecum$PID),0), ')', sep='')

p.but.deseq.D21.cecum <- ggplot(sigtab.but.D21.cecum, aes(x=reorder(rownames(sigtab.but.D21.cecum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) +
  geom_text_sciname(aes(x=rownames(sigtab.but.D21.cecum), y=0, sci = ital2, nonsci=paste(nonital, PID, sep = ' '), important=imp), size=3)+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_blank())+ ggtitle('Cecal mucosa')+ coord_flip()

p.but.deseq.D21.cecum


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

p.but.deseq.D21.colon <- ggplot(sigtab.but.D21.colon, aes(x=reorder(rownames(sigtab.but.D21.colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.but.D21.colon), y=0, label = name), size=3)+ labs(x="Genus")+
  scale_fill_brewer(palette = "Dark2") + theme(axis.text.x=element_text(color = 'black', size = 12),
                                             axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
                                             axis.title.x=element_text(size = 10),
                                             axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant but OTUs: D21.colon')+ coord_flip()
p.but.deseq.D21.colon

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






fec_cec <- sigtab.but.D21[sigtab.but.D21$otu %in% sigtab.but.D21.cecum$otu,]$otu

fec_cec_il <- fec_cec[fec_cec %in% sigtab.but.D21.ileum$otu]


# THESE ARE THE OTUS TO BOLD!!!!!! #######
all_tiss <- fec_cec_il[fec_cec_il %in% sigtab.but.D21.cec_cont_RNA$otu]


# bold_but <- c(grep('Otu007', sigtab.but.D21.cec_cont_RNA$otu),
#                        grep('Otu034', sigtab.but.D21.cec_cont_RNA$otu),
#                        grep('Otu045', sigtab.but.D21.cec_cont_RNA$otu),
#                        grep('Otu067', sigtab.but.D21.cec_cont_RNA$otu))

# all this garbage is to get the proper italics in the figure #
sigtab.but.D21.cec_cont_RNA$PID <- sub('.* ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.cec_cont_RNA$name)
sigtab.but.D21.cec_cont_RNA$ital <- sub('(.*) ([0-9]*\\.?[0-9]+)','\\1' ,sigtab.but.D21.cec_cont_RNA$name)
sigtab.but.D21.cec_cont_RNA$ital <- sub('uncultured ', '', sigtab.but.D21.cec_cont_RNA$ital)
sigtab.but.D21.cec_cont_RNA$ital <- sub('\\[', '', sigtab.but.D21.cec_cont_RNA$ital)
sigtab.but.D21.cec_cont_RNA$ital <- sub('\\]', '', sigtab.but.D21.cec_cont_RNA$ital)
sigtab.but.D21.cec_cont_RNA$ital2 <- sub('(.*) sp\\. CAG:[0-9]+', '\\1', sigtab.but.D21.cec_cont_RNA$ital)
sigtab.but.D21.cec_cont_RNA$nonital <- sub('(.*) (sp\\. CAG:[0-9]+)', '\\2', sigtab.but.D21.cec_cont_RNA$ital)
sigtab.but.D21.cec_cont_RNA$nonital[sigtab.but.D21.cec_cont_RNA$ital2 == sigtab.but.D21.cec_cont_RNA$nonital] <- ''
sigtab.but.D21.cec_cont_RNA$ital2[15] <- 'Eubacterium hallii'
sigtab.but.D21.cec_cont_RNA$nonital[15] <- 'DSM 3353'
sigtab.but.D21.cec_cont_RNA$ital2[12] <- 'Clostridium'
sigtab.but.D21.cec_cont_RNA$nonital[12] <- 'sp.'
#sigtab.but.D21.cec_cont_RNA$lab <- paste('italic(', sigtab.but.D21.cec_cont_RNA$ital2, ') plain(',sigtab.but.D21.cec_cont_RNA$nonital, ' ', sigtab.but.D21.cec_cont_RNA$PID, ')',sep = '')

# end garbage #
sigtab.but.D21.cec_cont_RNA$PID <- paste('(', round(as.numeric(sigtab.but.D21.cec_cont_RNA$PID),0), ')', sep='')




sigtab.but.D21.cec_cont_RNA$imp <- sigtab.but.D21.cec_cont_RNA$otu %in% c('Otu007','Otu034', 'Otu045','Otu067')

p.but.deseq.D21.cec_cont_RNA <- ggplot(sigtab.but.D21.cec_cont_RNA, aes(x=reorder(rownames(sigtab.but.D21.cec_cont_RNA), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity', show.legend = FALSE) +
  geom_text_sciname(aes(x=rownames(sigtab.but.D21.cec_cont_RNA), y=0, sci = ital2, nonsci=paste(nonital, PID, sep = ' '), important=imp), size=3)+
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) +
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_blank())+ ggtitle('Cecal contents (RNA)')+ coord_flip()
p.but.deseq.D21.cec_cont_RNA


#####



######### use this for butfigs ######

dif_ab_but <- rbind(sigtab.but.D12,
                    sigtab.but.D15,
                    sigtab.but.D19,
                    sigtab.but.D21[,-c(12:16)],
                    sigtab.but.D21.ileum[,-c(12:16)],
                    sigtab.but.D21.cecum[,-c(12:16)],
                    sigtab.but.D21.colon[,-c(12:16)],
                    sigtab.but.D21.cec_cont_RNA[,-c(12:16)])


write.table(dif_ab_but, './data/dif_ab_but.txt', row.names = TRUE, col.names = TRUE, quote = FALSE, sep='\t')

######################### Flow cytometry analysis  ###################

meta <- read.table('./data/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)

CD3neg <- read.table(file = './data/CD3-.txt.csv', header = TRUE, sep = '\t', as.is = TRUE, check.names = FALSE)

CD3pos <- read.table(file = './data/CD3+.txt.csv', header = TRUE, sep = '\t', as.is = TRUE, check.names = FALSE)

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

write.table(all.flow, file = './data/FS2_all_flow.txt', sep = '\t', row.names = FALSE)

#

#flow_counts <- read.table('FS2_all_flow.txt', header = TRUE, check.names = FALSE, sep = '\t', as.is = TRUE, comment.char = '')

meta <- read.table('./data/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)
colnames(all.flow)
flow_counts <- all.flow[,-c(2:17)]
colnames(flow_counts)

rowSums(flow_counts[,-c(1,18,19,20)]) # this tells us the total number of live single cells sorted for each sample.



flow_counts[,-c(1,18,19,20)] <- flow_counts[,-c(1,18,19,20)] / rowSums(flow_counts[,-c(1,18,19,20)]) # converts to relative abundance

# alpha div  #
meta.flow <- flow_counts


# This is for CD3 populations supplement #

flowgath <- flow_counts %>% 
  gather(key = cell_type, value = value, -Sample, -(tissue:treatment)) %>% 
  group_by(tissue, cell_type) %>% 
  summarise(av_abund=(sum(value)/n())*100) 
flowgath$cell_type <- gsub('/', '', flowgath$cell_type)
CD3supp <- flowgath %>% mutate(av_abund = round(av_abund, 1)) %>%  spread(cell_type, av_abund)
CD3supp <- as.data.frame(t(CD3supp))
colnames(CD3supp)
CD3supp$cell_type <- rownames(CD3supp)
colnames(CD3supp) <- c('cecum', 'lymph-node', 'PBMC', 'cell type')
CD3supp <- CD3supp[-1,]
rownames(CD3supp) <- NULL
CD3supp <- CD3supp[,c(4,1,2,3)]

write.table(CD3supp, './output/CD3supp.txt', sep = '\t', quote = FALSE, row.names = FALSE)
#end supplement stuff

### T cell Ecology stuff 

H <- vegan::diversity(flow_counts[,-c(1,18,19,20)])
J <- H/log(specnumber(flow_counts[,-c(1,18,19,20)])) # Pielou's J (evenness)
meta.flow$evenness <- J
meta.flow$Shannon <- H

# maybe this one for fig 7 as well ##########

meta.flow$tissue[meta.flow$tissue =='cec'] <- 'cecum'
meta.flow$tissue[meta.flow$tissue =='ln'] <- 'lymph node'

p.flow.even <- ggplot(meta.flow, aes(tissue, evenness, fill = treatment)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + geom_point(shape=21, size=2, stroke=1, position = position_jitterdodge(jitter.width = .2), show.legend = FALSE) + 
  annotate(x=.8, y=.68, xend=1.1, yend=.68, geom='segment') + annotate(x=.95, y=.685, label='*p=0.03*', geom = 'text') + theme(axis.text.x = element_text(c('cecum', 'lymph node', 'PBMC')))
p.flow.even
wilcox.test(filter(meta.flow, tissue=='cecum')$evenness~filter(meta.flow, tissue=='cecum')$treatment)

#ggplot(meta.flow) + geom_boxplot(aes(tissue, Shannon, fill = treatment))

#

barflow <- flow_counts
#barflow$treatment <- c(rep( c("control", "RPS"), each=7))
rownames(barflow) <- barflow$Sample
#barflow <- barflow[,-1]
barflow <- barflow %>% gather(attribute, value, -c(Sample, pig_num, treatment, tissue))

#barflow$pig_num <- factor(barflow$pig_num)

# some neat stacked bar graphs that didnt make the paper

ggplot(data=barflow, aes(x=pig_num, y=value, fill=attribute)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Proportion of live single cells') +
  facet_grid(~tissue) + geom_segment(x=7.5, y=0, xend=7.5, yend=1)

barflow[barflow$tissue =='PBMC'& barflow$treat == 'RPS',]$value <- barflow[barflow$tissue =='PBMC'& barflow$treat == 'RPS',]$value * 7/6 #B/C missing 1 obs in RPS treat

ggplot(data=barflow, aes(x=treatment, y=value, fill=attribute)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Proportion of CD3+ cells') +
  facet_grid(~tissue)

######### I REALLY LIKE THE ABOVE ONE #######

rownames(flow_counts) <- flow_counts$Sample

meta.flow$treatment[meta.flow$treatment == 'control'] <- 'CON'

rownames(flow_counts) == rownames(meta.flow)

#flow_counts <- flow_counts[,-c(1,18:20)]
meta.flow$tissueXtreat <- paste(meta.flow$tissue, meta.flow$treatment, sep = ' ')

flow.nmds <- NMDS_ellipse(metadata = meta.flow, OTU_table = flow_counts[,-c(1,18:20)], grouping_set = 'tissueXtreat')

flow.nmds[[1]]
flow.nmds[[2]]
flow.mds <- flow.nmds[[3]]

stressplot(flow.mds) # pretty good stressplot here....
nmds.stress <- flow.nmds[[3]]$stress


#meta <- flow_counts[,c(1,34,35,36)]

flow.nmds[[2]]$treatment <- gsub('([A-Za-z ]+) ([A-Z]+)', '\\2', flow.nmds[[2]]$group)
flow.nmds[[2]]$tissue <- gsub('([A-Za-z ]+) ([A-Z]+)', '\\1', flow.nmds[[2]]$group)

colnames(flow.nmds[[1]])[23] <- 'type'



p.floword.all <- ggplot(flow.nmds[[1]], aes(x=MDS1, y = MDS2)) +
  #geom_text_repel(data=spp.scrs, aes(MDS1, MDS2, label = cell_type), size=3, alpha=.7)+
  #geom_segment(data = spp.scrs, aes(x=0, xend=spp.scrs$MDS1, y=0, yend=spp.scrs$MDS2), alpha=.5)+
  geom_point(data=flow.nmds[[1]], aes(color=type ), size=3) +
  geom_path(data = flow.nmds[[2]], aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
  #ggtitle('Cecal T-cell community similarity: NMDS ordination using Bray-Curtis distances', subtitle = 'Treatments separated within Tissues') +
  geom_segment(data = flow.nmds[[1]],
               aes(x=flow.nmds[[1]]$MDS1,
                   xend=flow.nmds[[1]]$centroidX,
                   y=flow.nmds[[1]]$MDS2,
                   yend=flow.nmds[[1]]$centroidY,
                   color=type ),size=1) + scale_color_brewer(palette="Dark2")


p.floword.all


#  boxplots  #

fin <- data.frame()
# this for loop does wilcoxon tests for each cell type and saves the results into a new column that I can use to label??

barflow <- flow_counts
#barflow$treatment <- c(rep( c("control", "RPS"), each=7))
rownames(barflow) <- barflow$Sample
barflow <- barflow[,-1]
barflow <- barflow %>% gather(attribute, value, -c(pig_num, treatment, tissue))

barflow$pig_num <- factor(barflow$pig_num)


barflow$attribute <- factor(barflow$attribute)

barflow$tissue <- factor(barflow$tissue)
barflow$attribute <- gsub('CD3\\+/', '' ,barflow$attribute)

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


fin$pvalue2 <- round(fin$pvalue, 3)
fin$pvalue2 <- factor(paste('wilcox pvalue = ', fin$pvalue2, sep = ''))


fin$value <- fin$value * 100

# just the significantly different populations in all tissues
p.flow.allsig <- filter(fin, pvalue < 0.05) %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+
  facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+
  labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))


p.flow.allsig

fin$attribute <- gsub('/', ' ', fin$attribute)
# just the significantly different populations only in the cecum

p.flow.cecsig <- filter(fin, pvalue < 0.1 & tissue == 'cec') %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = 'total'))+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=1, width = .2)+
  facet_wrap(~attribute+pvalue2, scales = 'free')+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+
  labs(x="treatment", y= "Percent CD3+ cells") + theme(strip.text = element_text(size = 9, family ='bold')) + 
  ggtitle('Differentially abundant T-cell types in cecal tissues')

# THIS ONE ###
p.flow.cecsig


#####################  VFAs  ####################
vfa.cec <- read.table('./data/cec_redo_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfa.fec <- read.table('./data/fecalplasma_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfas2 <- read.table('./data/vfas_orig.txt', header = TRUE, sep = '\t') %>% filter(location == 'cecum' & day ==21)

#colnames(vfas)
colnames(vfas2)[c(2,5, 13)] <- c('tissue', 'pig_num', 'lactate')

vfas <- rbind(vfa.cec, vfa.fec)
vfas[,2:16] <- vfas[,2:16]*3
vfas$BCFA <- vfas$isobutyrate + vfas$isovalerate
vfas$total <- rowSums(vfas[,2:16])
colnames(vfas)[8] <- 'lactate'


vfas$sample <- paste(vfas$tissue, 21, vfas$pig_num, sep = '_')
vfas$treatment <- ifelse(vfas$pig_num %in% c(67,68,69,70,71,72,73,81,82,83,84,85,86,87), 'control', 'RPS')
vfas$design <- paste(vfas$tissue, vfas$treatment, sep = '_')

vfas2 <- vfas2[,colnames(vfas2) %in% colnames(vfas)]
vfas <- vfas[,colnames(vfas) %in% colnames(vfas2)]
vfas <- vfas[,match(colnames(vfas2), colnames(vfas))]

colnames(vfas) == colnames(vfas2)

vfas.final <- rbind(vfas, vfas2)
vfas.cec <- vfas.final %>% filter(tissue=='cecum') %>% select(c(4:20)) %>% group_by(pig_num) %>% summarise_all(funs(mean))
vfas.cec$tissue <- 'cecum'
vfas.cec[,17]
vfas.cec <- vfas.cec[,-17]
colnames(vfa.fec)[8] <- 'lactate'

colnames(vfa.fec) == colnames(vfas.cec)

vfas.final <- rbind(vfas.cec, vfa.fec)
vfas.final$BCFA <- vfas.final$isobutyrate + vfas.final$isovalerate
vfas.final$total <- rowSums(vfas.final[,c(2:11, 12:16)])
# colnames(vfas.final)[c(2:11, 12:16)]
#
colnames(vfas.final)[c(7,17)]

#vfas.final <- vfas.final[,-c(7,12)]
vfas.melt <- melt(data = vfas.final, id.vars = c(1,17))
control <- c(67:73, 81:87)
vfas.melt$treatment <- ifelse(vfas.melt$pig_num %in% control, 'control', 'RPS')
vfas.melt$design <- paste(vfas.melt$tissue, vfas.melt$treatment)
vfas.melt[vfas.melt$tissue == 'feces',]$value <- vfas.melt[vfas.melt$tissue == 'feces',]$value*3
#vfas.melt  %>% ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
#vfas.melt %>% filter(tissue == 'portal') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')
#vfas.melt %>% filter(tissue == 'cecum') %>% ggplot(aes(x=treatment, y=value, group=treatment)) + geom_boxplot() + facet_wrap(~variable, scales = 'free')

vfas.melt %>% filter(variable %in% c('propionate', 'butyrate', 'lactate','valerate', 'caproate', 'total') & tissue != 'portal') %>%
  ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total'),outlier.shape = NA) + facet_wrap(~variable, scales = 'free', nrow = 2, ncol = 3) + geom_point(shape=21, size=2, stroke=1, position = position_jitterdodge(jitter.width = .2)) + 
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('Concentration (mM)') + ggtitle('Cecal and fecal SCFA concentrations')



p.vfas.final <- vfas.melt %>% filter(variable %in% c('propionate', 'butyrate', 'lactate','valerate', 'caproate', 'total') & tissue != 'portal') %>%
  ggplot(aes(x=tissue, y=value, group=design, fill=treatment)) +
  geom_boxplot(position = position_dodge2(preserve = 'total'),outlier.shape = NA) + facet_wrap(~variable, scales = 'free', nrow = 2, ncol = 3) + geom_point(shape=21, size=2, stroke=1, position = position_jitterdodge(jitter.width = .2)) + 
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('Concentration (mM)') + ggtitle('Cecal and fecal SCFA concentrations')



p.vfas.final

#VFA STATISTICAL TESTS HERE #

# this for loop does wilcoxon tests for each cell type and saves the results into a new column that I can use to label??

vfas.melt$variable <- factor(vfas.melt$variable)

vfas.melt$tissue <- factor(vfas.melt$tissue)

vfas.melt$tisXatt <- paste(vfas.melt$tissue, vfas.melt$variable, sep = '_')

vfas.melt$tisXatt <- factor(vfas.melt$tisXatt)

vfas.WC <- data.frame()

for (lev in levels(vfas.melt$tisXatt)){
  temp <- vfas.melt[which(vfas.melt$tisXatt == lev),]
  twilc <- wilcox.test(temp$value ~ temp$treatment)
  temp$pvalue <- twilc$p.value
  vfas.WC <- rbind(vfas.WC,temp)
}


vfas.WC$pvalue2 <- round(vfas.WC$pvalue, 4)
vfas.WC$pvalue2 <- factor(paste('wilcox pvalue = ', vfas.WC$pvalue2, sep = ''))

# All
p.vfas.allsig <- filter(vfas.WC, pvalue < 1) %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
  geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+
  facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+
  labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))


p.vfas.allsig


# just the (almost) significantly different populations only in the cecum

# p.vfas.cecsig <- filter(vfas.WC, pvalue < 0.1 & tissue == 'cecum') %>% ggplot(aes(x=treatment, y=value, fill=treatment))+geom_boxplot()+
#   geom_jitter(shape=21, aes(fill=treatment), size =2, stroke=.75, width = .1)+ geom_text(aes(label = pig_num))+
#   facet_wrap(~tisXatt+pvalue2, scales = 'free')+ scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))+
#   labs(x="treatment", y= "Percent Live Single cells") + theme(strip.text = element_text(size = 9))
# 

# p.vfas.cecsig

############### BUTNET ###############

# but taxonomy #

# column names for the blast results
xx <- "qacc sacc sallseqid staxids sblastnames salltitles evalue qstart qend sstart send sscinames pident length"
xx <- unlist(strsplit(xx, split = ' '))

# reading in blast results
blast <- read.table('./data/butreps_blastx2.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
shared2 <- read.table('./data/but2.shared', header = TRUE, stringsAsFactors = FALSE)

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

butswap[blast$otu]  # this allows me to swap out the OTU names with a more informative name, kinda like a python dict


# Read in data #

tax <- extract_mothur_tax('./data/V4.final.taxonomy')
tax <- tax[,-c(2,3)]
swap <- otu_tax_labels(tax)

meta <- read.table('./data/16S_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)

otu <- read.table('./data/16S_shared_forcorr.txt', header = TRUE) %>% filter(group %in% meta$group)

rownames(otu) <- meta$sample

# differential abundance data for but and 16S otus, previously calculated by DeSeq2



dif_ab_16S <- read.table('./data/dif_ab_16s.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
dif_ab_but <- read.table('./data/dif_ab_but.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# HERE IS ENRICH 

dif_ab_16S$node <- swap[dif_ab_16S$otu]
dif_ab_but$node <- butswap[dif_ab_but$otu]

# limiting the differential abundance dataframes to only have the same columns
dif_ab_16S <- dif_ab_16S[,which(colnames(dif_ab_16S) %in% colnames(dif_ab_but))]
dif_ab_but <- dif_ab_but[,which(colnames(dif_ab_but) %in% colnames(dif_ab_16S))]

enrich <- rbind(dif_ab_16S, dif_ab_but) # This data frame has all the enrichment data in it


enrich$node

# but data #


but <- read.table('./data/but_shared_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
but$group
butmeta <- read.table('./data/but_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
but$group == butmeta$group
butmeta$sample <- paste(butmeta$tissue, butmeta$day, butmeta$pig_num, sep = '_')
rownames(but) <- butmeta$sample

# getting everything in order #
# Only feces used in this correlation

meta$tissue

# excluding day 0
meta <- meta %>% filter(tissue == 'feces' & day != 0)
butmeta <- butmeta %>% filter(tissue == 'feces' & day != 0)



meta <- filter(meta, sample %in% butmeta$sample)
meta <- meta[match(butmeta$sample, meta$sample),]

otu <- otu[otu$group %in% meta$group,]
otu <- otu[match(meta$group, otu$group),]

otu <- otu[,which(names(otu) != 'group')]
but <- but[,which(names(but) != 'group')]


sort(rowSums(otu))
otu <- otu[rowSums(otu)>1000,] # throwing out samples with fewer than 1000 reads
meta <- meta[which(meta$sample %in% rownames(otu)),]
butmeta <- butmeta[which(butmeta$sample %in% meta$sample),]
butmeta$sample == meta$sample

but <- but[which(rownames(but) %in% butmeta$sample),]

rownames(but) == butmeta$sample
rowSums(but)

rowSums(otu)

otu <- otu[,colMeans(otu) > 0]
but <- but[,colMeans(but) > 0]
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

colnames(otu) <- swap[colnames(otu)] # renaming OTUs with some taxonomic info here
colnames(but) <- butswap[colnames(but)]

rownames(but) == rownames(otu)

#### THIS WILL TAKE A WHILE ####

but_16S.all <- ccrepe(but, otu, verbose = FALSE, min.subj = 5)
but.all <- ccrepe(but, verbose = FALSE, min.subj = 5)
otu.all <-  ccrepe(otu, verbose = FALSE, min.subj = 5)

but_16S.all.sigs <- ccrepe_to_geomnet(but_16S.all, pcut=0.05, spearcut = 0.6)
but.all.sigs <- ccrepe_to_geomnet(but.all, pcut = 0.05, spearcut = 0.6)
otu.all.sigs <- ccrepe_to_geomnet(otu.all, pcut = 0.01, spearcut = 0.6)

nodes <- rbind(gather_nodes(otu, '16S'),
               gather_nodes(but, 'but'))


######## THIS IS WHERE ENRICH IS LIMITED ##############

enrich <- enrich[enrich$tissue %in% c('D21_feces', 'D19_feces', 'D15_feces', 'D12_feces'),]   

enrich$node


# this limits the butnet to only features that are diffabund

but_16S.all.sigs <- but_16S.all.sigs[but_16S.all.sigs$from %in% enrich$node & but_16S.all.sigs$to %in% enrich$node,] # one of the nodes in a but-otu connection must be diff abund between the two groups
but.all.sigs <- but.all.sigs[but.all.sigs$from %in% enrich$node & but.all.sigs$to %in% enrich$node,] # but-but connections must both be diff abund
otu.all.sigs <- otu.all.sigs[otu.all.sigs$from %in% enrich$node & otu.all.sigs$to %in% enrich$node,] # otu-otu connections must both be diff abund

all <- rbind(but.all.sigs, but_16S.all.sigs, otu.all.sigs)

all <- fortify(as.edgedf(all), nodes)

# removing small subclusters....

filtered4 <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 4)

p.butnet1 <- ggplot(filtered4, aes(from_id = from_id, to_id = to_id, label=from, color=type)) +
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)),
           aes(color = type, label = from_id),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
           labelgeom = 'text') +
  theme_net()

p.butnet1




net <- as.data.frame(ggplot_build(p.butnet1)$data)

#colnames(dif_ab_16S[,c(13,14,15,16,17)])
#colnames(dif_ab_but[,c(7,8,10,11,12)])

#colnames(enrich)[5] <- 'otu'
colnames(enrich)[11] <- 'from'

enrich.feces <- enrich

# all these little dataframes are for getting the plot labels and features just right
# there is probably a way better way of doing this.

net2 <- merge(enrich.feces, net, by = 'from', all = TRUE)

net2 <- net2[!is.na(net2$x),]
net2$type <- NA
net2$type[grep('Otu[0-9][0-9][0-9]', net2$from)] <- 'but'
net2$type[grep('Otu[0-9][0-9][0-9][0-9][0-9]', net2$from)] <- '16S'

net2$label <- butswap[net2$label]
net2$label[is.na(net2$label)] <- swap[net2$from[is.na(net2$label)] ]

colnames(net2)

net3 <- unique(net2[,c(1,9,14,15)])
net4 <- unique(net2[,c(14,15,16,17)])

net4 <- unique(na.omit(net4[net4 != net4[,c(3,4,1,2)],]))

net3$type <- NA
net3$type[grep('Otu.....', net3$from)] <- '16S'
net3$type[grep('Otu...$', net3$from)] <- 'but'

net3$from <- gsub('(.*)\\(.*\\).*','\\1',net3$from)
net3$from <- gsub('(.*) Otu...', '\\1', net3$from)

net3 <- net3[!(grepl('Roseburia inulinivorans', net3$from) & grepl('Control', net3$Treatment)),]
net3 <- net3[!(grepl('Prevotella_7', net3$from) & grepl('Control', net3$Treatment)),]


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
  geom_label_repel(data=net3, aes(label=from, fill=type, fontface= 'italic'), size=3, alpha=.7, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=c('grey50','grey80')) + theme(panel.background = element_blank(),
                                                axis.title = element_blank(),
                                                axis.text = element_blank(),
                                                axis.ticks = element_blank())
p.butnet.final


############################# Correlation Network ####################

# read in data #

# this has cd3- population in it as well #

allflow <- read.table('./data/FS2_all_flow.txt', header = TRUE, as.is = TRUE, check.names = FALSE)
allflow <- filter(allflow, tissue =='cec')
rownames(allflow) <- allflow$pig_num
allflow <- allflow[,-c(1,34:36)]
allflow <- allflow/rowSums(allflow)
sum((colSums(allflow)/14)*100 >.1) # 21 celltypes with average relative abundance over 0.1%

allflow <- allflow[,(colSums(allflow)/14)*100 >.1] # removes cell types with less than 0.05% abundance (% total live single cells)

flowr <- allflow

colMeans(flowr) *100

flownetsupp <- data.frame(cell_type = names(colMeans(flowr) *100), 
                          Percent_tot_live = colMeans(flowr) *100)

flownetsupp$cell_type <- gsub('/', '', flownetsupp$cell_type)

rownames(flownetsupp) <- NULL

write.table(flownetsupp, file = './output/Flow_net_supp.txt', sep = '\t')


meta <- read.table('./data/V4.metadata.txt', header = TRUE, stringsAsFactors = FALSE)
meta <- meta[meta$day == 21 & meta$tissue == 'cecum',]

# misc #
##################### NEED TO CHANGE THIS SECTION HERE #################
misc <- read.table('./data/miscforcorr.txt', header = TRUE, stringsAsFactors = FALSE, sep = ',', as.is = TRUE, check.names = FALSE)

rownames(misc) <- misc$pig_num
misc$control <- ifelse(misc$pig_num %in% control, 1, 0)
misc$RPS <- ifelse(misc$pig_num %in% control, 0, 1)
misc$`ng IgA/mg dry contents`
misc$treatment <- c(rep('control', 7), c(rep('RPS', 7)))
#IgA fig
p.IgA <- ggplot(misc, aes(x=treatment, y=`ng IgA/mg dry contents`, fill=treatment)) +
  geom_boxplot() + geom_jitter(shape=21, size=2, stroke=1, width = .2) +
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS'))

# qPCR fig
qPCRrel <- read.table('./data/Cassidy_rel_expr.txt', header = TRUE, stringsAsFactors = FALSE, as.is = TRUE, check.names = FALSE)

qPCRrel$treatment <- c(rep('control', 7), rep('RPS', 7))

qPCR.m <- melt(qPCRrel)
qPCR.m$design <- paste(qPCR.m$variable, qPCR.m$treatment)

# wilcoxons for the qPCR genes
qPCR.m %>% group_by(variable) %>% do(W=wilcox.test(value~treatment, data = .data)) %>% summarise(variable, Wilcox=W$p.value)



p.qPCR.2 <- qPCR.m %>% filter(variable %in% c('DEFB1', 'IL-6', 'MUC2')) %>%
  ggplot(aes(x=variable, y=value, fill=treatment)) + geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(jitter.width = .2), shape=21, size=2, stroke=1, show.legend = FALSE) + 
  scale_fill_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('Relative expression (log2 fold change)') + xlab('gene')


qPCR.m <- qPCR.m %>% group_by(design) %>% summarise(mean=mean(value), sd=sd(value), se= sd(value)/sqrt(7)) %>% merge(qPCR.m, by = 'design', all = TRUE)
qPCR.sum <- qPCR.m %>% group_by(design) %>% summarise(mean=mean(value), sd=sd(value), se= sd(value)/sqrt(7)) #%>% merge(qPCR.m, by = 'design', all = TRUE)
qPCR.sum$treatment <- gsub('.* ([A-Za-z]+)', '\\1', qPCR.sum$design)
qPCR.sum$gene <- gsub('(.*) ([A-Za-z]+)', '\\1', qPCR.sum$design)

qPCR.sum$minn <- qPCR.sum$mean- qPCR.sum$se
qPCR.sum$maxx <- qPCR.sum$mean+ qPCR.sum$se

#qpcr.m %>% filter(gene %in% c('DEFB1', 'IL.6', 'MUC2'))


# p.qPCR <- filter(qPCR.sum, gene %in% c('DEFB1', 'IL-6', 'MUC2')) %>%
#   ggplot(aes(x=gene, y=mean, group=design, fill=treatment, color=treatment)) +
#   geom_point(size=10, position = position_dodge(.5), show.legend = FALSE)+
#   # stat_summary(fun.y=mean, geom="point", size=5, position = position_dodge(width = .5)) +
#   geom_errorbar(aes(ymin=minn, ymax=maxx), position = position_dodge(.5), width=.2, size=2, show.legend = FALSE) +
#   scale_color_brewer(palette = "Dark2", labels=c('CON', 'RPS')) + ylab('Relative expression (log2 fold change)') + xlab('gene')
# 
# p.qPCR
# 

##

tax <- extract_mothur_tax('./data/V4.final.taxonomy')
swap <- otu_tax_labels(tax)

# VFAS #
vfa.cec <- read.table('./data/cec_redo_R.txt', header = TRUE, stringsAsFactors = FALSE)
vfa.fec <- read.table('./data/fecalplasma_R.txt', header = TRUE, stringsAsFactors = FALSE)

vfas <- rbind(vfa.cec, vfa.fec)
vfas[,2:16] <- vfas[,2:16]*3  # correct for the dillution in our GC protocol
vfas$BCFA <- vfas$isobutyrate + vfas$isovalerate
vfas$total <- rowSums(vfas[,2:16])
colnames(vfas)[8] <- 'lactate'

vfas$sample <- paste(vfas$tissue, 21, vfas$pig_num, sep = '_')
vfas$treatment <- ifelse(vfas$pig_num %in% c(67,68,69,70,71,72,73,81,82,83,84,85,86,87), 'control', 'RPS')
vfas$design <- paste(vfas$tissue, vfas$treatment, sep = '_')
colnames(vfas)

cecbact <- read.table('./data/cecum.shared', header = TRUE)
rownames(cecbact) <- meta$pig_num
cecbact <- cecbact[,-c(1,2,3)]
rowSums(cecbact)  #31441 seqs per sample
cecbact <- cecbact/rowSums(cecbact)  # converts to relative abundance
cecbact <- cecbact[,colSums(cecbact)>.0001] #removes OTUs with less than 0.01% abundance across all cecal tissue samples
colnames(cecbact) <- swap[colnames(cecbact)]

########## trying to insert but otus ###########
# cut this
# 
# but <- read.table('./data/but_shared_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
# but$group
# butmeta <- read.table('./data/but_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE)
# but$group == butmeta$group
# butmeta$sample <- paste(butmeta$tissue, butmeta$day, butmeta$pig_num, sep = '_')
# rownames(but) <- butmeta$sample
# 
# butmeta <- butmeta %>% filter(butmeta$tissue == 'cecum')
# but <- but %>% filter(rownames(but) %in% butmeta$sample)
# but <- but %>% select(-group)
# rownames(but) <- butmeta$pig_num
# colnames(but) <- butswap[colnames(but)]

# qpcr deltaCTs

qPCR <- read.table('./data/Cec_tissue_forcorr.txt', header = TRUE)
qPCR <- aggregate(x=qPCR, data=qPCR, by=list(qPCR$pig_num), FUN=mean)
rownames(qPCR) <- qPCR$pig_num
qPCR <- qPCR[,-c(1,2)]

#qPCR$ACTb/qPCR$DEFb1
deltaCT <- qPCR - qPCR$ACTb

deltaCT.m <- as.matrix(deltaCT)
deltaCT.m <- 1/deltaCT.m[,-1]



##########

vfas
vfas.m <- as.matrix(filter(vfas, tissue == 'cecum')[,-c(1,7,12,16,17,18,19,20:22)])
#colnames(vfas)
flowrm <- as.matrix(flowr)
cecbactm <- as.matrix(cecbact)
butm <- as.matrix(but)
#deltaCT.m[,8]

#misc2$mucexpr <- deltaCT.m[,8]
misc$dummy <- 1
miscm <- as.matrix(misc)
miscm <- miscm[,c(5,9)]

# Correlation calculations #
set.seed(37)
# THIS IS THE ORIGINAL SWITCH BACK IF WEIRD

TxBut <- ccrepe(x=but, y=flowr, min.subj = 5)
TxB <- ccrepe(x = cecbact, y = flowr, min.subj = 5, verbose = FALSE)
butx16 <- ccrepe(x=cecbact, y=but, min.subj = 5)
#bac <- ccrepe(x = cecbact, min.subj = 7, verbose = FALSE, compare.within.x = TRUE)
im <- ccrepe(x = flowr, min.subj = 5, verbose = FALSE)


vfaVSbut <- rcorr(vfas.m, butm)
vfaVSflow <- rcorr(vfas.m, flowrm)
vfaVSbact <- rcorr(vfas.m, cecbactm)
qPCRvsFlow <- rcorr(deltaCT.m, flowrm)
#qPCRvsBac <- rcorr(deltaCT.m, cecbactm)
miscVSflow <- rcorr(miscm, flowrm)
miscVSqPCR <- rcorr(miscm, deltaCT.m)
#miscVSbac <- rcorr(miscm, cecbactm)
#miscVSbaccont <- rcorr(miscm, cecbactcontm)
vfaVSmisc <- rcorr(vfas.m, miscm)


vfaVSqPCR <- rcorr(vfas.m, deltaCT.m)

#

# Gathering node data #

nodes <- rbind(gather_nodes(flowrm, 'T-cell'),
               gather_nodes(cecbact, '16s'),
               gather_nodes(vfas.m, 'VFA'),
               gather_nodes(deltaCT.m, 'barrier'),
               gather_nodes(misc, 'barrier'))

nodes$type[grep('CD3-', nodes$node)] <- 'CD3neg'
nodes$type[grep('CD3\\+', nodes$node)] <- 'CD3pos'
# Convert to geomnet format #

# muci.sigs <- ccrepe_to_geomnet(TxB, spearcut = -0.6, pcut = 1)
# library(funfuns)
# muci.sigs[grep('Muci', muci.sigs$to),]
# 
# TxBut.sigs <- ccrepe_to_geomnet(TxBut, spearcut = 0.6, pcut = 0.05)
# butx16.sigs <- ccrepe_to_geomnet(butx16)

TxB.sigs <- ccrepe_to_geomnet(TxB, spearcut = 0.6, pcut = 0.01)
im.sigs <- ccrepe_to_geomnet(im, spearcut = 0.6, pcut = 0.05)


#vf.sig <- rcorr_to_geomnet(vfaVSflow, pcut = 0.015)

vf.sig <- rcorr_to_geomnet(vfaVSflow, pcut = 0.05)

qPCR_flow.sig <- rcorr_to_geomnet(qPCRvsFlow, pcut = 0.05, spearcut = 0)
misc_qPCR.sig <- rcorr_to_geomnet(miscVSqPCR, pcut = 0.05)
misc_flow.sig <- rcorr_to_geomnet(miscVSflow, pcut = 0.05)
vfa_misc.sig <- rcorr_to_geomnet(vfaVSmisc, pcut = 0.05)
vfa_qPCR.sig <- rcorr_to_geomnet(vfaVSqPCR, pcut = 0.05)

# vf_but.sig <- rcorr_to_geomnet(vfaVSbut, pcut=0.05)
# 

# these below need to be pruned #

#vfbac.sig <- rcorr_to_ggnet(vfaVSbact, pcut = 0.001)
#qPCR_bact.sig <- rcorr_to_ggnet(qPCRvsBac, pcut = 0.001)
#misc_bac.sig <- rcorr_to_ggnet(miscVSbac, pcut = 0.05)
#bac.sigs <- ccrepe_to_ggnet(bac, pcut = 0.05, spearcut = 0.6)


# removing stuff from bac.sigs and vf.sigs #

# removes tcell-tcell correlations calculated by rcorr()
# vf_but.sig <- vf_but.sig[!(grepl('Otu', vf_but.sig$from) & grepl('Otu', vf_but.sig$to)),]
# vf_but.sig <- vf_but.sig[!(grepl('^[a-z]+$', vf_but.sig$from) & grepl('^[a-z]+', vf_but.sig$to)),]

vf.sig <- vf.sig[!(grepl('CD25', vf.sig$from) & grepl('CD25', vf.sig$to)),]
misc_flow.sig <- misc_flow.sig[!(grepl('CD25', misc_flow.sig$from) & grepl('CD25', misc_flow.sig$to)),]
#misc_bac.sig<- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),]

#vfbac.sig <- vfbac.sig[!(grepl('Otu', vfbac.sig$from) & grepl('Otu', vfbac.sig$to)),]

#qPCR_bact.sig <- qPCR_bact.sig[!(grepl('Otu', qPCR_bact.sig$from) & grepl('Otu', qPCR_bact.sig$to)),]
#misc_bac.sig <- misc_bac.sig[!(grepl('Otu', misc_bac.sig$from) & grepl('Otu', misc_bac.sig$to)),]

# this makes it so that bac to bac correlations are limited to those OTUs which also correlate with other features

#bac.sigs <- rbind(bac.sigs[bac.sigs$from %in% TxB.sigs$to,],
#                  bac.sigs[bac.sigs$to %in% TxB.sigs$to,])

qPCR_flow.sig <- qPCR_flow.sig[!(grepl('CD25', qPCR_flow.sig$from) & grepl('CD25', qPCR_flow.sig$to)),]

#qPCR_bact.sig <- qPCR_bact.sig[!(grepl('Otu', qPCR_bact.sig$from) & grepl('Otu', qPCR_bact.sig$to)),]

all <- rbind(TxB.sigs,
             vf.sig,
             im.sigs,
             qPCR_flow.sig,
             misc_flow.sig,
             vfa_misc.sig,
             vfa_qPCR.sig)

all <- fortify(as.edgedf(all), nodes)

# all$from_id


all$label <- gsub('(.*):Otu[0-9]+', '\\1', all$from_id)

filtered <- prune_graph(fortified.edgelist = all, node.dataframe = nodes, min.vert = 5)
# filtered$type

filtered[grep('CD3\\+', filtered$from),]$type <- 'CD3+'
filtered[grep('CD3-', filtered$from),]$type <- 'CD3-'


# Plotting

library(RColorBrewer)
set.seed(1)
palette()
mypal <- colorRampPalette(brewer.pal(8, "Dark2"))(22)

ggplot() + geom_col(aes(factor(1:30), y = 1:30, fill=factor(1:30))) + scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(30)) + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2))
#to see the colors:
colorRampPalette(brewer.pal(8, "Dark2"))(22)
# you can change 10 to any value
# I was using 22 for a bit







mypal <- mypal[c(18,11,6,13,21)]
#mypal[1] <- 'grey65'
# 16s, barrier, cd3-, cd3+, vfa
as.factor(filtered$type)

filtered

p <- ggplot(filtered, aes(from_id = from_id, to_id = to_id, label=from_id, color=type)) +
  geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)),
           aes(color = type),
           linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
           repel = FALSE, fontsize=3.5, singletons = FALSE,labelcolour="black",
           labelgeom = 'label') + scale_color_manual(values = mypal)+
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

newplot$type[which(is.na(newplot$type))] <- 'barrier'


CD25 <- newplot[grep('CD25\\+', newplot$from),]
FoxP3 <- newplot[grep('FoxP3\\+', newplot$from),]
CD4 <- newplot[grep('CD4\\+', newplot$from),]
CD8 <- newplot[grep('CD8a\\+', newplot$from),]

# throw enrigglies up here #


########### Shiny app stuff #############
# I made a shiny app to help explore this network the stuff that follows here
# is the data to help build this app


# Data for interactive boxplots #

rownames(vfas.m) <- c(67:80)

important_nodes <- unique(c(as.character(newplot$from), as.character(newplot$to)))

cecact_int <- cecbact[,colnames(cecbact) %in% important_nodes] %>% rownames_to_column(var = 'pig') %>%
  gather(key=node, value=value, -pig)

deltaCT_int <- deltaCT.m[,colnames(deltaCT.m) %in% important_nodes] %>%  as.data.frame( )%>% rownames_to_column(var = 'pig') %>%
  gather(key=node, value=value, -pig)

vfas_int <- vfas.m[,colnames(vfas.m) %in% important_nodes] %>% as.data.frame() %>% rownames_to_column(var = 'pig') %>%
  gather(key=node, value=value, -pig)


flowr_int <- flowr[,colnames(flowr) %in% important_nodes] %>% rownames_to_column(var = 'pig') %>%
  gather(key=node, value=value, -pig)

misc_int <- miscm[,colnames(miscm) %in% important_nodes] %>% as.data.frame %>% rownames_to_column(var = 'pig') %>%
  gather(key=node, value=value, -pig)


cecact_int$scale_lab <- 'proportion of total community'
flowr_int$scale_lab <- 'proportion of live single cells'
vfas_int$scale_lab <- 'mM'
deltaCT_int$scale_lab <- '1 / deltaCT'
misc_int$scale_lab <- 'ng IgA/mg dry contents'

misc_int$node <- 'ng IgA/mg dry contents'

int_data <- rbind(cecact_int, flowr_int, vfas_int, deltaCT_int, misc_int)
int_data$treatment <- ifelse(int_data$pig %in% FS2.RPS, 'RPS', 'control')

int_data$value <- as.numeric(gsub(' ','',int_data$value))

# colnames(int_data[focors,])

highlight_these <-int_data %>% group_by(node) %>%
  do(wilco=wilcox.test(value~treatment, data = ., paired=FALSE)) %>% 
  summarise(node, Wilcox=wilco$p.value) %>% filter(Wilcox < 0.1)

enrigglies <- int_data %>% group_by(node, treatment) %>% summarise(nodetreatmean=mean(value)) %>% 
  group_by(node) %>% mutate(enriched_in=ifelse(nodetreatmean[1]>nodetreatmean[2], treatment[1], treatment[2])) %>% 
  select(node, enriched_in) %>% unique()

enrigglies <- merge(highlight_these, enrigglies, by='node')


library(RColorBrewer)
enrigglies$color <- ifelse(enrigglies$enriched_in == 'control', "#1B9E77", "#D95F02")

# should be all the data I need for this shiny thing?

write.table(newplot, file = 'shiny_net.txt', quote = TRUE, sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(int_data, file = 'int_data.txt', quote = TRUE, sep = '\t', col.names = TRUE, row.names = FALSE)
write.table(enrigglies, file = 'enrich.txt', quote = TRUE, sep = '\t', col.names = TRUE, row.names = FALSE)
read.table('shiny_net.txt', header = TRUE, sep = '\t')
read.table('int_data.txt', header = TRUE,  sep = '\t')
read.table('enrich.txt', header = TRUE,  sep = '\t')
################# back to the network plot ###########





colnames(enrigglies)[1] <- 'from'
enrigglies <- merge(newplot, enrigglies, by='from') %>% select(from, color, x, y, enriched_in) %>% unique()
# write.table(enrigglies, file = 'enrich.txt', quote = TRUE, sep = '\t', col.names = TRUE, row.names = FALSE)

newplot$label <- gsub('CD3./(CD4./CD8a.)/FoxP3./CD25.', '\\1', newplot$label)
newplot$label <- gsub('(.*):Otu[0-9]+', '\\1', newplot$label)
labbes <- unique(newplot[,c(1,2,4,5,8,30)])

as.factor(newplot$type)

dif_ab_16S


# rev : brown, yellow, green, pink, blue
# norm: blue, pink, green, yellow, brown
# 16s, barrier, cd3-, cd3+, vfa
# 4, 2, 1, 3, 5
# c(4,2,1,3,5)
#c('grey57', )
######## legend needs work ###########

# have to do some wrangling to get the italics just right...

newplot$lab22 <- ifelse(grepl('\\([0-9]+\\)',newplot$label), paste("italic('", newplot$label, "')", sep = ''), paste("plain('", newplot$label, "')", sep = ''))

ital_wrang <- data.frame(l1 = sub('(.*)(\\(.*\\))','\\1',newplot[newplot$type == '16S',]$label), 
           l2 = sub('(.*)(\\(.*\\))','\\2',newplot[newplot$type == '16S',]$label))

ital_wrang$l3 <- paste("italic('", ital_wrang$l1, "')", sep = '')
ital_wrang$l4 <- paste("plain('", ital_wrang$l2, "')", sep = '')
ital_wrang$l5 <- paste(ital_wrang$l3, ital_wrang$l4, sep = '~')

newplot[newplot$type == '16S',]$lab22 <- ital_wrang$l5

labbes <- unique(newplot[,c(1,2,4,5,8,30, 31)])

pl <- ggplot(newplot) +
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=46, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=44, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=42, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=40, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=38, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=36, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=34, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=32, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=30, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=28, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=26, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=24, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  geom_point(data=enrigglies, aes(x=x, y=y, color=color), size=22, alpha=.1, show.legend = FALSE, inherit.aes = FALSE)+
  
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color='black') +
  geom_point(data=CD25, aes(x=x, y=y), color='red', size=12, show.legend = FALSE)+
  geom_point(aes(x=x, y=y, fill=colour), shape=21, size=8, show.legend = FALSE, inherit.aes = FALSE) +
  geom_rect(data=CD8, aes(xmin=x-.0013,ymin=y-.03, xmax=x+.0013, ymax=y+.03), fill='black', color='black', show.legend = FALSE)+
  geom_rect(data=CD4, aes(xmin=x-.03,ymin=y-.0013, xmax=x+.03, ymax=y+.0013), fill='gold', color='gold', show.legend = FALSE)+
  geom_point(data=FoxP3, aes(x=x, y=y), color='blue', size=3.5)+
  geom_label_repel(data=labbes,aes(x=x, y=y, label=lab22, fill=colour),point.padding = unit(.75, 'lines'), fontface="bold", color='black',
                   alpha=.75, show.legend=FALSE, size=3, parse = TRUE) + 
  scale_color_identity()+
  scale_fill_identity()+
  
  
  # this is all for the legend
  annotate('point', x=1.027, y=1.08, color="#7B6B4D", size=8) + annotate('text', x=1.072, y=1.08 , label='VFA', size=3, color='black') +
  annotate('point', x=1.027, y=1.04, color="#BB8714", size=8) + annotate('text', x=1.072, y=1.04 , label='16S', size=3, color='black') +
  annotate('point', x=1.027, y=1.00, color="#BC5266", size=8) + annotate('text', x=1.072, y=1.00, label='barrier', size=3, color='black') +
  annotate('point', x=1.027, y=.96, color="#66A61E", size=8) + annotate('text', x=1.072, y=.96 , label='CD3+', size=3, color='black') +
  annotate('point', x=1.027, y=.92, color="#966A78", size=8) + annotate('text', x=1.072, y=.92 , label='CD3-', size=3, color='black') +
  annotate('point', x=1.027, y=.88, color="grey", size=8) + annotate('text', x=1.072, y=.88 , label='CD4+', size=3, color='black') +
  annotate('point', x=1.027, y=.84, color="grey", size=8) + annotate('text', x=1.072, y=.84 , label='CD8+', size=3, color='black') +
  annotate('point', x=1.027, y=.80, color="red", size=11) +
  annotate('point', x=1.027, y=.80, color="grey", size=8) + annotate('text', x=1.072, y=.80 , label='CD25+', size=3, color='black') +
  annotate('point', x=1.027, y=.76, color="grey", size=8) + annotate('text', x=1.072, y=.76 , label='FoxP3+', size=3, color='black') +
  annotate('rect', xmin=1.004, xmax=1.05, ymin=.878, ymax=.882, color="gold", fill = 'gold') +
  annotate('rect', xmin=1.024, xmax=1.030, ymin=.818, ymax=.862, color="black", fill = 'black') +
  annotate('point', x=1.027, y=.76, color="blue", size=4) +
  
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1)

pl

pl + annotate(geom='text', x=.05,y=.01, label='Subnetwork A', size=7) + annotate(geom='text', x=.05,y=1, label='Subnetwork B', size=7)

plotdata <- pl$data

# CD3- #966A78
# CD3+ #66A61E
# VFA  #7B6B4D
#barri #BC5266
# 16S  #BB8714




# OTU_blast <- data.frame(OTU=sub('.*:(Otu.*)', '\\1', unique(cecact_int$node)),descrip=NA) # these are the OTUs I want blast data for
# was going to add BLAST info for these OTUs but ive run out of time...


#

proc.time() - ptm
sessionInfo()
writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

# I think there is a memory issue with CCREPE, even after the calculations are done 
# all 16 G of my laptop's memory is being used.  I even removed some R objects but it is still
# not available.  I can't call sessionInfo() because of this... grumble...
#















############
library(cowplot)


##### fig 1 ######

fig.1 <- ggdraw() +
  draw_plot(p.alladon + xlim(c(0,22)), 0,0,.5,1)+
  draw_plot(p.16s.disp.time + ylab('distance from group mean'), .5,0,.5,1)+
  draw_plot_label(c('A', 'B'), c(0,.5), c(1,1))
fig.1

p.16s.disp.time

p.alladon + xlim(c(0,22))
p.16s.disp.time + ylab('distance from group mean')

##### Fig 2 #####

fig.2 <- ggdraw() +
  draw_plot(p.tissord, 0,0,.3,1)+
  draw_plot(p.but.tiss.ord, .3,0,.3,1)+
  draw_plot(p.vfas.final, .6,0,.4,1)+
  draw_plot_label(c('A', 'B','C'), c(0,.3,.6), c(1,1,1))
fig.2

p.tissord
p.but.tiss.ord


library(cowplot)
### fig 3 ###
fig.3 <- ggdraw() +
  draw_plot(p.deseq.ileum,0,.75,.45,.25)+
  draw_plot(p.deseq.cecum, 0,.24,.45,.51)+
  draw_plot(p.deseq.colon, 0,0,.45,.24)+
  draw_plot(p.deseq.D21.feces, .45,0,.55,1) +
  draw_plot_label(y=c(1,.75, .24,1), x=c(0,0,0,.45), c('A', 'B', 'C', 'D'))

fig.3

##

fig.4 <- ggdraw() +
  draw_plot(p.otu87.fec, 0,0,.3,1)+
  draw_plot(p.otu87.tiss, .3,0,.7,1)+
  draw_plot_label(c('A', 'B'), c(0,.3), c(1,1), size=15)

fig.4

######## butfig

fig.7 <- ggdraw() +
  draw_plot(p.flow.even, 0,0,.35,.5)+
  draw_plot(p.flow.cecsig,.35, 0,.65,.5 )+
  draw_plot(p.floword.all, 0,.5, 1,.5) +
  draw_plot_label(c('A', 'B', 'C'), c(0,0,.3), c(1,.5,.5), size = 15)

fig.7

#fig.1
fig.5 <- ggdraw() +
  draw_plot(p.but.deseq.D21.ileum,0,.75,.45,.25)+
  draw_plot(p.but.deseq.D21.cecum, 0,.39,.45,.36)+
  draw_plot(p.but.deseq.D21.cec_cont_RNA, 0,0,.45,.39)+
  draw_plot(p.but.deseq.D21, .45,0,.55,1) +
  draw_plot_label(y=c(1,.75, .39,1), x=c(0,0,0,.45), c('A', 'B', 'C', 'D'))
fig.5
fig.7



#
fig.9 <- ggdraw()+
  draw_plot(p.qPCR.2, 0,0,.6,1)+
  draw_plot(p.IgA, .6,0,.4,1)+
  draw_plot_label(x=c(0,.6), y=c(1,1), label = c('A', 'B'))
fig.9

#


sessionInfo()
