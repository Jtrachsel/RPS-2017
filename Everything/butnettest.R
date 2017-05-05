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

system('ls')

################# ALL FUNCTIIONS HERE ###############

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


############### BUTNET 1###############

########## but taxonomy ############

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


################# Read in data ####################

tax <- extract_mothur_tax('V4.final.taxonomy')
tax <- tax[,-c(2,3)]
swap <- otu_tax_labels(tax)

meta <- read.table('16S_meta_forcorr.txt', header = TRUE, stringsAsFactors = FALSE) 

otu <- read.table('16S_shared_forcorr.txt', header = TRUE) %>% filter(group %in% meta$group)

rownames(otu) <- meta$sample

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

########## getting everything in order ############
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
# poorly so their abundance is underestimated.  A sqrt transformation helps correct this by reducing the abundance
# of highly abundant features and increasing the abundance of less abundant features, yet the ranks abundances are unaffected.
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

######## fecals only ###########
but_16S.all <- ccrepe(but, otu, verbose = TRUE, min.subj = 7)
but.all <- ccrepe(but, verbose = TRUE, min.subj = 7)
otu.all <-  ccrepe(otu, verbose = TRUE, min.subj = 7)


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

p.butnet1


enrich$from2 <- swap[enrich$otu]
net <- as.data.frame(ggplot_build(p.butnet1)$data)
enrich$otu
en


colnames(dif_ab_16S[,c(13,14,15,16,17)])
colnames(dif_ab_but[,c(7,8,10,11,12)])

colnames(enrich)[5] <- 'from'
enrich <- enrich[!(enrich$tissue %in% c('ileum', 'cecum', 'cec_cont_RNA', 'colon')),]
unique(enrich)


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

nodes
library(ggrepel)
net3$type <- NA
net3$type[grep('Otu.....', net3$from)] <- '16S'
net3$type[grep('Otu...$', net3$from)] <- 'but'

net3$from <- gsub('(.*)\\(.*\\).*','\\1',net3$from)
net3$from <- gsub('(.*) Otu...', '\\1', net3$from)

ggplot(net2, aes(x=x, y=y)) +  
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

###################################### end butnet ##################################

proc.time() - ptm