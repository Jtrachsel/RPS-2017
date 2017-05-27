setwd('H:/FS2/NitroPhenEnz/feces_D21/')
install.packages('tidyverse')
install.packages('reshape2')
install.packages('lubridate')

library(tidyverse)
library(reshape2)
library(lubridate)

### D0 stuff ###
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


### D21 stuff ###
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

############# A-D-glu ############
plateinfo.ADglu <- read.table('A-D-glu_layout.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plateinfo.ADglu <- plateinfo.ADglu[,-2]
plate.ADglu <- read.table('D0_21_4ADglucorrectPROT.txt', header = TRUE, sep = '\t', check.names = FALSE, as.is = TRUE, stringsAsFactors = FALSE)
plate.ADglu.melt <- melt(data = plate.ADglu, id.vars = 'Time', variable.name = 'well', value.name = 'abs405') %>% filter(well!='Blank')
plate.ADglu.melt$substrate <- 'A-D-glucopyranoside'

ADglu <- merge(plate.ADglu.melt,plateinfo.ADglu, by = 'well')

#D0_21_4ADglucorrectPROT.txt
#A-D-glu_layout.txt

####### ############

master <- rbind(plate2_4.melt, plate3_5.melt, plate6_7.melt, platesulf.melt)


master <- merge(plateinfo, master, by = 'well', all = TRUE)

master <- rbind(D0, master, ADglu)

master$seconds <- period_to_seconds(hms(master$Time))
#master$seconds


#####

control <- c(67,68,69,70,71,72,73,81,82,83,84,85,86,87)

master$treatment <- ifelse(master$pig %in% control, 'control', 'RPS')



master <- master[order(master$seconds),]
master$seconds

master$treatment[master$pig == 'NTC'] <- 'NC'
master$day <- factor(master$day)
master <- master[master$pig != '',]

master$group <- paste(master$substrate, master$day, master$tissue, master$pig)

filter(master, seconds >4000 & seconds<10000 & tissue == 'feces' & pig != 'NTC') %>%
  ggplot(aes(x=seconds, y=abs405, group=group, color=treatment, linetype=day)) + #geom_text(aes(label=pig))+
  geom_path()+ facet_wrap(~substrate, scales = 'free')

### calc rates with lm and check on r2 to see how linear my data really is.  

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
comm_mat$`4-nitrophenyl sulfate`[comm_mat$`4-nitrophenyl sulfate` < 0] <- 0
comm_mat$muc_deg <- rowSums(comm_mat[,c(2,4,5,8)])

write.table(comm_mat, file='glycoside_hydrolase_matrix.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

scale(comm_mat[,2:8])
###############


MUC2dCT <- c(11.5094865,9.9108935,11.2395315,11.4639575,11.571964,9.6729365,11.6880475,8.1974225,9.222601,9.666652,9.9558245,10.776997,9.433513,11.492316)

cec <- filter(comm_mat, tissue == 'cec_cont')
cec$MUC2dCT <- MUC2
cec$mucDvmucR <- cec$muc_deg * cec$MUC2dCT
t.test(cec$mucDv~cec$treatment)
wilcox.test(cec$mucDv~cec$treatment)

##############################


filter(comm_mat, tissue == 'cec_cont') %>%
  ggplot(aes(x=treatment, y=mucDvmucR, fill = treatment))+
  geom_boxplot()

filter(comm_mat, tissue == 'feces') %>%
  ggplot(aes(x=day, y=muc_deg, fill = treatment))+
  geom_boxplot()


#master.results$day <- as.numeric(master.results$day)
filter(master.results, pig != 'NTC' & tissue == 'feces' & r.squared > 0.7) %>%
  ggplot(aes(x=day, y=rate, group=design)) +
  geom_boxplot(aes(fill=treatment)) + 
  scale_fill_brewer(palette = 'Dark2') + facet_wrap(~substrate, scales = 'free') + expand_limits(y=0)


filter(master.results, pig != 'NTC' & tissue != 'feces' & r.squared > 0.7) %>%
  ggplot(aes(x=day, y=rate, group=design)) +
  geom_boxplot(aes(fill=treatment)) + 
  scale_fill_brewer(palette = 'Dark2') + facet_wrap(~substrate, scales = 'free') + expand_limits(y=0)

master.results.filter <- filter(master.results, pig != 'NTC' & tissue == 'feces' & r.squared > 0.7)
master.results.filter2 <- filter(master.results, pig != 'NTC' & tissue == 'cec_cont' & r.squared > 0.7)
#master.results.filter %>% group_by(substrate) %>% pairwise.wilcox.test(x=.$rate, g=.$design, paired = FALSE)

#spread(master.results.filter)
# need to check that NAs arent removing dups

tests <- data.frame()
for (x in levels(factor(master.results.filter$substrate))){
  print(x)
  temp <- filter(master.results.filter, substrate == x)
  temp <- pairwise.wilcox.test(temp$rate, temp$design)$p.value
  temp <- melt(temp)
  tests <- rbind(tests, temp) 
}

tests <- tests[!is.na(tests$value),]
tests <- tests[tests$value < 0.1,]

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




#####

hmm <- master.results %>% filter(substrate == 'N-acetyl-B-D-glucosamine' & tissue == 'cec_cont')
master.results %>% filter(substrate == 'N-acetyl-B-D-glucosamine'& tissue != 'cec_cont' & treatment != 'NTC') %>% ggplot(aes(x=day, y=rate, group=design)) +
  geom_boxplot(aes(fill=treatment)) + geom_jitter()+
  facet_wrap(~substrate, scales = 'free')




