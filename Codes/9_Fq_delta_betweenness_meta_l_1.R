##########################################################################
### Calculation of interaction flexibility for each species
### Beta_OS_prime
### Hirokazu Toju
###	2022.10.31
##########################################################################

library(betalink)
library(igraph)
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(bipartite)
library(foreach)
library(doParallel)
library(dplyr)

library(palmerpenguins)
library(ggrepel)
library(tidyverse)


#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
detectCores()
registerDoParallel(cores=6)

setwd("../Originaldata")

##########################################################################
# loading indices

delta_OS_u <- readRDS("../Output/delta_u_OS_Fq.rds")
delta_OS_l <- readRDS("../Output/delta_l_OS_Fq.rds")

delta_ST_u <- readRDS("../Output/delta_u_ST_Fq.rds")
delta_ST_l <- readRDS("../Output/delta_l_ST_Fq.rds")

delta_OSpWN_u <- readRDS("../Output/delta_u_OSpWN_Fq.rds")
delta_OSpWN_l <- readRDS("../Output/delta_l_OSpWN_Fq.rds")

delta_WN_u <- readRDS("../Output/delta_u_WN_Fq.rds")
delta_WN_l <- readRDS("../Output/delta_l_WN_Fq.rds")

delta_S_u <- readRDS("../Output/delta_u_S_Fq.rds")
delta_S_l <- readRDS("../Output/delta_l_S_Fq.rds")

dprime_u_each <- readRDS("../Output/dprime_u_each.rds")
dprime_l_each <- readRDS("../Output/dprime_l_each.rds")

dprime_u_meta <- readRDS("../Output/dprime_u_meta.rds")
dprime_l_meta <- readRDS("../Output/dprime_l_meta.rds")

betweenness_v_each <- readRDS("../Output/betweenness_v_each.rds")
eigenvectorcent_each <- readRDS("../Output/eigenvectorcent_each.rds")
degree_each <- readRDS("../Output/degree_each.rds")
closeness_each <- readRDS("../Output/closeness_each.rds")

betweenness_v_meta <- readRDS("../Output/betweenness_v_meta.rds")
eigenvectorcent_meta <- readRDS("../Output/eigenvectorcent_meta.rds")
degree_meta <- readRDS("../Output/degree_meta.rds")
closeness_meta <- readRDS("../Output/closeness_meta.rds")

month <- c("April", "May", "June", "July", "August", "September", "Otober", "November")

##########################################################################
# loading node information

spider.info <- read.table('Taxonomy_spider_1.txt',header=T,sep ='\t',row.names=1)	
hexa.info <- read.table('03_taxonomy_table_mer.Hexa_ID_notUn.2.txt',header=T,sep ='\t',row.names=1)	

hinfo <- data.frame(Name=rownames(hexa.info), hexa.info)
bet_l_meta <- right_join(data.frame(Name=names(betweenness_v_meta), Betweenness=betweenness_v_meta), hinfo[,1:11], by="Name")
il <- data.frame(bet_l_meta, col=rep("grey70", times=nrow(bet_l_meta)))

il[il$order=='Diptera',"col"] <- "#FF99CC"
il[il$order=='Hemiptera',"col"] <-"#05A820"
il[il$order=='Hymenoptera',"col"] <-"#001DD7"
il[il$order=='Orthoptera',"col"] <-"#08FAD2"
il[il$order=='Collembola',"col"] <-"#FF0000"
il[il$order=='Lepidoptera',"col"] <-"#FFD92F"
il[il$order=='Coleoptera',"col"] <-"#AB9B0E"

##########################################################################

g <- list()

##########################################################################
# delta OS vs. Betweenness: metaweb

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_l_max=apply(dOS_l[,-1], 1, max))
dOS_bet <- left_join(dOS_max_l, il, by="Name")

g[[1]] <- ggplot(dOS_bet, aes(x=Betweenness, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta ST vs. Betweenness: metaweb

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)
dST_max_l <- data.frame(Name=dST_l$Name, dST_l_max=apply(dST_l[,-1], 1, max))
dST_bet <- left_join(dST_max_l, il, by="Name")

g[[2]] <- ggplot(dST_bet, aes(x=Betweenness, y= dST_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta OSpWN vs. Betweenness: metaweb

dOSpWN_l <- data.frame(Name=rownames(delta_OSpWN_l), delta_OSpWN_l)
dOSpWN_max_l <- data.frame(Name=dOSpWN_l$Name, dOSpWN_l_max=apply(dOSpWN_l[,-1], 1, max))
dOSpWN_bet <- left_join(dOSpWN_max_l, il, by="Name")

g[[3]] <- ggplot(dOSpWN_bet, aes(x=Betweenness, y= dOSpWN_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta S vs. Betweenness: metaweb

dS_l <- data.frame(Name=rownames(delta_S_l), delta_S_l)
dS_max_l <- data.frame(Name=dS_l$Name, dS_l_max=apply(dS_l[,-1], 1, max))
dS_bet <- left_join(dS_max_l, il, by="Name")

g[[4]] <- ggplot(dS_bet, aes(x=Betweenness, y= dS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta WN vs. Betweenness: metaweb

dWN_l <- data.frame(Name=rownames(delta_WN_l), delta_WN_l)
dWN_max_l <- data.frame(Name=dWN_l$Name, dWN_l_max=apply(dWN_l[,-1], 1, max))
dWN_bet <- left_join(dWN_max_l, il, by="Name")

g[[5]] <- ggplot(dWN_bet, aes(x=Betweenness, y= dWN_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta OS/delta ST vs. Betweenness: metaweb

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_l_max=apply(dOS_l[,-1], 1, max))

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)
dST_max_l <- data.frame(Name=dST_l$Name, dST_l_max=apply(dST_l[,-1], 1, max))

d <- full_join(dOS_max_l, dST_max_l, by="Name")
d2 <- data.frame(d, dOS_dST=d$dOS_l_max/d$dST_l_max)
d3 <- left_join(d2, il, by="Name")

g[[6]] <- ggplot(d3, aes(x=Betweenness, y=log10(dOS_dST), color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# output

ge <- grid.arrange(g[[2]], g[[5]], g[[4]], g[[1]], g[[3]], g[[6]], nrow=2, ncol=3)

ggsave(plot=ge, filename="../Output/delta_bet_meta_l_Fq.pdf", w=18, h=11)
     
##########################################################################
# Corrlation

v <- c("OS", "ST", "OSpWN", "S", "WN", "dOS/dST")
r <- matrix(NA, nrow=6, ncol=2)

r[1,1] <- cor.test(dOS_bet$Betweenness, dOS_bet$dOS_l_max, method="spearman")$estimate
r[2,1] <- cor.test(dST_bet$Betweenness, dST_bet$dST_l_max, method="spearman")$estimate
r[3,1] <- cor.test(dOSpWN_bet$Betweenness, dOSpWN_bet$dOSpWN_l_max, method="spearman")$estimate
r[4,1] <- cor.test(dS_bet$Betweenness, dS_bet$dS_l_max, method="spearman")$estimate
r[5,1] <- cor.test(dWN_bet$Betweenness, dWN_bet$dWN_l_max, method="spearman")$estimate
r[6,1] <- cor.test(d3$Betweenness, d3$dOS_dST, method="spearman")$estimate

r[1,2] <- cor.test(dOS_bet$Betweenness, dOS_bet$dOS_l_max, method="spearman")$p.value
r[2,2] <- cor.test(dST_bet$Betweenness, dST_bet$dST_l_max, method="spearman")$p.value
r[3,2] <- cor.test(dOSpWN_bet$Betweenness, dOSpWN_bet$dOSpWN_l_max, method="spearman")$p.value
r[4,2] <- cor.test(dS_bet$Betweenness, dS_bet$dS_l_max, method="spearman")$p.value
r[5,2] <- cor.test(dWN_bet$Betweenness, dWN_bet$dWN_l_max, method="spearman")$p.value
r[6,2] <- cor.test(d3$Betweenness, d3$dOS_dST, method="spearman")$p.value

FDR <- p.adjust(r[,2], method="fdr")

r2 <- data.frame(Index=v, rho=r[,1], P=r[,2], FDR)

write.table(r2, file="../Output/correlation_delta_bet_l_Fq.txt", sep='\t', quote=F, row.names=F)

