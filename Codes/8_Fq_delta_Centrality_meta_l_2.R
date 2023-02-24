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

deg <- data.frame(Name=names(degree_meta), degree_meta)
bet <- data.frame(Name=names(betweenness_v_meta), betweenness_meta=betweenness_v_meta)
eigen <- data.frame(Name=names(eigenvectorcent_meta[[1]]), eigenvectorcent_meta=eigenvectorcent_meta[[1]])
close <- data.frame(Name=names(closeness_meta), closeness_meta)

deg <- deg[grep('H_', deg$Name),]
bet <- bet[grep('H_', bet$Name),]
eigen <- eigen[grep('H_', eigen$Name),]
close <- close[grep('H_', close$Name),]

d <- right_join(hinfo[,1:11], deg, by="Name")
d <- left_join(d, bet, by="Name")
d <- left_join(d, eigen, by="Name")
d <- left_join(d, close, by="Name")

il <- data.frame(d, col=rep("grey70", times=nrow(d)))

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
# delta OS: metaweb

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_l_max=apply(dOS_l[,-1], 1, max))
dOS_cent <- full_join(dOS_max_l, il, by="Name")

g[[1]] <- ggplot(dOS_cent, aes(x= degree_meta, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none") + xlab('Network degree')

g[[2]] <- ggplot(dOS_cent, aes(x= betweenness_meta, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none") + xlab('Betweenness centrality')

g[[3]] <- ggplot(dOS_cent, aes(x= eigenvectorcent_meta, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none") + xlab('Eigenvector centrality')

g[[4]] <- ggplot(dOS_cent, aes(x= closeness_meta, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none") + xlab('Closeness centrality')

##########################################################################

ge2 <- grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], nrow=4, ncol=1)

ggsave(plot=ge2, filename="../Output/delta_centrality_meta_l_Fq.pdf", w=5, h=18)

     
##########################################################################

