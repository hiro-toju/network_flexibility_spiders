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

delta_OS_u <- readRDS("../Output/delta_u_OS_Pq.rds")
delta_OS_l <- readRDS("../Output/delta_l_OS_Pq.rds")

delta_ST_u <- readRDS("../Output/delta_u_ST_Pq.rds")
delta_ST_l <- readRDS("../Output/delta_l_ST_Pq.rds")

delta_OSpWN_u <- readRDS("../Output/delta_u_OSpWN_Pq.rds")
delta_OSpWN_l <- readRDS("../Output/delta_l_OSpWN_Pq.rds")

delta_WN_u <- readRDS("../Output/delta_u_WN_Pq.rds")
delta_WN_l <- readRDS("../Output/delta_l_WN_Pq.rds")

delta_S_u <- readRDS("../Output/delta_u_S_Pq.rds")
delta_S_l <- readRDS("../Output/delta_l_S_Pq.rds")

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
dp_l_meta <- left_join(dprime_l_meta, hinfo[,1:11], by="Name")
il <- data.frame(dp_l_meta, col=rep("grey70", times=nrow(dp_l_meta)))

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
# delta OS vs. interaction generality: metaweb

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_l_max=apply(dOS_l[,-1], 1, max))
dOS_dp <- full_join(dOS_max_l, il, by="Name")

g[[1]] <- ggplot(dOS_dp, aes(x=1-dprime_meta, y= dOS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta ST vs. interaction generality: metaweb

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)
dST_max_l <- data.frame(Name=dST_l$Name, dST_l_max=apply(dST_l[,-1], 1, max))
dST_dp <- full_join(dST_max_l, il, by="Name")

g[[2]] <- ggplot(dST_dp, aes(x=1-dprime_meta, y= dST_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta OSpWN vs. interaction generality: metaweb

dOSpWN_l <- data.frame(Name=rownames(delta_OSpWN_l), delta_OSpWN_l)
dOSpWN_max_l <- data.frame(Name=dOSpWN_l$Name, dOSpWN_l_max=apply(dOSpWN_l[,-1], 1, max))
dOSpWN_dp <- full_join(dOSpWN_max_l, il, by="Name")

g[[3]] <- ggplot(dOSpWN_dp, aes(x=1-dprime_meta, y= dOSpWN_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta S vs. interaction generality: metaweb

dS_l <- data.frame(Name=rownames(delta_S_l), delta_S_l)
dS_max_l <- data.frame(Name=dS_l$Name, dS_l_max=apply(dS_l[,-1], 1, max))
dS_dp <- full_join(dS_max_l, il, by="Name")

g[[4]] <- ggplot(dS_dp, aes(x=1-dprime_meta, y= dS_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta WN vs. interaction generality: metaweb

dWN_l <- data.frame(Name=rownames(delta_WN_l), delta_WN_l)
dWN_max_l <- data.frame(Name=dWN_l$Name, dWN_l_max=apply(dWN_l[,-1], 1, max))
dWN_dp <- full_join(dWN_max_l, il, by="Name")

g[[5]] <- ggplot(dWN_dp, aes(x=1-dprime_meta, y= dWN_l_max, color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################
# delta OS vs. interaction generality: metaweb

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_l_max=apply(dOS_l[,-1], 1, max))

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)
dST_max_l <- data.frame(Name=dST_l$Name, dST_l_max=apply(dST_l[,-1], 1, max))

d <- full_join(dOS_max_l, dST_max_l, by="Name")
d2 <- data.frame(d, dOS_dST=d$dOS_l_max/d$dST_l_max)
d3 <- left_join(d2, il, by="Name")

g[[6]] <- ggplot(d3, aes(x=1-dprime_meta, y=log10(dOS_dST), color=col)) + geom_point(size=3) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

##########################################################################

ge <- grid.arrange(g[[2]], g[[5]], g[[4]], g[[1]], g[[3]], g[[6]], nrow=2, ncol=3)

ggsave(plot=ge, filename="../Output/delta_denerality_meta_l_Pq.pdf", w=18, h=11)
     
##########################################################################
# Corrlation

v <- c("OS", "ST", "OSpWN", "S", "WN", "dOS/dST")
r <- matrix(NA, nrow=6, ncol=2)

r[1,1] <- cor.test(1-dOS_dp$dprime_meta, dOS_dp$dOS_l_max, method="spearman")$estimate
r[2,1] <- cor.test(1-dST_dp$dprime_meta, dST_dp$dST_l_max, method="spearman")$estimate
r[3,1] <- cor.test(1-dOSpWN_dp$dprime_meta, dOSpWN_dp$dOSpWN_l_max, method="spearman")$estimate
r[4,1] <- cor.test(1-dS_dp$dprime_meta, dS_dp$dS_l_max, method="spearman")$estimate
r[5,1] <- cor.test(1-dWN_dp$dprime_meta, dWN_dp$dWN_l_max, method="spearman")$estimate
r[6,1] <- cor.test(1-d3$dprime_meta, d3$dOS_dST, method="spearman")$estimate

r[1,2] <- cor.test(1-dOS_dp$dprime_meta, dOS_dp$dOS_l_max, method="spearman")$p.value
r[2,2] <- cor.test(1-dST_dp$dprime_meta, dST_dp$dST_l_max, method="spearman")$p.value
r[3,2] <- cor.test(1-dOSpWN_dp$dprime_meta, dOSpWN_dp$dOSpWN_l_max, method="spearman")$p.value
r[4,2] <- cor.test(1-dS_dp$dprime_meta, dS_dp$dS_l_max, method="spearman")$p.value
r[5,2] <- cor.test(1-dWN_dp$dprime_meta, dWN_dp$dWN_l_max, method="spearman")$p.value
r[6,2] <- cor.test(1-d3$dprime_meta, d3$dOS_dST, method="spearman")$p.value

FDR <- p.adjust(r[,2], method="fdr")

r2 <- data.frame(Index=v, rho=r[,1], P=r[,2], FDR)

write.table(r2, file="../Output/correlation_delta_generality_l_Pq.txt", sep='\t', quote=F, row.names=F)
   
