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

z <- data.frame(Name=names(degree_meta), degree_meta)

g <- list()

##########################################################################
# delta OS vs. interaction degree: metaweb

dOS_u <- data.frame(Name=rownames(delta_OS_u), delta_OS_u)
dOS_max_u <- data.frame(Name=dOS_u$Name, dOS_u_max=apply(dOS_u[,-1], 1, max))
dOS_dp <- full_join(dOS_max_u, z[-grep('H_', z$Name),], by="Name")

g[[1]] <- ggplot(dOS_dp, aes(x=degree_meta, y= dOS_u_max)) + geom_point(colour="dodgerblue3", size=3) + geom_label_repel(aes(label =Name))

##########################################################################
# delta ST vs. interaction degree: metaweb

dST_u <- data.frame(Name=rownames(delta_ST_u), delta_ST_u)
dST_max_u <- data.frame(Name=dST_u$Name, dST_u_max=apply(dST_u[,-1], 1, max))
dST_dp <- full_join(dST_max_u, z[-grep('H_', z$Name),], by="Name")

g[[2]] <- ggplot(dST_dp, aes(x= degree_meta, y= dST_u_max)) + geom_point(colour="chartreuse4", size=3) + geom_label_repel(aes(label =Name))

##########################################################################
# delta OSpWN vs. interaction degree: metaweb

dOSpWN_u <- data.frame(Name=rownames(delta_OSpWN_u), delta_OSpWN_u)
dOSpWN_max_u <- data.frame(Name=dOSpWN_u$Name, dOSpWN_u_max=apply(dOSpWN_u[,-1], 1, max))
dOSpWN_dp <- full_join(dOSpWN_max_u, z[-grep('H_', z$Name),], by="Name")

g[[3]] <- ggplot(dOSpWN_dp, aes(x= degree_meta, y= dOSpWN_u_max)) + geom_point(colour="firebrick2", size=3) + geom_label_repel(aes(label =Name))

##########################################################################
# delta S vs. interaction degree: metaweb

dS_u <- data.frame(Name=rownames(delta_S_u), delta_S_u)
dS_max_u <- data.frame(Name=dS_u$Name, dS_u_max=apply(dS_u[,-1], 1, max))
dS_dp <- full_join(dS_max_u, z[-grep('H_', z$Name),], by="Name")

g[[4]] <- ggplot(dS_dp, aes(x= degree_meta, y= dS_u_max)) + geom_point(colour="darkorchid4", size=3) + geom_label_repel(aes(label =Name))

##########################################################################
# delta WN vs. interaction degree: metaweb

dWN_u <- data.frame(Name=rownames(delta_WN_u), delta_WN_u)
dWN_max_u <- data.frame(Name=dWN_u$Name, dWN_u_max=apply(dWN_u[,-1], 1, max))
dWN_dp <- full_join(dWN_max_u, z[-grep('H_', z$Name),], by="Name")

g[[5]] <- ggplot(dWN_dp, aes(x= degree_meta, y= dWN_u_max)) + geom_point(colour="burlywood4", size=3) + geom_label_repel(aes(label =Name))

##########################################################################
# delta OS/delta ST vs. interaction degree: metaweb

dOS_u <- data.frame(Name=rownames(delta_OS_u), delta_OS_u)
dOS_max_u <- data.frame(Name=dOS_u$Name, dOS_u_max=apply(dOS_u[,-1], 1, max))

dST_u <- data.frame(Name=rownames(delta_ST_u), delta_ST_u)
dST_max_u <- data.frame(Name=dST_u$Name, dST_u_max=apply(dST_u[,-1], 1, max))
dST_dp <- full_join(dST_max_u, z[-grep('H_', z$Name),], by="Name")

d <- full_join(dOS_max_u, dST_dp, by="Name")
d2 <- data.frame(d, dOS_dST=d$dOS_u_max/d$dST_u_max)
#d3 <- left_join(d2, data.frame(Name=names(degree_meta), degree_meta), by="Name")

g[[6]] <- ggplot(d2, aes(x=degree_meta, y=log10(dOS_dST))) + geom_point(colour="goldenrod2", size=3) + geom_label_repel(aes(label =Name))

##########################################################################

ge <- grid.arrange(g[[2]], g[[5]], g[[4]], g[[1]], g[[3]], g[[6]], nrow=2, ncol=3)

ggsave(plot=ge, filename="../Output/delta_degree_meta_u_Fq.pdf", w=18, h=11)
   
ge2 <- grid.arrange(g[[2]], g[[5]], g[[4]], nrow=3, ncol=1)

ggsave(plot=ge2, filename="../Output/delta_degree_meta_u_Fq_2.pdf", w=5, h=14.5)

     
##########################################################################
# Corrlation

v <- c("OS", "ST", "OSpWN", "S", "WN", "dOS/dST")
r <- matrix(NA, nrow=6, ncol=2)

r[1,1] <- cor.test(dOS_dp$degree_meta, dOS_dp$dOS_u_max, method="kendall")$estimate
r[2,1] <- cor.test(dST_dp$degree_meta, dST_dp$dST_u_max, method="kendall")$estimate
r[3,1] <- cor.test(dOSpWN_dp$degree_meta, dOSpWN_dp$dOSpWN_u_max, method="kendall")$estimate
r[4,1] <- cor.test(dS_dp$degree_meta, dS_dp$dS_u_max, method="kendall")$estimate
r[5,1] <- cor.test(dWN_dp$degree_meta, dWN_dp$dWN_u_max, method="kendall")$estimate
r[6,1] <- cor.test(d2$degree_meta, d2$dOS_dST, method="kendall")$estimate

r[1,2] <- cor.test(dOS_dp$degree_meta, dOS_dp$dOS_u_max, method="kendall")$p.value
r[2,2] <- cor.test(dST_dp$degree_meta, dST_dp$dST_u_max, method="kendall")$p.value
r[3,2] <- cor.test(dOSpWN_dp$degree_meta, dOSpWN_dp$dOSpWN_u_max, method="kendall")$p.value
r[4,2] <- cor.test(dS_dp$degree_meta, dS_dp$dS_u_max, method="kendall")$p.value
r[5,2] <- cor.test(dWN_dp$degree_meta, dWN_dp$dWN_u_max, method="kendall")$p.value
r[6,2] <- cor.test(d2$degree_meta, d2$dOS_dST, method="kendall")$p.value

FDR <- p.adjust(r[,2], method="fdr")

r2 <- data.frame(Index=v, rho=r[,1], P=r[,2], FDR)

write.table(r2, file="../Output/correlation_delta_degree_u_Fq.txt", sep='\t', quote=F, row.names=F)
   
   
