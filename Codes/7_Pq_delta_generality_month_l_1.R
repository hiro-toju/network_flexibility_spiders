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

##########################################################################
# delta OS vs. interaction generality: each month

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)

dOS_dp <- data.frame()

for (i in month) {
	x <- full_join(dOS_l[, c("Name", i)], dprime_l_each[, c("Name", i)], by="Name")
	x2 <- left_join(x, hinfo[,1:11], by="Name")
	x3 <- data.frame(x2, col=rep("grey70", times=nrow(x2)))
	y <- data.frame(x3, rep(i, times=nrow(x3)))
	colnames(y)[c(2,3,15)] <- c("deltaOS", "dprime", "Month")
	dOS_dp <- rbind(dOS_dp, y)
}

dOS_dp[dOS_dp$order=='Diptera',"col"] <- "#FF99CC"
dOS_dp[dOS_dp$order=='Hemiptera',"col"] <-"#05A820"
dOS_dp[dOS_dp$order=='Hymenoptera',"col"] <-"#001DD7"
dOS_dp[dOS_dp$order=='Orthoptera',"col"] <-"#08FAD2"
dOS_dp[dOS_dp$order=='Collembola',"col"] <-"#FF0000"
dOS_dp[dOS_dp$order=='Lepidoptera',"col"] <-"#FFD92F"
dOS_dp[dOS_dp$order=='Coleoptera',"col"] <-"#AB9B0E"

pdf('../Output/generality_dOS_eachmonth_l_Pq.pdf', w= 15, h= 7)

g <- ggplot(dOS_dp, aes(x=1-dprime, y= deltaOS, color=col)) + geom_point(size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

plot(g)
dev.off()

##########################################################################
# delta ST vs. interaction generality: each month

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)

dST_dp <- data.frame()

for (i in month) {
	x <- full_join(dST_l[, c("Name", i)], dprime_l_each[, c("Name", i)], by="Name")
	x2 <- left_join(x, hinfo[,1:11], by="Name")
	x3 <- data.frame(x2, col=rep("grey70", times=nrow(x2)))
	y <- data.frame(x3, rep(i, times=nrow(x3)))
	colnames(y)[c(2,3,15)] <- c("deltaST", "dprime", "Month")
	dST_dp <- rbind(dST_dp, y)
}

dST_dp[dST_dp$order=='Diptera',"col"] <- "#FF99CC"
dST_dp[dST_dp$order=='Hemiptera',"col"] <-"#05A820"
dST_dp[dST_dp$order=='Hymenoptera',"col"] <-"#001DD7"
dST_dp[dST_dp$order=='Orthoptera',"col"] <-"#08FAD2"
dST_dp[dST_dp$order=='Collembola',"col"] <-"#FF0000"
dST_dp[dST_dp$order=='Lepidoptera',"col"] <-"#FFD92F"
dST_dp[dST_dp$order=='Coleoptera',"col"] <-"#AB9B0E"

pdf('../Output/generality_dST_eachmonth_l_Pq.pdf', w= 15, h= 7)

g <- ggplot(dST_dp, aes(x=1-dprime, y= deltaST, color=col)) + geom_point(size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

plot(g)
dev.off()

##########################################################################
# delta OSpWN vs. interaction generality: each month

dOSpWN_l <- data.frame(Name=rownames(delta_OSpWN_l), delta_OSpWN_l)

dOSpWN_dp <- data.frame()

for (i in month) {
	x <- full_join(dOSpWN_l[, c("Name", i)], dprime_l_each[, c("Name", i)], by="Name")
	x2 <- left_join(x, hinfo[,1:11], by="Name")
	x3 <- data.frame(x2, col=rep("grey70", times=nrow(x2)))
	y <- data.frame(x3, rep(i, times=nrow(x3)))
	colnames(y)[c(2,3,15)] <- c("deltaOSpWN", "dprime", "Month")
	dOSpWN_dp <- rbind(dOSpWN_dp, y)
}

dOSpWN_dp[dOSpWN_dp$order=='Diptera',"col"] <- "#FF99CC"
dOSpWN_dp[dOSpWN_dp$order=='Hemiptera',"col"] <-"#05A820"
dOSpWN_dp[dOSpWN_dp$order=='Hymenoptera',"col"] <-"#001DD7"
dOSpWN_dp[dOSpWN_dp$order=='Orthoptera',"col"] <-"#08FAD2"
dOSpWN_dp[dOSpWN_dp$order=='Collembola',"col"] <-"#FF0000"
dOSpWN_dp[dOSpWN_dp$order=='Lepidoptera',"col"] <-"#FFD92F"
dOSpWN_dp[dOSpWN_dp$order=='Coleoptera',"col"] <-"#AB9B0E"

pdf('../Output/generality_dOSpWN_eachmonth_l_Pq.pdf', w= 15, h= 7)

g <- ggplot(dOSpWN_dp, aes(x=1-dprime, y= deltaOSpWN, color=col)) + geom_point(size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

plot(g)
dev.off()

##########################################################################
# delta S vs. interaction generality: each month

dS_l <- data.frame(Name=rownames(delta_S_l), delta_S_l)

dS_dp <- data.frame()

for (i in month) {
	x <- full_join(dS_l[, c("Name", i)], dprime_l_each[, c("Name", i)], by="Name")
	x2 <- left_join(x, hinfo[,1:11], by="Name")
	x3 <- data.frame(x2, col=rep("grey70", times=nrow(x2)))
	y <- data.frame(x3, rep(i, times=nrow(x3)))
	colnames(y)[c(2,3,15)] <- c("deltaS", "dprime", "Month")
	dS_dp <- rbind(dS_dp, y)
}

dS_dp[dS_dp$order=='Diptera',"col"] <- "#FF99CC"
dS_dp[dS_dp$order=='Hemiptera',"col"] <-"#05A820"
dS_dp[dS_dp$order=='Hymenoptera',"col"] <-"#001DD7"
dS_dp[dS_dp$order=='Orthoptera',"col"] <-"#08FAD2"
dS_dp[dS_dp$order=='Collembola',"col"] <-"#FF0000"
dS_dp[dS_dp$order=='Lepidoptera',"col"] <-"#FFD92F"
dS_dp[dS_dp$order=='Coleoptera',"col"] <-"#AB9B0E"

pdf('../Output/generality_dS_eachmonth_l_Pq.pdf', w= 15, h= 7)

g <- ggplot(dS_dp, aes(x=1-dprime, y= deltaS, color=col)) + geom_point(size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

plot(g)
dev.off()

##########################################################################
# delta WN vs. interaction generality: each month

dWN_l <- data.frame(Name=rownames(delta_WN_l), delta_WN_l)

dWN_dp <- data.frame()

for (i in month) {
	x <- full_join(dWN_l[, c("Name", i)], dprime_l_each[, c("Name", i)], by="Name")
	x2 <- left_join(x, hinfo[,1:11], by="Name")
	x3 <- data.frame(x2, col=rep("grey70", times=nrow(x2)))
	y <- data.frame(x3, rep(i, times=nrow(x3)))
	colnames(y)[c(2,3,15)] <- c("deltaWN", "dprime", "Month")
	dWN_dp <- rbind(dWN_dp, y)
}

dWN_dp[dWN_dp$order=='Diptera',"col"] <- "#FF99CC"
dWN_dp[dWN_dp$order=='Hemiptera',"col"] <-"#05A820"
dWN_dp[dWN_dp$order=='Hymenoptera',"col"] <-"#001DD7"
dWN_dp[dWN_dp$order=='Orthoptera',"col"] <-"#08FAD2"
dWN_dp[dWN_dp$order=='Collembola',"col"] <-"#FF0000"
dWN_dp[dWN_dp$order=='Lepidoptera',"col"] <-"#FFD92F"
dWN_dp[dWN_dp$order=='Coleoptera',"col"] <-"#AB9B0E"

pdf('../Output/generality_dWN_eachmonth_l_Pq.pdf', w= 15, h= 7)

g <- ggplot(dWN_dp, aes(x=1-dprime, y= deltaWN, color=col)) + geom_point(size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name)) + scale_colour_identity(guide="legend") + theme(legend.position="none")

plot(g)
dev.off()

##########################################################################

