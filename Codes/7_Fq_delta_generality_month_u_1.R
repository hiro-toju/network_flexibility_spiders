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
# delta OS vs. interaction generality: each month


dOS_u <- data.frame(Name=rownames(delta_OS_u), delta_OS_u)

dOS_dp <- data.frame()

for (i in month) {
	x <- full_join(dOS_u[, c("Name", i)], dprime_u_each[, c("Name", i)], by="Name")
	y <- data.frame(x, rep(i, times=nrow(x)))
	colnames(y)[2:4] <- c("deltaOS", "dprime", "Month")
	dOS_dp <- rbind(dOS_dp, y)
}

pdf('../Output/generality_dOS_eachmonth_u_Fq.pdf', w= 15, h= 7)

g <- ggplot(dOS_dp, aes(x=1-dprime, y= deltaOS)) + geom_point(colour="dodgerblue3", size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name))

plot(g)
dev.off()

##########################################################################
# delta ST vs. interaction generality: each month

dST_u <- data.frame(Name=rownames(delta_ST_u), delta_ST_u)

dST_dp <- data.frame()

for (i in month) {
	x <- full_join(dST_u[, c("Name", i)], dprime_u_each[, c("Name", i)], by="Name")
	y <- data.frame(x, rep(i, times=nrow(x)))
	colnames(y)[2:4] <- c("deltaST", "dprime", "Month")
	dST_dp <- rbind(dST_dp, y)
}

pdf('../Output/generality_dST_eachmonth_u_Fq.pdf', w= 15, h= 7)

g <- ggplot(dST_dp, aes(x=1-dprime, y= deltaST)) + geom_point(colour="chartreuse4", size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name))

plot(g)
dev.off()

##########################################################################
# delta OSpWN vs. interaction generality: each month

dOSpWN_u <- data.frame(Name=rownames(delta_OSpWN_u), delta_OSpWN_u)

dOSpWN_dp <- data.frame()

for (i in month) {
	x <- full_join(dOSpWN_u[, c("Name", i)], dprime_u_each[, c("Name", i)], by="Name")
	y <- data.frame(x, rep(i, times=nrow(x)))
	colnames(y)[2:4] <- c("deltaOSpWN", "dprime", "Month")
	dOSpWN_dp <- rbind(dOSpWN_dp, y)
}

pdf('../Output/generality_dOSpWN_eachmonth_u_Fq.pdf', w= 15, h= 7)

g <- ggplot(dOSpWN_dp, aes(x=1-dprime, y= deltaOSpWN)) + geom_point(colour="firebrick2", size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name))

plot(g)
dev.off()

##########################################################################
# delta S vs. interaction generality: each month


dS_u <- data.frame(Name=rownames(delta_S_u), delta_S_u)

dS_dp <- data.frame()

for (i in month) {
	x <- full_join(dS_u[, c("Name", i)], dprime_u_each[, c("Name", i)], by="Name")
	y <- data.frame(x, rep(i, times=nrow(x)))
	colnames(y)[2:4] <- c("deltaS", "dprime", "Month")
	dS_dp <- rbind(dS_dp, y)
}

pdf('../Output/generality_dS_eachmonth_u_Fq.pdf', w= 15, h= 7)

g <- ggplot(dS_dp, aes(x=1-dprime, y= deltaS)) + geom_point(colour="darkorchid4", size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name))

plot(g)
dev.off()

##########################################################################
# delta WN vs. interaction generality: each month

dWN_u <- data.frame(Name=rownames(delta_WN_u), delta_WN_u)

dWN_dp <- data.frame()

for (i in month) {
	x <- full_join(dWN_u[, c("Name", i)], dprime_u_each[, c("Name", i)], by="Name")
	y <- data.frame(x, rep(i, times=nrow(x)))
	colnames(y)[2:4] <- c("deltaWN", "dprime", "Month")
	dWN_dp <- rbind(dWN_dp, y)
}

pdf('../Output/generality_dWN_eachmonth_u_Fq.pdf', w= 15, h= 7)

g <- ggplot(dWN_dp, aes(x=1-dprime, y= deltaWN)) + geom_point(colour="burlywood4", size=3) + facet_wrap(. ~ factor(Month, levels=month), ncol=4) + geom_label_repel(aes(label =Name))

plot(g)
dev.off()

##########################################################################

