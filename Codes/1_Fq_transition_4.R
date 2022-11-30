##########################################################################
### Calculation of interaction flexibility for each species
### Beta_OS_prime
### Hirokazu Toju
###	2022.10.31
##########################################################################

library(bipartite)
library(foreach)
library(doParallel)
library(igraph)
library(ggplot2)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(gridExtra)

#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
detectCores()
registerDoParallel(cores=6)

setwd("../Originaldata")

##########################################################################
# loading data

month <- c("April", "May", "June", "July", "August", "September", "Otober", "November")

mat <- list()

for (i in 1:8) { 
	z <- i + 3	
	x <- t(as.matrix(read.table(sprintf('prey.id.prey.matrix.agg.%s.txt',z), header=T, sep ='\t', row.names=1)))
	mat[[i]] <- x[rowSums(x) > 0, colSums(x) > 0]
	}

names(mat) <- month

##########################################################################
# Beta-diversity: comparison with consecutive months

beta_cons <- data.frame(foreach(i=1:7, .combine="rbind") %dopar% betalinkr(webs2array(list(mat[[i]], mat[[i+1]])), partitioning="commondenom", binary=FALSE, partition.st=TRUE, partition.rr=FALSE))

Transition <- 1:7
rownames(beta_cons) <- Transition

beta_cons2 <- data.frame(Transition, beta_cons)

saveRDS(beta_cons, file="../Output/beta_cons_Fq.rds")
write.table(beta_cons2, file="../Output/beta_cons_Fq.txt", quote=F, sep='\t', row.names=F)

##########################################################################
# Figure Transitions: assembled with labels

col1 <- c("dodgerblue3", "chartreuse4", "firebrick2", "darkorchid4", "burlywood4")

beta_cons3 <- data.frame(beta_cons2[, 1:5], OSpWN=beta_cons2$OS/beta_cons2$WN)

g <- list()

beta_cons4 <- beta_cons3[,c(1,3,5,6,2,4)]
bc_m <- melt(beta_cons4, id="Transition")

g[[1]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1)

beta_cons5 <- beta_cons4[,1:3]
bc_m <- melt(beta_cons5, id="Transition")

g[[2]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1)

beta_cons6 <- beta_cons4[,c(1,4)]
bc_m <- melt(beta_cons6, id="Transition")

g[[3]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1[3])

ge <- grid.arrange(g[[1]], g[[2]], g[[3]], nrow=1, ncol=3)

ggsave(plot=ge, filename="../Output/transition_betadiversity_assembled_Fq.pdf", w=9, h=2)
     
##########################################################################
# Figure Transitions: assembled without labels

col1 <- c("dodgerblue3", "chartreuse4", "firebrick2", "darkorchid4", "burlywood4")

beta_cons3 <- data.frame(beta_cons2[, 1:5], OSpWN=beta_cons2$OS/beta_cons2$WN)

g <- list()

beta_cons4 <- beta_cons3[,c(1,3,5,6,2,4)]
bc_m <- melt(beta_cons4, id="Transition")

g[[1]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1) + theme(legend.position="none")

beta_cons5 <- beta_cons4[,1:3]
bc_m <- melt(beta_cons5, id="Transition")

g[[2]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1) + theme(legend.position="none")


beta_cons6 <- beta_cons4[,c(1,4)]
bc_m <- melt(beta_cons6, id="Transition")

g[[3]] <- ggplot(bc_m, aes(x=Transition, y=value, color=variable)) + geom_line() + labs(x= "Transition", y= "Beta-diversity") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1[3]) + theme(legend.position="none")

ge <- grid.arrange(g[[1]], g[[2]], g[[3]], nrow=1, ncol=3)

ggsave(plot=ge, filename="../Output/transition_betadiversity_assembled_nolabel_Fq.pdf", w=6, h=2)
     
##########################################################################

  