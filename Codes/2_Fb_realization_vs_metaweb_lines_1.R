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

month <- c("April", "May", "June", "July", "August", "September", "Otober", "November")

mat <- list()

for (i in 1:8) { 
	z <- i + 3	
	x <- t(as.matrix(read.table(sprintf('prey.id.prey.matrix.agg.%s.txt',z), header=T, sep ='\t', row.names=1)))
	mat[[i]] <- x[rowSums(x) > 0, colSums(x) > 0]
	}

names(mat) <- month

##########################################################################
# http://aoki2.si.gunma-u.ac.jp/lecture/mb-arc/arc040/03212.html

merge2 <- function(dfs, ...)
{
　 # res <- merge2(list(Df1, Df2, Df3), by="ID")
　 base <- dfs[1]
　 lapply(dfs[-1], function(i) base <<- merge(base, i, ...)) # [1]
  return(base)
}

##########################################################################
# Making species list
# Making metaweb

v <- list()

for (i in 1:8) {
	v[[i]] <- data.frame(Name=c(rownames(mat[[i]]), colnames(mat[[i]])))
}

vlist <- merge2(v, by="Name", all=T)

list_u <- vlist[-grep("H_", unlist(vlist)),]
list_l <- vlist[grep("H_", unlist(vlist)),]

meta <- matrix(0, nrow=length(list_l), ncol=length(list_u))
rownames(meta) <- list_l
colnames(meta) <- list_u

for (x in 1:8) {
m <- rownames(mat[[x]])
n <- colnames(mat[[x]])
for (i in 1:length(m)) {
	for (j in 1:length(n)) {	
		meta[rownames(meta)==m[i], colnames(meta)==n[j]] <- meta[rownames(meta)==m[i], colnames(meta)==n[j]] + mat[[x]][i,j]
	}
}}

##########################################################################
# Beta-diversity: comparison with metaweb

beta_meta <- data.frame(foreach(i=1:8, .combine="rbind") %dopar% betalinkr(webs2array(list(mat[[i]], meta)), partitioning="commondenom", binary=TRUE, partition.st=TRUE, partition.rr=FALSE))

beta_meta <- data.frame(beta_meta, OSpWN=beta_meta$OS/beta_meta$WN)

rownames(beta_meta) <- month

saveRDS(beta_meta, file="../Output/beta_meta_Fb.rds")
write.table(cbind(month, beta_meta), file="../Output/beta_meta_Fb.txt", quote=F, sep='\t', row.names=F)

##########################################################################
# Figure: realization vs metaweb

col1 <- c("dodgerblue3", "chartreuse4", "firebrick2", "darkorchid4", "burlywood4")

d <- data.frame(Month=4:11, beta_meta)

g <- list()

d1 <- d[,c(1,3,5,9,2,4)]
m1 <- melt(d1, id="Month")

g[[1]] <- ggplot(m1, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1)

d2 <- d1[,1:3]
m2 <- melt(d2, id="Month")

g[[2]] <- ggplot(m2, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1)

d3 <- d1[,c(1,4)]
m3 <- melt(d3, id="Month")

g[[3]] <- ggplot(m3, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1[3])

ge <- grid.arrange(g[[1]], g[[2]], g[[3]], nrow=1, ncol=3)

ggsave(plot=ge, filename="../Output/realization_vs_metaweb_Fb.pdf", w=9, h=2)

##########################################################################
# Figure: realization vs metaweb: without labels

col1 <- c("dodgerblue3", "chartreuse4", "firebrick2", "darkorchid4", "burlywood4")

d <- data.frame(Month=4:11, beta_meta)

g <- list()

d1 <- d[,c(1,3,5,9,2,4)]
m1 <- melt(d1, id="Month")

g[[1]] <- ggplot(m1, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1) + theme(legend.position="none")

d2 <- d1[,1:3]
m2 <- melt(d2, id="Month")

g[[2]] <- ggplot(m2, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1) + theme(legend.position="none")

d3 <- d1[,c(1,4)]
m3 <- melt(d3, id="Month")

g[[3]] <- ggplot(m3, aes(x=Month, y=value, color=variable)) + geom_line() + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=col1[3]) + theme(legend.position="none")

ge <- grid.arrange(g[[1]], g[[2]], g[[3]], nrow=1, ncol=3)

ggsave(plot=ge, filename="../Output/realization_vs_metaweb_nolabel_Fb.pdf", w=6, h=2)
     
##########################################################################
  