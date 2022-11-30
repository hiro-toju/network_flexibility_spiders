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

beta_meta <- data.frame(foreach(i=1:8, .combine="rbind") %dopar% betalinkr(webs2array(list(mat[[i]], meta)), partitioning="commondenom", binary=FALSE, partition.st=TRUE, partition.rr=FALSE))

beta_meta <- data.frame(beta_meta, OSpWN=beta_meta$OS/beta_meta$WN)

rownames(beta_meta) <- month

saveRDS(beta_meta, file="../Output/beta_meta_Fq.rds")
write.table(cbind(month, beta_meta), file="../Output/beta_meta_Fq.txt", quote=F, sep='\t', row.names=F)

##########################################################################
# Figure: realization vs metaweb

#col1 <- c("dodgerblue3", "chartreuse4", "firebrick2", "darkorchid4", "burlywood4")

d <- data.frame(Month=4:11, beta_meta)

d1 <- d[,c(1,5,3)]
m1 <- melt(d1, id="Month")

g <- ggplot(m1, aes(x=factor(Month), y=value, fill=variable)) + geom_bar(stat="identity") + labs(x="Month", y= "Beta-diversity_prime") + scale_y_continuous(limits=c(0,1)) + scale_fill_manual(values=c("chartreuse4", "dodgerblue3"))

ggsave(plot=g, filename="../Output/realization_vs_metaweb_bar_Fq.pdf", w=3.5, h=3)

##########################################################################
