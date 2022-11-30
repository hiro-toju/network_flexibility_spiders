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
# Beta-diversity: comparison with consecutive months

beta_cons <- data.frame(foreach(i=1:7, .combine="rbind") %dopar% betalinkr(webs2array(list(mat[[i]], mat[[i+1]])), partitioning="commondenom", binary=FALSE, partition.st=TRUE, partition.rr=FALSE))

Transition <- c("4->5", "5->6", "6->7", "7->8", "8->9", "9->10", "10->11")
rownames(beta_cons) <- Transition

saveRDS(beta_cons, file="../Output/beta_cons_Fq.rds")
write.table(cbind(Transition, beta_cons), file="../Output/beta_cons_Fq.txt", quote=F, sep='\t', row.names=F)


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
# Targeted (selected) romoval of a species
# upper trophic level

each_sel <- list()

for (j in 1:8) {
	each_sel[[j]] <- foreach(i = 1:length(list_u)) %dopar%  as.matrix(mat[[j]][,colnames(mat[[j]])!=list_u[i]])
}

meta_sel <- foreach(i = 1:length(list_u)) %dopar%  as.matrix(meta[,colnames(meta)!=list_u[i]])


beta_removal <- list()

for (j in 1:8) {
	beta_removal[[j]] <- data.frame(foreach(i=1:length(list_u), .combine="rbind") %dopar% betalinkr(webs2array(list(each_sel[[j]][[i]], meta_sel[[i]])), partitioning="commondenom", binary=FALSE, partition.st=TRUE, partition.rr=FALSE))
	rownames(beta_removal[[j]]) <- list_u
}

for (j in 1:8) {
	beta_removal[[j]] <- data.frame(beta_removal[[j]], OSpWN=beta_removal[[j]]$OS/beta_removal[[j]]$WN)
}

saveRDS(beta_removal, file="../Output/beta_removal_u_Fq.rds")

options(scipen=10)

delta <- list()

for (k in 1:ncol(beta_meta)) {
	delta[[k]] <- matrix(NA, nrow=length(list_u), ncol=8)
for (i in 1:length(list_u)) {
for (j in 1:8) {
delta[[k]][i,j] <- beta_meta[j,k] - beta_removal[[j]][i,k]
rownames(delta[[k]]) <- list_u
colnames(delta[[k]]) <- month
}}}

names(delta) <- colnames(beta_meta)

for (i in 1:ncol(beta_meta)) {
	saveRDS(delta[[i]], file=sprintf("../Output/delta_u_%s_Fq.rds", colnames(beta_meta)[i]))
}



##########################################################################
# Targeted (selected) romoval of a species
# lower trophic level


each_sel <- list()

for (j in 1:8) {
	each_sel[[j]] <- foreach(i = 1:length(list_l)) %dopar%  as.matrix(mat[[j]][rownames(mat[[j]])!=list_l[i],])
}

meta_sel <- foreach(i = 1:length(list_l)) %dopar%  as.matrix(meta[rownames(meta)!=list_l[i],])


beta_removal <- list()

for (j in 1:8) {
	beta_removal[[j]] <- data.frame(foreach(i=1:length(list_l), .combine="rbind") %dopar% betalinkr(webs2array(list(each_sel[[j]][[i]], meta_sel[[i]])), partitioning="commondenom", binary=FALSE, partition.st=TRUE, partition.rr=FALSE))
	rownames(beta_removal[[j]]) <- list_l
}

for (j in 1:8) {
	beta_removal[[j]] <- data.frame(beta_removal[[j]], OSpWN=beta_removal[[j]]$OS/beta_removal[[j]]$WN)
}

saveRDS(beta_removal, file="../Output/beta_removal_l_Fq.rds")

options(scipen=10)

delta <- list()

for (k in 1:ncol(beta_meta)) {
	delta[[k]] <- matrix(NA, nrow=length(list_l), ncol=8)
for (i in 1:length(list_l)) {
for (j in 1:8) {
delta[[k]][i,j] <- beta_meta[j,k] - beta_removal[[j]][i,k]
rownames(delta[[k]]) <- list_l
colnames(delta[[k]]) <- month
}}}

names(delta) <- colnames(beta_meta)

for (i in 1:ncol(beta_meta)) {
	saveRDS(delta[[i]], file=sprintf("../Output/delta_l_%s_Fq.rds", colnames(beta_meta)[i]))
}


##########################################################################


presence_u <- matrix(NA, nrow=length(list_u), ncol=8)
presence_l <- matrix(NA, nrow=length(list_l), ncol=8)

for (i in 1:8) {
	presence_u[,i] <- list_u %in% colnames(mat[[i]])
	presence_l[,i] <- list_l %in% rownames(mat[[i]])
}

saveRDS(presence_u, file="../Output/presence_u_Fq.rds")
saveRDS(presence_l, file="../Output/presence_l_Fq.rds")

