
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

#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
registerDoParallel(cores=6)

setwd("../Originaldata")

##########################################################################

month <- c("April", "May", "June", "July", "August", "September", "Otober", "November")

mat <- list()
gd <- list()

for (i in 1:8) { 
	z <- i + 3	
	x <- t(as.matrix(read.table(sprintf('prey.id.prey.matrix.agg.%s.txt',z), header=T, sep ='\t', row.names=1)))
	mat[[i]] <- x[rowSums(x) > 0, colSums(x) > 0]
	mat[[i]][which(mat[[i]]>0)] <- 1
	m <- nrow(mat[[i]])
	n <- ncol(mat[[i]])
	upper <- cbind(matrix(0,m,m), mat[[i]])
	colnames(upper) <- c(rownames(mat[[i]]), colnames(mat[[i]]))
	lower <- cbind(t(mat[[i]]), matrix(0,n,n))
	colnames(lower) <- c(rownames(mat[[i]]), colnames(mat[[i]]))
	sq <- as.matrix(rbind(upper, lower))
	gd[[i]] <- graph.adjacency(sq, mode="undirected", diag=FALSE)
}


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

meta[which(meta > 0)] <- 1
m <- nrow(meta)
n <- ncol(meta)
upper <- cbind(matrix(0,m,m), meta)
	colnames(upper) <- c(rownames(meta), colnames(meta))
lower <- cbind(t(meta), matrix(0,n,n))
	colnames(lower) <- c(rownames(meta), colnames(meta))
sq_meta <- as.matrix(rbind(upper, lower))
gd_meta <- graph.adjacency(sq_meta, mode="undirected", diag=FALSE)
write.table(data.frame(rownames(sq_meta), sq_meta), file="../Output/sq_meta.txt", quote=F, sep='\t', row.names=F)

##########################################################################
# centrality: each month

bet_e <- foreach (i=1:8) %dopar% igraph::edge_betweenness(gd[[i]], directed=FALSE)
bet_v <- foreach (i=1:8) %dopar% igraph::betweenness(gd[[i]], directed=FALSE, normalized=TRUE)
eig <- foreach (i=1:8) %dopar% eigen_centrality(gd[[i]], directed=FALSE, scale=TRUE)
deg <- foreach (i=1:8) %dopar% igraph::degree(gd[[i]], normalized=TRUE)
clo <- foreach (i=1:8) %dopar% igraph::closeness(gd[[i]], normalized=TRUE)

names(bet_e) <- month
names(bet_v) <- month
names(eig) <- month
names(deg) <- month
names(clo) <- month

saveRDS(bet_e, file="../Output/betweenness_e_each.rds")
saveRDS(bet_v, file="../Output/betweenness_v_each.rds")
saveRDS(eig, file="../Output/eigenvectorcent_each.rds")
saveRDS(deg, file="../Output/degree_each.rds")
saveRDS(clo, file="../Output/closeness_each.rds")

##########################################################################
# centrality: metaweb

bet_e_meta <- igraph::edge_betweenness(gd_meta, directed=FALSE)
bet_v_meta <- igraph::betweenness(gd_meta, directed=FALSE, normalized=TRUE)
eig_meta <- eigen_centrality(gd_meta, directed=FALSE, scale=TRUE)
deg_meta <- igraph::degree(gd_meta, normalized=TRUE)
clo_meta <- igraph::closeness(gd_meta, normalized=TRUE)

saveRDS(bet_e_meta, file="../Output/betweenness_e_meta.rds")
saveRDS(bet_v_meta, file="../Output/betweenness_v_meta.rds")
saveRDS(eig_meta, file="../Output/eigenvectorcent_meta.rds")
saveRDS(deg_meta, file="../Output/degree_meta.rds")
saveRDS(clo_meta, file="../Output/closeness_meta.rds")

