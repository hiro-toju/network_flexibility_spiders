
##########################################################################
### Calculation of interaction flexibility for each species
### Beta_OS_prime
### Hirokazu Toju
###	2022.10.31
##########################################################################


library(bipartite)
library(foreach)
library(doParallel)

#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
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
# dprime: each month

r_l <- foreach(i=1:8) %dopar% dfun(mat[[i]])$dprime
r_u <- foreach(i=1:8) %dopar% dfun(t(mat[[i]]))$dprime

dp_l_each <- list()
dp_u_each <- list()

for (i in 1:8) {
	Name <- names(r_l[[i]])
	dp_l_each[[i]]  <- data.frame(Name, r_l[[i]])
	Name <- names(r_u[[i]])
	dp_u_each[[i]]  <- data.frame(Name, r_u[[i]])
}

dprime_l_each <- merge2(dp_l_each, by="Name", all=T)
dprime_u_each <- merge2(dp_u_each, by="Name", all=T)
colnames(dprime_l_each)[2:9] <- month
colnames(dprime_u_each)[2:9] <- month

saveRDS(dprime_l_each, file="../Output/dprime_l_each.rds")
saveRDS(dprime_u_each, file="../Output/dprime_u_each.rds")


##########################################################################
# dprime: metaweb

dprime_meta <- dfun(meta)$dprime
dprime_l_meta <- data.frame(Name=rownames(meta), dprime_meta)
saveRDS(dprime_l_meta, file="../Output/dprime_l_meta.rds")

dprime_meta <- dfun(t(meta))$dprime
dprime_u_meta <- data.frame(Name=colnames(meta), dprime_meta)
saveRDS(dprime_u_meta, file="../Output/dprime_u_meta.rds")


