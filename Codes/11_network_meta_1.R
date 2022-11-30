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

#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
detectCores()
registerDoParallel(cores=6)

setwd("../Originaldata")

##########################################################################
# loading node information

spider.info <- read.table('Taxonomy_spider_1.txt',header=T,sep ='\t',row.names=1)	

hexa.info <- read.table('03_taxonomy_table_mer.Hexa_ID_notUn.2.txt',header=T,sep ='\t',row.names=1)	

##########################################################################
# loading matrices

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

sp.info <- left_join(data.frame(species=colnames(meta)), spider.info, by="species")

#data.frame(colnames(meta), sp.info$species)

colnames(meta) <- sp.info$abbreviation

##########################################################################
# bipartite network: metaweb: prey order

il <- hexa.info[rownames(hexa.info) %in% rownames(meta), 1:10]
il <- data.frame(il, col=rep("grey70", times=nrow(il)))

il[il[,6]=='Diptera',"col"] <- "#FF99CC"
il[il[,6]=='Hemiptera',"col"] <-"#05A820"
il[il[,6]=='Hymenoptera',"col"] <-"#001DD7"
il[il[,6]=='Orthoptera',"col"] <-"#08FAD2"
il[il[,6]=='Collembola',"col"] <-"#FF0000"
il[il[,6]=='Lepidoptera',"col"] <-"#FFD92F"
il[il[,6]=='Coleoptera',"col"] <-"#AB9B0E"

meta2 <- meta[order(il$order),]
il2 <- il[order(il$order),]

pdf('../Output/network_meta.pdf', w= 20, h= 7)

plotweb(t(meta2),
	method="normal",
	text.rot=90,
	labsize=0.4, empty=F,
	col.interaction="grey40",
	bor.col.interaction="grey40",
	col.high=il2$col, 
	col.low="dodgerblue3",
	bor.col.high=NA,
	bor.col.low=NA)

dev.off()

##########################################################################
