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
	sp.info <- left_join(data.frame(species=colnames(mat[[i]])), spider.info, by="species")
	colnames(mat[[i]]) <- sp.info$abbreviation
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
# bipartite network: each network: prey order

for (i in 1:8) {
d <- mat[[i]]
il <- hexa.info[rownames(hexa.info) %in% rownames(d), 1:10]
il <- data.frame(il, col=rep("grey70", times=nrow(il)))

il[il[,6]=='Diptera',"col"] <- "#FF99CC"
il[il[,6]=='Hemiptera',"col"] <-"#05A820"
il[il[,6]=='Hymenoptera',"col"] <-"#001DD7"
il[il[,6]=='Orthoptera',"col"] <-"#08FAD2"
il[il[,6]=='Collembola',"col"] <-"#FF0000"
il[il[,6]=='Lepidoptera',"col"] <-"#FFD92F"
il[il[,6]=='Coleoptera',"col"] <-"#AB9B0E"

d2 <- d[order(il$order),]
il2 <- il[order(il$order),]

pdf(sprintf('../Output/network_%s.pdf', month[i]), w= 12, h= 5)

plotweb(t(d2),
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
}

##########################################################################
