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
library(ggplot2)
library(ggraph)
library(ggnetwork)

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

colnames(meta) <- sp.info$species

##########################################################################

il <- hexa.info[rownames(hexa.info) %in% rownames(meta), 1:10]
il <- data.frame(il, col=rep("#B7B7B7", times=nrow(il)))

il[il[,6]=='Diptera',"col"] <- "#FFD2E0" # light pink
il[il[,6]=='Hemiptera',"col"] <-"#05A820"
il[il[,6]=='Hymenoptera',"col"] <-"#001DD7"
il[il[,6]=='Orthoptera',"col"] <-"#08FAD2"
il[il[,6]=='Collembola',"col"] <-"#ae30ff" # purple
il[il[,6]=='Lepidoptera',"col"] <-"#FFD92F"
il[il[,6]=='Coleoptera',"col"] <-"#AB9B0E"


meta2 <- meta[order(il$order),]
il2 <- il[order(il$order),]

upper <- cbind(matrix(0, nrow=nrow(meta2), ncol=nrow(meta2)), meta2)
lower <- cbind(t(meta2), matrix(0, nrow=ncol(meta2), ncol=ncol(meta2)))
square <- rbind(upper, lower)
label <- c(rownames(meta2), colnames(meta2))
rownames(square) <- label
colnames(square) <- label

gob <- graph.adjacency(square, mode='undirected', diag=FALSE, weighted=TRUE)

sp.list <- data.frame(cbind(matrix("Spider", nrow=ncol(meta2), ncol=10), rep("#FC0053", times=ncol(meta2))))
rownames(sp.list) <- colnames(meta2)
colnames(sp.list) <- colnames(il2)
info2 <- rbind(il2, sp.list)

##########################################################################

dOS_l <- data.frame(Name=rownames(delta_OS_l), delta_OS_l)
dOS_max_l <- data.frame(Name=dOS_l$Name, dOS_max=apply(dOS_l[,-1], 1, max))
dOS_u <- data.frame(Name=rownames(delta_OS_u), delta_OS_u)
dOS_max_u <- data.frame(Name=dOS_u$Name, dOS_max=apply(dOS_u[,-1], 1, max))
dOS <- rbind(dOS_max_l, dOS_max_u)

dST_l <- data.frame(Name=rownames(delta_ST_l), delta_ST_l)
dST_max_l <- data.frame(Name=dST_l$Name, dST_max=apply(dST_l[,-1], 1, max))
dST_u <- data.frame(Name=rownames(delta_ST_u), delta_ST_u)
dST_max_u <- data.frame(Name=dST_u$Name, dST_max=apply(dST_u[,-1], 1, max))
dST <- rbind(dST_max_l, dST_max_u)

dWN_l <- data.frame(Name=rownames(delta_WN_l), delta_WN_l)
dWN_max_l <- data.frame(Name=dWN_l$Name, dWN_max=apply(dWN_l[,-1], 1, max))
dWN_u <- data.frame(Name=rownames(delta_WN_u), delta_WN_u)
dWN_max_u <- data.frame(Name=dWN_u$Name, dWN_max=apply(dWN_u[,-1], 1, max))
dWN <- rbind(dWN_max_l, dWN_max_u)

dS_l <- data.frame(Name=rownames(delta_S_l), delta_S_l)
dS_max_l <- data.frame(Name=dS_l$Name, dS_max=apply(dS_l[,-1], 1, max))
dS_u <- data.frame(Name=rownames(delta_S_u), delta_S_u)
dS_max_u <- data.frame(Name=dS_u$Name, dS_max=apply(dS_u[,-1], 1, max))
dS <- rbind(dS_max_l, dS_max_u)

##########################################################################

shape <- c(rep(16, times=nrow(meta2)), rep(18, times=ncol(meta2)))

info3 <- cbind(info2, shape)
colnames(info3)[1] <- "Name"
info3[,1]  <- rownames(info3)
info4 <- left_join(info3, dOS, by='Name')
info5 <- left_join(info4, dST, by='Name')
info6 <- left_join(info5, dWN, by='Name')
info7 <- left_join(info6, dS, by='Name')

##########################################################################

write.table(info7, file="Metanetwork.Info.txt", quote=F, sep='\t', row.names=F)
write.table(meta2, file="Metanetwork.Topology.txt", quote=F, sep='\t', row.names=F)

##########################################################################
# Node name

rownames(square) <- label
colnames(square) <- label

arrow <- graph.adjacency(square, mode="lower", weighted=T)
weight <- E(arrow)$weight
edge <- get.edgelist(arrow)
el <- data.frame(cbind(edge, weight))

write.table(el, file="Metanetwork.edge.weight_name.txt", col.names=F, row.names=F, sep="\t", quote=F)

##########################################################################

colnames(el) <- c('Prey', 'Spider', 'Count')
colnames(dOS_max_l) <- c('Prey', 'dOS_max_Prey')
colnames(dOS_max_u) <- c('Spider', 'dOS_max_Spider')

el2 <- left_join(el, dOS_max_l, by="Prey")
el3 <- left_join(el2, dOS_max_u, by="Spider")

x <- el3$dOS_max_Prey * el3$dOS_max_Spider
dOS_x_dOS <- (x-min(x)) / (max(x) - min(x))  # standardization 0 - 1

x <- edge_betweenness(gob, directed=FALSE, weights=NULL)
Edge_Betweenness <- (x-min(x)) / (max(x) - min(x))  # standardization 0 - 1

el4 <- data.frame(cbind(el3, dOS_x_dOS), Edge_Betweenness, Label=paste(el4$Prey, el4$Spider, sep='--'))

write.table(el4, file="Edge_dOS_x_dOS.txt", col.names=T, row.names=F, sep="\t", quote=F)

el5 <- subset(el4, as.numeric(el4$Count) > 9)
#plot(el5$Edge_Betweenness, el5$dOS_x_dOS)


g <- ggplot(el5, aes(x=Edge_Betweenness, y=dOS_x_dOS, col=as.numeric(Count))) + geom_point(size=5) + geom_label_repel(aes(label=Label)) + scale_color_gradient(low = "lightpink  ", high = "violetred ")

ggsave(plot=g, filename="dOS_x_dOS_Pq.pdf", w=8, h=7)
ggsave(plot=g, filename="dOS_x_dOS_Pq_label.pdf", w=12, h=10)

##########################################################################
# Node number

rownames(square) <- 1:nrow(square)
colnames(square) <- 1:ncol(square)

arrow <- graph.adjacency(square, mode="lower", weighted=T)
weight <- E(arrow)$weight
edge <- get.edgelist(arrow)
data <- cbind(edge, weight)

write.table(data, file="Metanetwork.edge.weight.txt", col.names=F, row.names=F, sep="\t", quote=F)

##########################################################################
