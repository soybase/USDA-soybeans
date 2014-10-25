setwd("~/Documents/R Projects/Soybeans/Shiny/Genealogy")

library(phyViz)
library(igraph)
tree <- read.csv("soybeanGenealogy.csv", stringsAsFactors=FALSE)
tree$year[tree$child=="IA 3023"] <- 2003
treeGraph <- processTreeGraph(tree)
# mygraph <- graph.data.frame(treeGraph, directed=F)

varieties <- unique(c(tree$child,tree$parent))
# varieties <- subset(varieties, !grepl(" x ", varieties) & !is.na(varieties))
varieties <- varieties[order(varieties)]

# kevinbacon <- matrix(1, nrow=length(varieties), ncol=length(varieties))
# row.names(kevinbacon) <- varieties
# dimnames(kevinbacon)[[2]] <- varieties
# kevinbacon[lower.tri(kevinbacon)] <- NA
# kevinbacon[diag(kevinbacon)] <- 0
# 
# kevinbacon <- kevinbacon[order(row.names(kevinbacon)), order(dimnames(kevinbacon)[[2]])]
# 
# library(reshape2)
# kevinbacon <- melt(kevinbacon, value.name = "degree", stringsAsFactors=FALSE)
# kevinbacon$Var1 <- as.character(kevinbacon$Var1)
# kevinbacon$Var2 <- as.character(kevinbacon$Var2)
# 
# kevinbacon$degree <- 1000
# library(parallel)
# kevinbacon$degree <- unlist(mclapply(1:nrow(kevinbacon), function(i){
#   try(getDegree(kevinbacon$Var1[i], kevinbacon$Var2[i], mygraph))
# }, mc.cores=8, mc.preschedule=TRUE))
# # 
# # kevinbacon$degree[kevinbacon$degree<0] <- NA
# # 
# # qplot(data=kevinbacon, x=Var1, y=Var2, geom="tile", fill=degree) + theme(axis.text.x=element_text(angle=90, hjust = 1), axis.title=element_blank()) + scale_fill_continuous("Generational\nDistance")
# 
# write.csv(kevinbacon, "kevinbaconTree.csv", row.names=FALSE)
kevinbacon <- read.csv("kevinbaconTree.csv", stringsAsFactors=FALSE)
