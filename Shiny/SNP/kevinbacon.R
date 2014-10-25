
library(phyViz)

load("./tree-large.rda")
treeGraph <- processTreeGraph(tree)
# mygraph <- graph.data.frame(treeGraph, directed=F)

varieties <- unique(c(tree$child,tree$parent))
varieties <- subset(varieties, !grepl(" x ", varieties) & !is.na(varieties))
varieties <- varieties[order(varieties)]

kevinbacon <- matrix(1, nrow=length(varieties), ncol=length(varieties))
row.names(kevinbacon) <- varieties
dimnames(kevinbacon)[[2]] <- varieties
kevinbacon[lower.tri(kevinbacon)] <- NA
kevinbacon[diag(kevinbacon)] <- 0

kevinbacon <- kevinbacon[order(row.names(kevinbacon)), order(dimnames(kevinbacon)[[2]])]

library(reshape2)
kevinbacon <- melt(kevinbacon, value.name = "degree", stringsAsFactors=FALSE)
kevinbacon$Var1 <- as.character(kevinbacon$Var1)
kevinbacon$Var2 <- as.character(kevinbacon$Var2)

kevinbacon$degree <- 1000
library(parallel)
library(dplyr)

library(doMC)
registerDoMC(cores=14)

kevinbaconPaths <- do.call("rbind.fill", mclapply(1:nrow(kevinbacon), function(i){
  path <- try(as.data.frame(getPath(kevinbacon$Var1[i], kevinbacon$Var2[i], treeGraph, tree=tree), stringsAsFactors=FALSE))
  if(nrow(path)==0){
    path <- data.frame(pathVertices="No Path Found", yearVertices = NA)
  }
  data.frame(Var1=kevinbacon$Var1[i], Var2=kevinbacon$Var2[i], path, stringsAsFactors=FALSE)
}))


kevinbaconPaths2 <- do(kevinbacon2, function(i){
  path <- as.data.frame(getPath(.$Var1[i], .$Var2[i], treeGraph, tree=tree), stringsAsFactors=FALSE)
  data.frame(Var1=.$Var1[i], Var2=.$Var2[i], path, stringsAsFactors=FALSE)
})

kevinbaconPaths <- rbind(kevinbaconPaths, kevinbaconPaths2)
rm(kevinbaconPaths2)
rm(kevinbacon1, kevinbacon2)


degree <- kevinbaconPaths %>% group_by(Var1, Var2) %>% summarize(nrow)

kevinbacon$degree <- unlist(mclapply(1:nrow(kevinbacon), function(i){
  getDegree(kevinbacon$Var1[i], kevinbacon$Var2[i], treeGraph, tree=tree)
}, mc.cores=12, mc.preschedule=TRUE))

# kevinbacon$degree[kevinbacon$degree==-1] <- 
#   sapply(which(kevinbacon$degree==-1), function(i){
#     getDegree(kevinbacon$Var1[i], kevinbacon$Var2[i], mygraph)
#   })

# kevinbacon <- kevinbacon[order(kevinbacon$Var1, kevinbacon$Var2),]
# 
# kevinbacon$Var1 <- factor(kevinbacon$Var1, levels=varieties[order(varieties)])
# kevinbacon$Var2 <- factor(kevinbacon$Var2, levels=varieties[order(varieties, decreasing=T)])

kevinbacon$degree[kevinbacon$degree<0] <- NA

qplot(data=kevinbacon, x=Var1, y=Var2, geom="tile", fill=degree) + theme(axis.text.x=element_text(angle=90, hjust = 1), axis.title=element_blank()) + scale_fill_continuous("Generational\nDistance")

write.csv(kevinbacon, "kevinbaconTree.csv", row.names=FALSE)
kevinbacon <- read.csv("kevinbaconTree.csv", stringsAsFactors=FALSE)

baconmat <- acast(kevinbacon, Var1 ~ Var2, value.var="degree")
baconmat[is.na(baconmat) | baconmat<0] <- 1000


dd.col <- rev(as.dendrogram(hclust(as.dist(baconmat), method="ward.D2")))
dd.row <- as.dendrogram(hclust(as.dist(baconmat), method="ward.D2"))
col.ord <- labels(dd.col)
row.ord <- labels(dd.row)

library(ggdendro)
library(grid)

ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)

dendro.multiplier <- .25
offset <- length(col.ord)
cutoff <- 250

kevinbacon$Var1fac <- factor(kevinbacon$Var1, levels=col.ord)
kevinbacon$Var2fac <- factor(kevinbacon$Var2, levels=row.ord)

kevinbacon$degree[kevinbacon$degree<0] <- NA

qplot(x=as.numeric(Var1fac), y=as.numeric(Var2fac), geom="tile", fill=degree, data=kevinbacon) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) +  
  xlab("Variety 1") + 
  ylab("Variety 2") + 
  scale_x_continuous(expand=c(0, 1), breaks=1:offset-.5, labels=col.ord, limits=c(0, cutoff)) + 
  scale_y_continuous(expand=c(0, 1), breaks=1:offset-.5, labels=row.ord, limits=c(0, cutoff)) + 
  scale_fill_gradient("Generational\nDistance",low="#ffffff", high="#374f6b") + 
  ggtitle('"Kevin Bacon" Distance between Soybean Varieties') + 
  coord_equal()+ 
  geom_segment(data=segment(ddata_y), aes(x=x, y=y*dendro.multiplier+offset+3, xend=xend, yend=yend*dendro.multiplier+offset+3), inherit.aes=F) + 
  geom_segment(data=segment(ddata_x), aes(x=y*dendro.multiplier+offset+3, y=x, xend=yend*dendro.multiplier+offset+3, yend=xend), inherit.aes=F) 
# ggsave("KevinBaconDistance.png", width=12, height=12)
