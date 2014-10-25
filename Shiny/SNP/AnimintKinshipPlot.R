library(plyr)
library(reshape2)
library(animint)
library(ggplot2) # need to use fork of ggplot2: install_github("tdhock/ggplot2") for animint compatibility
library(ggdendro)
library(grid)


# Symmetric kinship matrix except for the first column, which contains the variety names
tmp <- read.csv("kinshipMatrix.txt", sep="\t", skip=1, head=F, stringsAsFactors=F)
names(tmp) <- c("Variety", tmp[,1])

# Fix names to correspond with actual usage
rownames <- tmp$Variety
rownames <- gsub("506-13640-2", "S06-13640-2", rownames)
rownames <- gsub("901-G0\\d{1}_RS_", "", rownames)
rownames <- gsub("_NuGEN", "", rownames)
rownames <- gsub("Gasoy17", "Gasoy 17", rownames)
rownames <- gsub("RS_", "", rownames)
rownames <- gsub("Williams82", "Williams 82", rownames)
rownames <- gsub("AK_004", "A.K. (Harrow)", rownames)
rownames <- gsub("NCR", "R", rownames)

# Enforce row and column names
tmp$Variety <- rownames
names(tmp) <- c("Variety", rownames)

# Create numerical matrix for use in clustering algorithms
meanMatrix <- data.matrix(tmp[,-1])
# Ensure dimnames are accurate
dimnames(meanMatrix)[[1]] <- rownames
dimnames(meanMatrix)[[2]] <- rownames

# Convert kinship (similarity) to distance matrix by subtracting relatedness from 2 (maximum relatedness)
dd.col <- rev(as.dendrogram(hclust(as.dist(2-meanMatrix), method="ward.D")))
dd.row <- as.dendrogram(hclust(as.dist(2-meanMatrix), method="ward.D"))
col.ord <- labels(dd.col)
row.ord <- labels(dd.row)
# Get x and y coordinates
ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)


# Melt kinship matrix into long form
matrixLong <- melt(tmp, id.vars = 1)
names(matrixLong) <- c("variety1", "variety2", "value")
matrixLong$value <- as.numeric(matrixLong$value)
write.csv(matrixLong, "LongKinshipMatrix.csv")

# Make long matrix factor names match clustered ordering (so tree matches x and y coordinates)
matrixLong$variety1 <- factor(matrixLong$variety1, levels=col.ord)
matrixLong$variety2 <- factor(matrixLong$variety2, levels=row.ord)

# Multiplier to ensure the tree has reasonable range compared to the heatmap
dendro.multiplier <- 1

# Set segment boundaries for bounding lines to show selected row/column more clearly
matrixLong$Varieties <- paste(matrixLong$variety1, ", ", matrixLong$variety2, " = ", round(matrixLong$value, 2))
matrixLong$xmin <- as.numeric(matrixLong$variety1)-.5
matrixLong$xmax <- as.numeric(matrixLong$variety1)+.5
matrixLong$ymin <- as.numeric(matrixLong$variety2)-.5
matrixLong$ymax <- as.numeric(matrixLong$variety2)+.5

# Static Plot
qplot(x=as.numeric(variety1), y=as.numeric(variety2), geom="tile", fill=value, data=matrixLong) + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) +  
  xlab("Variety 1") + 
  ylab("Variety 2") + 
  scale_x_continuous(expand=c(0, 1), breaks=1:79, labels=col.ord) + 
  scale_y_continuous(expand=c(0, 1), breaks=1:79, labels=row.ord) + 
  scale_fill_gradient("Relatedness",low="#ffffff", high="#374f6b") + 
  ggtitle("Kinship Matrix-based Relatedness") + 
  coord_equal()+ 
  geom_segment(data=segment(ddata_y), aes(x=x, y=y*dendro.multiplier+80, xend=xend, yend=yend*dendro.multiplier+80), inherit.aes=F) + 
  geom_segment(data=segment(ddata_x), aes(x=y*dendro.multiplier+80, y=x, xend=yend*dendro.multiplier+80, yend=xend), inherit.aes=F)

# Animated, interactive plot using Animint
p <- ggplot() + 
  geom_tile(data=matrixLong, aes(x=as.numeric(variety1), y=as.numeric(variety2), fill=value, color=value, clickSelects=Varieties)) + 
  geom_segment(aes(xend=xmin, x=xmin, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(xend=xmax, x=xmax, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(y=ymin, yend=ymin, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(y=ymax, yend=ymax, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong) + 
  theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) +  
  xlab("Variety 1") + 
  ylab("Variety 2") + 
  scale_x_continuous(expand=c(0, 1), breaks=1:79, labels=col.ord) + 
  scale_y_continuous(expand=c(0, 1), breaks=1:79, labels=row.ord) + 
  scale_fill_gradient("Relatedness",low="#ffffff", high="#374f6b") + 
  scale_color_gradient("Relatedness",low="#ffffff", high="#374f6b") + 
  ggtitle("Kinship Matrix Relatedness") + 
  geom_segment(data=segment(ddata_y), aes(x=x, y=y*dendro.multiplier+80, xend=xend, yend=yend*dendro.multiplier+80), inherit.aes=F) + 
  geom_segment(data=segment(ddata_x), aes(x=y*dendro.multiplier+80, y=x, xend=yend*dendro.multiplier+80, yend=xend), inherit.aes=F) + 
  geom_text(data=matrixLong, aes(x=40, y=-3, label=Varieties, showSelected=Varieties), size=14) +
  theme_animint(width=1000, height=1000)

animint2dir(list(kinship=p), out.dir="Plot")
