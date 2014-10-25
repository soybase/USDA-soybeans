library(stringr)
library(reshape2)
# library(ggdendro)
library(animint)
library(ggplot2)
library(doMC)
registerDoMC()
library(plyr)
library(dplyr)
library(grid)



setwd("~/Documents/R Projects/Soybeans/Shiny/SNP")

# library(sqldf)
# # Create a connection to the sqlite database
# db <- dbConnect(SQLite(), dbname="SNPdb.sqlite")

# Correct column names
col.names <- c("Chromosome", "Position", "id", "Reference", "Alternate", 
               "qual", "filter", "Allele.Freq", "AlleleRSquared", "DosageRSquared",
               "Variety", "Alt_Allele_Freq", "Genotype_Probability", "Gene_State")
# Read in the data from the TSV
vcfTable <- read.table("MeanPlus0p150SDimputed.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)
names(vcfTable) <- col.names


vcfTable$Alt_Allele_Count <- round(vcfTable$Alt_Allele_Freq)

varieties <- as.character(unique(vcfTable$Variety))
seqnames <- unique(vcfTable$Chromosome[grepl("Chr", vcfTable$Chromosome)])


# remove scaffolds for now to make display easier
vcfTable <- filter(vcfTable, grepl("Chr", Chromosome))

# Change names according to Jim's suggestions.
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("901-G04_RS_5601T"), "5601T")
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("901-G06_RS_Bonus"), "Bonus")
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("RS_"), "")
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("AK_004"), "A.K.")
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("Clark_NuGEN"), "Clark")
vcfTable$Variety <- str_replace(vcfTable$Variety, fixed("NCRaleigh"), "Raleigh")

 
snpList <- vcfTable[,c(1, 2, 4, 5, 8, 11, 12, 13, 14, 15)] %>% group_by(Chromosome, Variety)
snpList <- filter(snpList, Alt_Allele_Freq>0)
save(snpList, file="snpList.rda")
# dbWriteTable(conn = db, name="snpList", value=as.data.frame(snpList), row.names=FALSE, overwrite=T)

nsnips <- snpList%>% group_by(Chromosome, Position) %>% summarize(n=length(Position))
nsnips <- nrow(nsnips)
chr.summary <- snpList %>% group_by(Chromosome) %>% summarize(start=min(Position), end=max(Position))
varieties <- unique(snpList$Variety)
varieties <- varieties[order(varieties)]
varieties <- c(varieties[!grepl("^[1-9]", varieties)], varieties[grepl("^[1-9]", varieties)])
seqnames <- unique(snpList$Chromosome)
save(nsnips, chr.summary, varieties, seqnames, file="ShinyStart.rda")

snp.summary <- table(snpList$Variety)
save(snp.summary, file="snpSummary.rda")

snp.density <- group_by(snpList, Chromosome, Variety) %>%
  do(as.data.frame(density(.$Position, n=2048*4, adjust=0.1, from=1, to=max(.$Position), weights=(.$Alt_Allele_Count)/sum(.$Alt_Allele_Count))[1:2]))
save(snp.density, file="SNPDensity.RData")
# dbWriteTable(conn = db, name="snpDensity", value=as.data.frame(snp.density), row.names=FALSE, overwrite=T)


# num.vars <- 10

n <- length(unique(varieties))
snp <- snpList %>% ungroup %>% group_by(Chromosome, Position, Reference) %>% 
  summarise(A=sum(Alt_Allele_Count*(Alternate=="A")), 
            G=sum(Alt_Allele_Count*(Alternate=="G")), 
            C=sum(Alt_Allele_Count*(Alternate=="C")), 
            T=sum(Alt_Allele_Count*(Alternate=="T")), 
            total=sum(Alt_Allele_Count)) %>%
  mutate(A = A + (Reference=="A")*(2*n-total), 
         G = G + (Reference=="G")*(2*n-total),
         C = C + (Reference=="C")*(2*n-total),
         T = T + (Reference=="T")*(2*n-total)) 
snp <- snp[,1:7]


# library(reshape2)
# snp2 <- as.data.frame(snp)
# snp.counts <- melt(snp2, id.vars=1:4, value.name="Count", variable.name="Nucleotide") %>% arrange(Chromosome, Position)

library(tidyr)
snp.counts <- snp %>% gather(Nucleotide, Count, 4:7)
snp.counts <- filter(snp.counts, Count>0)
save(snp.counts, file="SNPCounts.RData")
# dbWriteTable(conn = db, name="SNPCounts", value=as.data.frame(snp.counts), row.names=FALSE, overwrite=T)

rm(chr.summary, snp, snp.counts, snp.density, snpList, vcfTable)
gc()


## Gmax Annotation
segments.full <- read.table(file="./Gmax_275_Wm82.a2.v1.gene_exons.gff3", sep="\t")
names(segments.full) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group")
segments <- segments.full
# keep non-scaffold segments
segments <- filter(segments.full, !grepl("scaffold", seqname))
names(segments)[1] <- "seqnames"
segments$seqnames <- gsub("Gm", "Chr", segments$seqnames)
segments$group <- as.character(segments$group)
segments$ID <- gsub("ID=", "", word(segments$group, sep=";"))
segments$Name <- gsub("Name=", "", word(segments$group, 2, sep=";"))
segments$Parent <- str_extract(as.character(segments$group), "Parent=.*;")

GlymaIDList <- filter(segments, feature=="gene")
GlymaIDList$numid <- gsub("G", "", str_extract(GlymaIDList$ID, "G[[:digit:]]+{6}"))
GlymaIDList$chrnum <- gsub("Chr", "", GlymaIDList$seqnames)
GlymaIDList$searchstr <- gsub("[Gg]lyma", "", gsub("[wW]m82a2v1", "", gsub(".", "", tolower(GlymaIDList$ID), fixed=TRUE)))

save(GlymaIDList, file="./GlymaID.rda")
# dbWriteTable(conn = db, name="GlymaIDList", value=as.data.frame(GlymaIDList), row.names=FALSE, overwrite=T)
dbDisconnect(db)


rm(GlymaIDList, segments.full, segments)
gc()
# 
# ########################## Geneological information ############################
tree2 <- read.csv("soybean_ped.csv", stringsAsFactors=FALSE)
tree <- melt(tree2, id.vars=c(1, 4:7), measure.vars=2:3, variable.name="parent.type", value.name="parent")
tree$year[tree$child=="IA 3023"] <- 2003

# tree <- subset(tree, !child%in%c("Hack", "Century 84", "Logan", "Miami", "Keller", "Wells II", "Beeson x L15", "Weber", "K1028", "OX383", "Clark 63", "C1453", "Beeson 80"))
# tree$parent <- gsub("Amsoy-71", "Amsoy 71", tree$parent)
# tree$child <- gsub("Amsoy-71", "Amsoy 71", tree$child)
# tree$parent <- gsub("No.171", "No. 171", tree$parent)
# tree$child <- gsub("No.171", "No. 171", tree$child)
# tree[tree$child=="L37-1355" & tree$parent.type=="mother", "parent"] <- "PI 81041"
# tree <- subset(tree, child!="Manchu")
# tree[tree$child=="PI 70502","child"] <- "Richland"

save(tree, file="tree.rda")
rm("tree2")

# ############################## Kinship By SNPs ###############################
tmp <- read.csv("kinshipMatrix.txt", sep="\t", skip=1, head=F, stringsAsFactors=F)
names(tmp) <- c("Variety", tmp[,1])

rownames <- tmp$Variety

# Fix variety names
rownames <- gsub("506-13640-2", "S06-13640-2", rownames)
rownames <- gsub("901-G0\\d{1}_RS_", "", rownames)
rownames <- gsub("_NuGEN", "", rownames)
rownames <- gsub("Gasoy17", "Gasoy 17", rownames)
rownames <- gsub("RS_", "", rownames)
rownames <- gsub("Williams82", "Williams 82", rownames)
rownames <- gsub("AK_004", "A.K. (Harrow)", rownames)
rownames <- gsub("NCR", "R", rownames)

tmp$Variety <- rownames
names(tmp) <- c("Variety", rownames)

# Ensure row and column names are appropriate
meanMatrix <- data.matrix(tmp[,-1])
dimnames(meanMatrix)[[1]] <- rownames
dimnames(meanMatrix)[[2]] <- rownames

# Compute dendrogram from the meanMatrix
dd.col <- rev(as.dendrogram(hclust(as.dist(2-meanMatrix), method="ward.D2")))
dd.row <- as.dendrogram(hclust(as.dist(2-meanMatrix), method="ward.D2"))
col.ord <- labels(dd.col)
row.ord <- labels(dd.row)

ddata_x <- ggdendro::dendro_data.dendrogram(dd.row)
ddata_y <- ggdendro::dendro_data.dendrogram(dd.col)

# Melt matrix into a long form plot-able version
matrixLong.snp <- melt(tmp, id.vars = 1)
names(matrixLong.snp) <- c("variety1", "variety2", "value")
matrixLong.snp$value <- 2-as.numeric(matrixLong.snp$value)

# Order factors correctly
matrixLong.snp$variety1 <- factor(matrixLong.snp$variety1, levels=col.ord)
matrixLong.snp$variety2 <- factor(matrixLong.snp$variety2, levels=row.ord)

# Multiplier to make the tree display appropriately proportioned
dendro.multiplier <- 2

# Tile boundaries correctly offset, and set mouseover to display information as well.
matrixLong.snp$Varieties <- paste0(matrixLong.snp$variety1, ", ", matrixLong.snp$variety2)#, " = ", round(matrixLong.snp$value, 2))
matrixLong.snp$xmin <- as.numeric(matrixLong.snp$variety1)-.5
matrixLong.snp$xmax <- as.numeric(matrixLong.snp$variety1)+.5
matrixLong.snp$ymin <- as.numeric(matrixLong.snp$variety2)-.5
matrixLong.snp$ymax <- as.numeric(matrixLong.snp$variety2)+.5


p1 <- ggplot() +
  theme_animint(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), width=750, height=750) + 
  geom_tile(data=matrixLong.snp, aes(x=as.numeric(variety1), y=as.numeric(variety2), fill=value, color=value, clickSelects=Varieties)) + 
  geom_segment(aes(xend=xmin, x=xmin, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong.snp) + 
  geom_segment(aes(xend=xmax, x=xmax, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong.snp) + 
  geom_segment(aes(y=ymin, yend=ymin, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong.snp) + 
  geom_segment(aes(y=ymax, yend=ymax, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong.snp) + 
  xlab("Variety 1") + 
  ylab("Variety 2") + 
  scale_x_continuous(expand=c(0, 1), breaks=1:79, labels=col.ord) + 
  scale_y_continuous(expand=c(0, 1), breaks=1:79, labels=row.ord) + 
  scale_fill_gradient("Distance", high="#ffffff", low="#374f6b") + 
  scale_color_gradient("Distance", high="#ffffff", low="#374f6b") + 
  ggtitle("SNP Distance") + 
  #   coord_equal()+ 
  geom_segment(data=ggdendro::segment(ddata_y), aes(x=x, y=y*dendro.multiplier+80, xend=xend, yend=yend*dendro.multiplier+80), inherit.aes=F) + 
  geom_segment(data=ggdendro::segment(ddata_x), aes(x=y*dendro.multiplier+80, y=x, xend=yend*dendro.multiplier+80, yend=xend), inherit.aes=F) + 
  geom_text(data=matrixLong.snp, aes(x=40, y=-3, label=paste0(Varieties, " = ", round(value, 3)), showSelected=Varieties), size=14)

# animint2dir(list(kinship=p), out.dir="www", open.browser = F)

# ########################## SNPs By Field Trial Data ###########################

fieldtrials <- read.csv("CombinedFieldTrialsData.csv", header=T, stringsAsFactors=FALSE)
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "IA 3023", "IA3023")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Williams 82", "Williams82")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Amsoy 71", "Amsoy")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Beeson 80", "Beeson")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Corsoy 79", "Corsoy")

fieldtrials <- filter(fieldtrials, fieldtrials$Cultivar%in%varieties)

fieldtrialsmatrix <- expand.grid(variety1=unique(fieldtrials$Cultivar), variety2=unique(fieldtrials$Cultivar), stringsAsFactors = F)
fieldtrialsmatrix <- merge(fieldtrialsmatrix, unique(fieldtrials[,c("Cultivar", "AvgYield")]), by.x="variety1", by.y="Cultivar")
names(fieldtrialsmatrix)[3] <- "yield1"
fieldtrialsmatrix <- merge(fieldtrialsmatrix, unique(fieldtrials[,c("Cultivar", "AvgYield")]), by.x="variety2", by.y="Cultivar")
names(fieldtrialsmatrix)[4] <- "yield2"
fieldtrialsmatrix$Varieties <- paste0(fieldtrialsmatrix$variety1, ", ", fieldtrialsmatrix$variety2)
fieldtrialsmatrix$yielddiff <- with(fieldtrialsmatrix, yield1-yield2)
fieldtrialsmatrix <- merge(fieldtrialsmatrix, matrixLong.snp[,c("Varieties", "value")], all.x=T, all.y=F)


# ########################### Kinship By Family Tree ############################
load("./tree-large.rda")
tree.large <- tree
load("tree.rda")
load("KevinBaconDistance.rda")
load("AllSoybeanPathDF.rda")
kevinbaconPaths <- subset(kevinbaconPaths, Var1%in%rownames & Var2 %in%rownames)
kevinbaconPaths <- ddply(kevinbaconPaths, .(Var1, Var2), transform, order=1:length(Var1), .parallel=T)
kevinbacon <- subset(kevinbacon, Var1%in%rownames & Var2 %in%rownames)[,c("Var1", "Var2", "degree")]
kevinbaconPaths <- merge(kevinbaconPaths, kevinbacon)


names(kevinbacon) <- c("variety1", "variety2", "degree")
names(kevinbaconPaths) <- c("variety1", "variety2", "vertex", "year", "order", "degree")

# fix missing years
kevinbaconPaths$year[kevinbaconPaths$vertex=="PI 88.788"] <- 1930
kevinbaconPaths$year[kevinbaconPaths$vertex=="Peking"] <- 1954
kevinbaconPaths$year[kevinbaconPaths$vertex=="IA3023"] <- 2003
kevinbaconPaths$year <- as.numeric(kevinbaconPaths$year)


# Create matrix of values for similarity matrix
meanMatrix <- acast(kevinbacon, variety1 ~ variety2, value.var="degree")

# Remove NAs - make them extremely dissimilar. 
# Assuming yearly breeding, varieties we're looking at could be at most 75 years apart
#range(tree$year, na.rm=T)
meanMatrix[is.na(meanMatrix) | meanMatrix<0] <- 75

# Compute dendrogram from the meanMatrix
dd.col2 <- rev(as.dendrogram(hclust(as.dist(meanMatrix), method="ward.D2")))
dd.row2 <- as.dendrogram(hclust(as.dist(meanMatrix), method="ward.D2"))
col.ord2 <- labels(dd.col2)
row.ord2 <- labels(dd.row2)

ddata_x2 <- ggdendro::dendro_data.dendrogram(dd.row2)
ddata_y2 <- ggdendro::dendro_data.dendrogram(dd.col2)

# Melt matrix into a long form plot-able version
matrixLong <- kevinbacon
names(matrixLong) <- c("variety1", "variety2", "value")
matrixLong$value <- as.numeric(matrixLong$value)

# Order factors correctly
matrixLong$variety1 <- factor(matrixLong$variety1, levels=col.ord2)
matrixLong$variety2 <- factor(matrixLong$variety2, levels=row.ord2)

# Multiplier to make the tree display appropriately proportioned
dendro.multiplier2 <- .25

# Set mouseover to display information for both matrices.
matrixLong$Varieties <- paste0(matrixLong$variety1, ", ", matrixLong$variety2) #, " = ", round(matrixLong$value, 2))
kevinbaconPaths$Varieties <- paste0(kevinbaconPaths$variety1, ", ", kevinbaconPaths$variety2) #, " = ", round(kevinbaconPaths$degree, 2))

# Dataset of segments for paths
segmentPaths <- ddply(kevinbaconPaths, .(Varieties), function(df){
  if(nrow(df)==1){
    data.frame()
  } else {
    with(df, data.frame(Varieties = unique(Varieties), x=year[-length(Varieties)], y=order[-length(Varieties)],
               xend=year[-1], yend=order[-1], stringsAsFactors=FALSE))
  }
})

# Tile boundaries correctly offset
matrixLong$xmin <- as.numeric(matrixLong$variety1)-.5
matrixLong$xmax <- as.numeric(matrixLong$variety1)+.5
matrixLong$ymin <- as.numeric(matrixLong$variety2)-.5
matrixLong$ymax <- as.numeric(matrixLong$variety2)+.5


p2 <- ggplot() +
  theme_animint(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), width=750, height=750) + 
  geom_tile(data=matrixLong, aes(x=as.numeric(variety1), y=as.numeric(variety2), fill=value, color=value, clickSelects=Varieties)) + 
  geom_segment(aes(xend=xmin, x=xmin, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(xend=xmax, x=xmax, 
                   y=.5, yend=max(as.numeric(variety2))+.5, showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(y=ymin, yend=ymin, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong) + 
  geom_segment(aes(y=ymax, yend=ymax, 
                   xend=.5, x=max(as.numeric(variety1)+.5), showSelected=Varieties), data=matrixLong) + 
  xlab("Variety 1") + 
  ylab("Variety 2") + 
  scale_x_continuous(expand=c(0, 1), breaks=1:length(col.ord2), labels=col.ord2) + 
  scale_y_continuous(expand=c(0, 1), breaks=1:length(row.ord2), labels=row.ord2) + 
  scale_fill_gradient("Generations", high="#ffffff", low="#374f6b") + 
  scale_color_gradient("Generations", high="#ffffff", low="#374f6b") + 
  ggtitle("Generational (Kevin Bacon) Distance") + 
  geom_segment(data=ggdendro::segment(ddata_y2), aes(x=x, y=y*dendro.multiplier2+length(col.ord2)+1, xend=xend, yend=yend*dendro.multiplier2+length(col.ord2)+1), inherit.aes=F) + 
  geom_segment(data=ggdendro::segment(ddata_x2), aes(x=y*dendro.multiplier2+length(row.ord2)+1, y=x, xend=yend*dendro.multiplier2+length(row.ord2)+1, yend=xend), inherit.aes=F) + 
  geom_text(data=matrixLong, aes(x=40, y=-3, label=paste0(Varieties, " = ", value, " generations apart"), showSelected=Varieties), size=14)


p3 <- ggplot() + 
  theme_animint(width=750, height=400) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, showSelected=Varieties), data=segmentPaths) + 
  geom_rect(aes(xmin=year-1.5, xmax=year+1.5, ymin=order-.25, ymax=order+.25, showSelected=Varieties), data=kevinbaconPaths, fill="white", colour="white") + 
  geom_text(aes(x=year, y=order, label=vertex, showSelected=Varieties), data=kevinbaconPaths) + 
  theme_bw() + 
  scale_y_continuous(name="Distance", breaks=c(0, 2, 4, 6, 8, 10, 12), limits=c(-1, 12)) + 
  xlab("Year") + 
  ggtitle("Shortest Path between Varieties")


AllVars <- merge(kevinbacon, matrixLong.snp[,c("variety1", "variety2", "value", "Varieties")])
AllVars$degree.jit <- jitter(AllVars$degree, amount=.4)  
AllVars$order1 <- order(AllVars$variety1)
AllVars$order2 <- order(AllVars$variety2)
AllVars <- subset(AllVars, order1 <= order2)
AllVars <- AllVars[,-c(7:8)]
AllVars2 <- AllVars
names(AllVars2) <- c("variety2", "variety1", names(AllVars2)[-c(1:2)])
AllVars2$Varieties <- with(AllVars2, paste0(variety1, ", ", variety2))
AllVars <- unique(rbind(AllVars, AllVars2))
rm(AllVars2)


p4 <- ggplot() + 
  theme_animint(width=375, height=400) +
  geom_point(aes(y=value, x=degree.jit, showSelected=Varieties), shape=1, colour="black", fill="white", size=3, data=AllVars) + 
  geom_point(aes(y=value, x=degree.jit, clickSelects=Varieties), alpha=1, colour="#374f6b", fill="#374f6b", data=AllVars) +
  scale_x_continuous(name="Generational Distance") + 
  ylab("SNP Distance") + 
  ggtitle("Generations and SNP Distance")

p5 <- ggplot() + 
  theme_animint(width=375, height=400) +
  geom_point(aes(y=value, x=yielddiff, showSelected=Varieties), shape=1, size=3, data=fieldtrialsmatrix) + 
  geom_point(aes(y=value, x=yielddiff, clickSelects=Varieties), alpha=1, colour="#374f6b", fill="#374f6b", data=fieldtrialsmatrix) +
  scale_x_continuous(name="Yield 1 - Yield 2") + 
  ylab("SNP Distance") + 
  ggtitle("Yield Difference and SNP Distance") + 

animint2dir(list(heatmap=p1, kinshipHeatmap=p2, kinship=p3, rel=p4, yield=p5), out.dir="www", open.browser = F)


#save.image(file="ShinyData.RData")
