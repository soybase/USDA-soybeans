# library(BiocInstaller)
# library(cn.mops)


# to generate animint plots:
# --------------------------
# library(devtools)
# install_github("tdhock/ggplot2") # required for theme_animint
# install_github("tdhock/animint") 
library(ggplot2)
library(animint)

library(stringr)
library(reshape2)
library(plyr)

library(ggbio)

library(doMC)
registerDoMC(8)


setwd("~/Documents/R Projects/Soybeans/Shiny/CNV")

# cn.mops data files now compiled on server. All we have to do is load them.
# To get cn.mops run on the server, use 
#   /data006a/GIF_2a/tengfei/data/cnv/scripts/cnMopsFullGenome.R
#   That script will call DataObjectsCnMops.R, 
#   which produces the ShinyDAtaObjects-____________.rda files. 
date <- "100414%H5228"
load(paste("./ShinyDataObjects-", date, ".rda", sep=""))

# res.df$Variety <- gsub("Sample_", "", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("Sample_", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("Sample_", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("Sample_", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub(".", "-", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub(".", "-", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub(".", "-", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub(".", "-", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("RS_", "", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("RS_", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("RS_", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("RS_", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("_", " ", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("_", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("_", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("_", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub(" ", "", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub(" ", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub(" ", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub(" ", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X901-G06 ", "", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X901-G06 ", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X901-G06 ", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X901-G06 ", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("AK004", "A.K. (Harrow)", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("AK004", "A.K. (Harrow)", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("AK004", "A.K. (Harrow)", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("AK004", "A.K. (Harrow)", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("ClarkNuGEN", "Clark", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("ClarkNuGEN", "Clark", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("ClarkNuGEN", "Clark", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("ClarkNuGEN", "Clark", glymaIDs$Variety, fixed=TRUE)
# 
# res.df$Variety <- gsub("S-100", "S 100", res.df$Variety, fixed=TRUE)
# segranges.df$Variety <- gsub("S-100", "S 100", segranges.df$Variety, fixed=TRUE)
# segranges.df2$Variety <- gsub("S-100", "S 100", segranges.df2$Variety, fixed=TRUE)
# glymaIDs$Variety <- gsub("S-100", "S 100", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X901-G04", "", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X901-G04", "", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X901-G04", "", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X901-G04", "", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("IA3023x", "IA 3023 (NAM)", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("IA3023x", "IA 3023 (NAM)", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("IA3023x", "IA 3023 (NAM)", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("IA3023x", "IA 3023 (NAM)", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X180501", "PI 180.501", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X180501", "PI 180.501", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X180501", "PI 180.501", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X180501", "PI 180.501", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X506-13640-2", "S06-130640-2", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X506-13640-2", "S06-130640-2", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X506-13640-2", "S06-130640-2", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X506-13640-2", "S06-130640-2", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X88788", "PI 88.788", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X88788", "PI 88.788", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X88788", "PI 88.788", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X88788", "PI 88.788", glymaIDs$Variety, fixed=TRUE)
# 
# res.df$Variety <- gsub("5601T", "TN 5601T", res.df$Variety, fixed=TRUE)
# segranges.df$Variety <- gsub("5601T", "TN 5601T", segranges.df$Variety, fixed=TRUE)
# segranges.df2$Variety <- gsub("5601T", "TN 5601T", segranges.df2$Variety, fixed=TRUE)
# glymaIDs$Variety <- gsub("5601T", "TN 5601T", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("X901-G06Bonus", "Bonus", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("X901-G06Bonus", "Bonus", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("X901-G06Bonus", "Bonus", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("X901-G06Bonus", "Bonus", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("NCRaleigh", "NC Raleigh", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("NCRaleigh", "NC Raleigh", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("NCRaleigh", "NC Raleigh", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("NCRaleigh", "NC Raleigh", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("NCRoy", "NC Roy", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("NCRoy", "NC Roy", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("NCRoy", "NC Roy", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("NCRoy", "NC Roy", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("Gasoy17", "Gasoy 17", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("Gasoy17", "Gasoy 17", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("Gasoy17", "Gasoy 17", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("Gasoy17", "Gasoy 17", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("Williams82", "Williams 82", res.df$Variety, fixed=TRUE)
segranges.df$Variety <- gsub("Williams82", "Williams 82", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("Williams82", "Williams 82", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("Williams82", "Williams 82", glymaIDs$Variety, fixed=TRUE)

# res.df$Variety <- gsub("(PI){1}( {0,})(\\d{2,3}?)(\\.{0,})(\\d{3})(-{0,}\\w{0,})", "\\1 \\3.\\5\\6", res.df$Variety)
segranges.df$Variety <- gsub("(PI){1}( {0,})(\\d{2,3}?)(\\.{0,})(\\d{3})(-{0,}\\w{0,})", "\\1 \\3.\\5\\6", segranges.df$Variety, fixed=TRUE)
segranges.df2$Variety <- gsub("(PI){1}( {0,})(\\d{2,3}?)(\\.{0,})(\\d{3})(-{0,}\\w{0,})", "\\1 \\3.\\5\\6", segranges.df2$Variety, fixed=TRUE)
glymaIDs$Variety <- gsub("(PI){1}( {0,})(\\d{2,3}?)(\\.{0,})(\\d{3})(-{0,}\\w{0,})", "\\1 \\3.\\5\\6", glymaIDs$Variety, fixed=TRUE)


varieties <- as.character(unique(glymaIDs$Variety))
seqnames <- as.character(unique(glymaIDs$Chromosome))
save(varieties, seqnames, file="ShinyStart.rda")

save(glymaIDs, segranges.df, segranges.df2, chr.summary, file="ChrPlot.rda")
save(glymaIDs, file="GlymaIDs.rda")

########################## Geneological information ############################
tree <- read.csv("soybeanGenealogy.csv", stringsAsFactors=FALSE)
# tree <- melt(tree2, id.vars=c(1, 4:7), measure.vars=2:3, variable.name="parent.type", value.name="parent")
# load("tree.rda")
fulltree <- tree

relatedness <- read.csv("kevinbaconTree.csv", stringsAsFactors=FALSE)
relatedness <- subset(relatedness, (Var1%in%varieties | Var2%in%varieties) & degree>=0 & !is.na(degree))

varlist <- unique(c(relatedness$Var1, relatedness$Var2))

tree <- subset(fulltree, child%in%varlist | parent %in% varlist)

tree$child <- gsub("IA3023x", "IA 3023", tree$child)
tmp <- subset(tree, child=="IA 3023" | parent=="IA 3023")
tmp$child <- gsub("IA 3023", "IA 3023 (NAM)", tmp$child)
tmp$parent <- gsub("IA 3023", "IA 3023 (NAM)", tmp$parent)
tree <- rbind.fill(tree, tmp)
tmp <- subset(tree, child=="Raleigh" | parent=="Raleigh")
tmp$child <- gsub("Raleigh", "NC Raleigh", tmp$child)
tmp$parent <- gsub("Raleigh", "NC Raleigh", tmp$parent)
tree <- rbind.fill(tree, tmp)
tmp <- subset(tree, child=="Roy" | parent=="Roy")
tmp$child <- gsub("Roy", "NC Roy", tmp$child)
tmp$parent <- gsub("Roy", "NC Roy", tmp$parent)
tree <- rbind.fill(tree, tmp)
tmp <- subset(tree, child=="Roy" | parent=="Roy")
tmp$child <- gsub("Roy", "NC Roy", tmp$child)
tmp$parent <- gsub("Roy", "NC Roy", tmp$parent)
tree <- rbind.fill(tree, tmp)
rm(tmp)

save(tree, file="tree.rda")
# rm("tree2")

############################# Yield information ################################
yield <- read.csv("CombinedFieldTrialsData.csv", stringsAsFactors=FALSE)
yield$Cultivar <- gsub("AK", "A.K.", yield$Cultivar)
load("tree.rda")
tree$data <- (tree$child %in% yield$Cultivar) | (tree$child %in% yield$name2)

# Initially, trim yield dataset down to average of important metrics
yield.sum <- ddply(yield, .(Cultivar), summarize, 
                   Year=mean(Year,na.rm=T), 
                   Yield=mean(YieldAdjForMat,na.rm=T), 
                   Protein=mean(AvgProtein,na.rm=T), 
                   Oil=mean(AvgOil,na.rm=T), 
                   Maturity=mean(Maturity,na.rm=T), 
                   Lodging=mean(Lodging,na.rm=T), 
                   SeedSize=mean(SeedSize,na.rm=T), 
                   SeedQuality=mean(SeedQuality,na.rm=T), 
                   PlCount=mean(PlCount), 
                   sequenced=sequenced[1])

# Get data frame with parent names and Cultivar name, as well as parent yield values
parent.data <- tree[,c("child", "parent", "parent.type", "year")]
parent.data$y <- as.numeric(parent.data$parent.type=="mother")/2
parent.data$x <- 0
parent.data <- merge(parent.data, yield.sum, by.x="parent", by.y="Cultivar", 
                     all.x=TRUE, all.y=FALSE)
parent.data$parent.type <- as.character(parent.data$parent.type)
names(parent.data)[2] <- "Cultivar"
# Where parent year is unknown, use imputed year. 
parent.data$Year[is.na(parent.data$Year)] <- parent.data$year[is.na(parent.data$Year)]
parent.data <- parent.data[,-which(names(parent.data)=="year")]
parent.data$sequenced[is.na(parent.data$sequenced)] <- FALSE

parent.data$parent.type <- as.character(parent.data$parent.type)
parent.data$parent[which(is.na(parent.data$parent))] <- "Unknown"

parent.data <- ddply(parent.data, .(Cultivar), function(df) {
  if(sum(df$parent=="Unknown")>=2){
    dfnew <- df[1,]
    dfnew$y <- mean(df$y)
    dfnew$x <- mean(df$x)
  } else {
    dfnew <- df
  } 
  return(dfnew)
})

parent.sub <- subset(parent.data, parent!="Unknown" & !is.na(Yield) & !is.na(Protein) & !is.na(Oil))

yieldplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Yield) & !is.na(Year)), aes(x=Year, y=Yield, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Yield, clickSelects=Cultivar, colour=sequenced), size=3, guide="none") + 
  geom_point(data=yield.sum, aes(x=Year, y=Yield, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 68, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Yield)), aes(x=1965.5+x, y=24+y*6, label=sprintf("Parent - %s: %.2f", parent, Yield), showSelected=Cultivar)) + 
  scale_colour_manual(values=c("red", "blue")) + 
  scale_fill_manual(values=c("red", "blue")) + 
  ggtitle("Yield") + 
  theme_animint(width=400, height=250)

proteinplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Protein) & !is.na(Year)), aes(x=Year, y=Protein, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Protein, clickSelects=Cultivar, colour=sequenced), size=3) +
  geom_point(data=yield.sum, aes(x=Year, y=Protein, showSelected=Cultivar, colour=sequenced), size=5) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  make_text(yield.sum, 1965.5, 39, "Cultivar") + 
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Protein)), aes(x=1965.5+x, y=30+y*1.5, label=sprintf("Parent - %s: %.2f", parent, Protein), showSelected=Cultivar)) + 
  ggtitle("Protein") + 
  theme(legend.position="none") + 
  theme_animint(width=400, height=250)

oilplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Oil) & !is.na(Year)), aes(x=Year, y=Oil, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Oil, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Year, y=Oil, showSelected=Cultivar, colour=sequenced), size=5) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  make_text(yield.sum, 1965.5, 21, "Cultivar") + 
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Oil)), aes(x=1965.5+x, y=16+y/1.5, label=sprintf("Parent - %s: %.2f", parent, Oil), showSelected=Cultivar)) + 
  ggtitle("Oil") + 
  theme(legend.position="none") + 
  theme_animint(width=400, height=250)

proteinOil <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Oil) & !is.na(Protein)), aes(x=Protein, y=Oil, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Protein, y=Oil, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Protein, y=Oil, showSelected=Cultivar, colour=sequenced), size=5) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  make_text(yield.sum, 35, 21, "Cultivar") + 
  geom_text(data=parent.sub, aes(x=35+x, y=16.5+y/2, label=sprintf("Parent - %s: Protein = %.2f     Oil = %.2f", parent, Protein, Oil), showSelected=Cultivar)) + 
  ggtitle("Protein and Oil Content") + 
  theme(legend.position="none") + 
  theme_animint(width=300, height=300)

yieldOil <- ggplot(data=yield.sum) + 
  geom_point(data=subset(parent.data,!is.na(Yield) & !is.na(Oil)), aes(x=Yield, y=Oil, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Yield, y=Oil, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Yield, y=Oil, showSelected=Cultivar, colour=sequenced), size=5) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  make_text(yield.sum, 47, 21, "Cultivar") + 
  geom_text(data=parent.sub, aes(x=47+x, y=16.5+y/2, label=sprintf("Parent - %s: Yield = %.2f     Oil = %.2f", parent, Yield, Oil), showSelected=Cultivar)) + 
  ggtitle("Yield and Oil Content") + 
  theme(legend.position="none") + 
  theme_animint(width=300, height=300)

yieldProtein <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Yield) & !is.na(Protein)), aes(x=Yield, y=Protein, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Yield, y=Protein, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Yield, y=Protein, showSelected=Cultivar, colour=sequenced), size=5) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  make_text(yield.sum, 47, 39, "Cultivar") + 
  geom_text(data=parent.sub, aes(x=47+x, y=30.5+y, label=sprintf("Parent - %s: Yield = %.2f     Protein = %.2f", parent, Yield, Protein), showSelected=Cultivar)) + 
  ggtitle("Yield and Protein Content") + 
  theme(legend.position="none") + 
  theme_animint(width=300, height=300)

maturityplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Maturity) & !is.na(Year)), aes(x=Year, y=Maturity, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Maturity, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Year, y=Maturity, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 938, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Maturity)), aes(x=1965.5+x, y=907+y*4, label=sprintf("Parent - %s: %.2f", parent, Yield), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Maturity")  + 
  theme_animint(width=400, height=250)

maturityplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Maturity) & !is.na(Year)), aes(x=Year, y=Maturity, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Maturity, clickSelects=Cultivar, colour=sequenced), size=3, guide="none") + 
  geom_point(data=yield.sum, aes(x=Year, y=Maturity, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 938, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Maturity)), aes(x=1965.5+x, y=907+y*4, label=sprintf("Parent - %s: %.2f", parent, Maturity), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Maturity")  + 
  theme_animint(width=400, height=250)

lodgingplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(Lodging) & !is.na(Year)), aes(x=Year, y=Lodging, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=Lodging, clickSelects=Cultivar, colour=sequenced), size=3, guide="none") + 
  geom_point(data=yield.sum, aes(x=Year, y=Lodging, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 3.4, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(Lodging)), aes(x=1965.5+x, y=1+y/3, label=sprintf("Parent - %s: %.2f", parent, Lodging), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Lodging")  + 
  theme_animint(width=400, height=250)

seedsizeplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(SeedSize) & !is.na(Year)), aes(x=Year, y=SeedSize, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=SeedSize, clickSelects=Cultivar, colour=sequenced), size=3, guide="none") + 
  geom_point(data=yield.sum, aes(x=Year, y=SeedSize, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 21, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(SeedSize)), aes(x=1965.5+x, y=10+y*1.5, label=sprintf("Parent - %s: %.2f", parent, SeedSize), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Seed Size")  + 
  theme_animint(width=400, height=250)

seedqualityplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(SeedQuality) & !is.na(Year)), aes(x=Year, y=SeedQuality, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=SeedQuality, clickSelects=Cultivar, colour=sequenced), size=3, guide="none") + 
  geom_point(data=yield.sum, aes(x=Year, y=SeedQuality, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 3.3, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(SeedQuality)), aes(x=1965.5+x, y=1.2+y/3, label=sprintf("Parent - %s: %.2f", parent, SeedQuality), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Seed Quality")  + 
  theme_animint(width=400, height=250)


plcountplot <- ggplot() + 
  geom_point(data=subset(parent.data,!is.na(PlCount) & !is.na(Year)), aes(x=Year, y=PlCount, showSelected=Cultivar), alpha=.75, size=5, colour="black", fill="white") +
  geom_point(data=yield.sum, aes(x=Year, y=PlCount, clickSelects=Cultivar, colour=sequenced), size=3) + 
  geom_point(data=yield.sum, aes(x=Year, y=PlCount, showSelected=Cultivar, colour=sequenced), size=5, guide="none") + 
  make_text(yield.sum, 1965.5, 40, "Cultivar") +
  geom_text(data=subset(parent.data, parent!="Unknown" & !is.na(PlCount)), aes(x=1965.5+x, y=10+y*4, label=sprintf("Parent - %s: %.2f", parent, PlCount), showSelected=Cultivar)) + 
  scale_colour_manual("Sequenced", values=c("red", "blue")) + 
  scale_fill_manual("Sequenced", values=c("red", "blue")) + 
  ggtitle("Pl Count")  + 
  theme_animint(width=400, height=250)

animint2dir(plot.list=list(yield=yieldplot, 
                          protein=proteinplot, 
                          oil=oilplot, 
                          maturity=maturityplot,
                          lodging=lodgingplot,
                          seedsize=seedsizeplot,
                          seedquality=seedqualityplot,
                          plcount=plcountplot,
                          yieldOil=yieldOil, 
                          yieldProtein=yieldProtein, 
                          proteinOil=proteinOil),
           out.dir="./www/", open.browser=FALSE)

save.image(file="ShinyData.RData")