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


setwd("~/Documents/R Projects/Soybeans/Shiny/Genealogy")

########################## Geneological information ############################
tree <- read.csv("soybeanGenealogy.csv", stringsAsFactors=FALSE)

# tree$child <- gsub("IA3023x", "IA3023", tree$child)
save(tree, file="tree.rda")

varieties <- unique(c(tree$child, tree$parent))
varieties <- varieties[!grepl(" x ", varieties)]
save(varieties, file="varieties.rda")

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

################################# Kevin Bacon ##################################

source("./kevinbacon.R")
save(kevinbacon, treeGraph, file="kevinbacon.rda")
