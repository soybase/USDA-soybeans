library(stringr)

fieldtrials <- read.csv("CombinedFieldTrialsData.csv", header=T, stringsAsFactors=F)
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Williams 82", "Williams82")
# fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Amsoy 71", "Amsoy")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "IA 3023", "IA3023")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "506-13640-2", "S06-13640-2")
fieldtrials$Cultivar2 <- fieldtrials$Cultivar
idx <- fieldtrials$Source=="Pioneer" & nchar(as.character(fieldtrials$name2)==4)
fieldtrials$Cultivar[idx] <- with(fieldtrials[idx,], paste(Source, name2))
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, fixed("AK (Harrow)"), fixed("A.K. (Harrow)"))


load("tree-large.rda")
fulltree <- tree
fulltree$parent <- gsub("506-13640-2", "S06-13640-2", fulltree$parent)
fulltree$child <- gsub("506-13640-2", "S06-13640-2", fulltree$child)
fulltree$parent <- gsub("IA 3023", "IA3023", fulltree$parent)
fulltree$child <- gsub("IA 3023", "IA3023", fulltree$child)
fulltree$parent <- gsub("Williams 82", "Williams82", fulltree$parent)
fulltree$child <- gsub("Williams 82", "Williams82", fulltree$child)
fulltree$parent <- gsub("S ?100", "S-100", fulltree$parent)
fulltree$child <- gsub("S ?100", "S-100", fulltree$child)
fulltree$parent <- gsub("Gasoy ?17", "Gasoy17", fulltree$parent)
fulltree$child <- gsub("Gasoy ?17", "Gasoy17", fulltree$child)
fulltree$parent <- gsub("No\\. 18", "180501", fulltree$parent)
fulltree$child <- gsub("No\\. 18", "180501", fulltree$child)
fulltree$parent <- gsub("No\\. 18", "180501", fulltree$parent)
fulltree$child <- gsub("No\\. 18", "180501", fulltree$child)
fulltree$parent <- gsub("PI 88\\.788", "88788", fulltree$parent)
fulltree$child <- gsub("PI 88\\.788", "88788", fulltree$child)
fulltree$parent <- gsub("P9381", "Pioneer 9381", fulltree$parent)
fulltree$year[fulltree$child=="IA3023"] <- 2003
fulltree$year.imputed[fulltree$child=="IA3023"] <- FALSE
load("tree.rda")


maturity <- read.csv("MaturityGroup.csv")
maturity$Variety <- gsub("Gasoy ?17", "Gasoy17", maturity$Variety)
maturity$Variety <- gsub("Williams 82", "Williams82", maturity$Variety)
names(maturity)[3] <- "MaturityGroup"

load("ShinyStart.rda")
# Get years, etc from full tree
tree.fixed <- subset(fulltree, fulltree$child %in% varieties | fulltree$child %in% fieldtrials$Cultivar)
tmp <- nrow(tree.fixed)-1
while(nrow(tree.fixed)>tmp){
  # recursively search for ancestors
  tmp <- nrow(tree.fixed)
  tree.fixed <- unique(rbind(tree.fixed, subset(fulltree, fulltree$parent %in% tree.fixed$child)))
}
varieties <- str_replace(str_replace(varieties, "506-13640-2", "S06-13640-2"), "NCRoy", "Roy")

tree2 <- merge(tree.fixed, unique(fieldtrials[,c("Cultivar", "Year", "AvgYield", "MG")]), by.x=c("child"), by.y=c("Cultivar"), all.x=T, all.y=F)
tree2$year <- pmin(tree2$year, tree2$Year, na.rm=T)
tree2 <- unique(tree2[,-which(names(tree2)%in%c("Year", "year.imputed", "min.repro.year", "yield"))])


tree2$geneticdata <- tree2$child%in%varieties
tree2 <- merge(tree2, maturity[,2:3], by.x="child", by.y="Variety", all.x=T, all.y=T)
tree2$MaturityGroup[!is.na(tree2$MG) & is.na(tree2$MaturityGroup)] <- tree2$MG[!is.na(tree2$MG) & is.na(tree2$MaturityGroup)]

tree2 <- tree2[,-which(names(tree2)=="MG")]

# Convert AvgYield to kg/ha
# Reference: http://www.ncagr.gov/agronomi/obt419.htm
tree2$AvgYield <- tree2$AvgYield * 0.067 * 1000

predmg <- plyr::rbind.fill(
            data.frame(MaturityGroup=2, slope1=9.05, int1=-15040, slope2=31.1, int2=2770, change=1968),
            data.frame(MaturityGroup=3, slope1=12.4, int1=-21678, slope2=29.4, int2=2676, change=1964),
            data.frame(MaturityGroup=4, slope1=12.7, int1=-22693, slope2=26.5, int2=2339, change=1971)
          )

tree2$PredYield <- NA
tmp <- tree2$year[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) &  tree2$MaturityGroup==2]
tree2$PredYield[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) &  tree2$MaturityGroup==2] <- with(predmg[1,], I(tmp<=change)*(slope1*tmp + int1 ) + I(tmp>change)*(slope2*(tmp-change)+int2))
tmp <- tree2$year[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) &  tree2$MaturityGroup==3]
tree2$PredYield[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) &  tree2$MaturityGroup==3] <- with(predmg[2,], I(tmp<=change)*(slope1*tmp + int1 ) + I(tmp>change)*(slope2*(tmp-change)+int2))
tmp <- tree2$year[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) & tree2$MaturityGroup==4]
tree2$PredYield[!is.na(tree2$MaturityGroup) & !is.na(tree2$year) &  tree2$MaturityGroup==4] <- with(predmg[3,], I(tmp<=change)*(slope1*tmp + int1 ) + I(tmp>change)*(slope2*(tmp-change)+int2))

write.csv(subset(tree2, !is.na(yield)), "~/Dropbox/Shoemaker_Specht/yield-analysis/PredictedYieldTreeDB.csv", row.names=F)
# ggplot2::qplot(data=unique(tree2[,c("child", "AvgYield", "PredYield")]), x=AvgYield, y=PredYield, geom="point")
