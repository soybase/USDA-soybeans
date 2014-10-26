library(stringr)

fieldtrials <- read.csv("CombinedFieldTrialsData.csv", header=T, stringsAsFactors=F)
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Williams 82", "Williams82")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Amsoy 71", "Amsoy")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "IA 3023", "IA3023")


load("tree-large.rda")
fulltree <- tree
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
load("tree.rda")


load("ShinyStart.rda")
# Get years, etc from full tree
tree.fixed <- subset(fulltree, fulltree$child %in% varieties | fulltree$child %in% fieldtrials$Cultivar)
tmp <- nrow(tree.fixed)-1
while(nrow(tree.fixed)>tmp){
  # recursively search for ancestors
  tmp <- nrow(tree.fixed)
  tree.fixed <- unique(rbind(tree.fixed, subset(fulltree, fulltree$parent %in% tree.fixed$child)))
}

tree2 <- merge(tree.fixed, unique(fieldtrials[,c("Cultivar", "Year", "AvgYield", "MG")]), by.x=c("child"), by.y=c("Cultivar"), all.x=T, all.y=F)
tree2$year <- pmin(tree2$year, tree2$Year, na.rm=T)
tree2 <- unique(tree2[,-which(names(tree2)%in%c("Year", "year.imputed", "min.repro.year", "yield"))])


tree2$geneticdata <- tree2$child%in%varieties

sum(varieties%in%tree2$child)
