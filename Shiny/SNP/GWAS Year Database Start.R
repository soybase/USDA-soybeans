library(stringr)

fieldtrials <- read.csv("CombinedFieldTrialsData.csv", header=T, stringsAsFactors=F)
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Williams 82", "Williams82")
fieldtrials$Cultivar <- str_replace(fieldtrials$Cultivar, "Amsoy 71", "Amsoy")


load("tree.rda")

tree2 <- merge(tree, unique(fieldtrials[,c("Cultivar", "Year", "AvgYield")]), by.x=c("child"), by.y=c("Cultivar"), all.x=T, all.y=F)
tree2$year <- pmin(tree2$year, tree2$Year, na.rm=T)
tree2 <- unique(tree2[,-which(names(tree2)%in%c("Year", "year.imputed", "min.repro.year", "yield"))])

