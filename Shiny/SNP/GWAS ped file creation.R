setwd("~/Documents/Rprojects/USDAsoybeans/Shiny/SNP/")
source("GWAS Year Database Start.R")

library(plyr)
library(dplyr)
library(reshape2)
tree2 <- tree2[,-which(names(tree2)%in%c("src", "parent.type"))] %>% group_by(child) %>% unique()
tree2 <- ddply(tree2, .(child), transform, parent.type=1:length(parent), nparents=length(parent))
tree2 <- subset(tree2, nparents<3)
tree2 <- ddply(tree2[,-which(names(tree2)%in%c("nparrows", "nparents"))], .(child), transform, parent.type=1:length(parent))

varlist <- unique(c(tree2$child, tree2$parent))
varlist <- sort(varlist[!is.na(varlist)])

IDlist <- data.frame(Name = as.character(varlist), ID = as.numeric(factor(varlist, levels=varlist)))
oldped <- read.csv("soybean_ped.csv", header=T, stringsAsFactors=F)
# tree2$child <- as.numeric(factor(tree2$child, levels=varlist))
# tree2$parent <- as.numeric(factor(tree2$parent, levels=varlist))
treecast <- dcast(tree2, child ~ parent.type, value.var="parent")
names(treecast)[2:3] <- c("parent1", "parent2")
treecast$founder <- rowSums(is.na(treecast[,2:3]))==2
treecast$parent1[treecast$founder] <- 0
treecast$parent2[treecast$founder] <- 0
treecast[treecast$child%in%oldped$child, "fatID"] <- oldped[oldped$child%in%treecast$child, "father"]
treecast[treecast$child%in%oldped$child, "matID"] <- oldped[oldped$child%in%treecast$child, "mother"]

tree2$yield <- tree2$AvgYield
tree2$yield[is.na(tree2$yield) & !is.na(tree2$PredYield)] <- tree2$PredYield[is.na(tree2$yield) & !is.na(tree2$PredYield)]

ped.table <- data.frame(famID = 1, indID=treecast$child, fatID=treecast$parent1, matID=treecast$parent2, gender=0)
ped.table <- merge(ped.table, unique(tree2[,c("child", "yield", "geneticdata")]), by.x="indID", by.y="child")
names(ped.table)[which(names(ped.table)=="yield")] <- "phenotype"
sanitycheck <- merge(IDlist[IDlist$ID%in%subset(ped.table, geneticdata)$indID,], 
                     subset(ped.table, geneticdata), by.x="ID", by.y="indID")

ped.table$geneticdata[which(is.na(ped.table$geneticdata))] <- subset(IDlist, Name%in%ped.table$indID[which(is.na(ped.table$geneticdata))])$Name%in%varieties

sum(!is.na(ped.table$phenotype))
sum(ped.table$geneticdata)
sum(!is.na(ped.table$phenotype)&ped.table$geneticdata)

# replace NA phenotype with -9
ped.table$phenotype[is.na(ped.table$phenotype)] <- -9

# replace missing parental IDs with unique strings
idx <- which(is.na(ped.table$fatID) & !is.na(ped.table$matID))
ped.table[idx,"fatID"] <- ped.table[idx,"matID"]
ped.table[idx,"matID"] <- NA
ped.table[is.na(ped.table$matID),"matID"] <- paste0("NA_", 1:sum(is.na(ped.table$matID)))

write.table(ped.table[,c("famID", "indID", "fatID", "matID", "gender", "phenotype")], file="soybeanYield.ped", sep="\t")


tmp <- read.table("~/Documents/Rprojects/Soybeans/snp/MeanPlus0p150SD.ped", sep="\\s")
