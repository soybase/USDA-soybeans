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

library(phyViz)
ugh <- melt(ped.table, id.vars=c("indID", "famID", "gender", "phenotype", "geneticdata"), measure.vars=c("fatID", "matID"))
ugh[,1] <- as.character(ugh[,1])
names(ugh) <- c("child", "famID", "gender", "phenotype", "geneticdata", "variable", "parent")
ugh$parent[which(ugh$parent=="0")] <- NA

# remove varieties that aren't ancestors or descendants of varieties with genetic data
varlist <- unique(unlist(sapply(varieties, function(i) as.character(phyViz:::nodeToDF(phyViz::buildDesList(i, ugh))$label))))
varlist <- unique(c(varlist, unlist(sapply(varieties, function(i) as.character(phyViz:::nodeToDF(phyViz::buildAncList(i, ugh))$label)))))
ped.table <- subset(ped.table, indID%in%varlist)

# replace NA phenotype with -9
ped.table$phenotype[is.na(ped.table$phenotype)] <- -9

# replace missing parental IDs with unique strings
idx <- which(is.na(ped.table$fatID) & !is.na(ped.table$matID))
ped.table[idx,"fatID"] <- ped.table[idx,"matID"]
ped.table[idx,"matID"] <- NA
ped.table[is.na(ped.table$matID),"matID"] <- paste0("NA_", 1:sum(is.na(ped.table$matID)))

write.table(ped.table[,c("famID", "indID", "fatID", "matID", "gender", "phenotype")], file="~/Documents/Rprojects/Soybeans/snp/FullSoybeanYield.ped", sep="\t")
tmp <- ped.table[,c("famID", "indID", "fatID", "matID", "gender", "phenotype")]
# names(tmp) <- c("FID", "IID", "yield")

fixids <- function(x){
  x <- str_replace(x, fixed("5601T"), "901-G04_RS_5601T")
  x <- str_replace(x, fixed("Bonus"), "901-G06_RS_Bonus")
  x <- str_replace(x, "^Shelby", "RS_Shelby")
  x <- str_replace(x, "^A3127", "RS_A3127")
  x <- str_replace(x, fixed("S-100"), "RS_S-100")
  x <- str_replace(x, fixed("A.K. (Harrow)"), "AK_004")
  x <- str_replace(x, fixed("Clark 63"), "Clark_NuGEN")
  x <- str_replace(x, fixed("Raleigh"), "NCRaleigh")
  x <- str_replace_all(x, " ", "_")
  return(x)
}

tmp$indID <- fixids(tmp$indID)
tmp$fatID <- fixids(tmp$fatID)
tmp$matID <- fixids(tmp$matID)

# write.table(tmp[,1:4], file="~/Documents/Rprojects/Soybeans/snp/soybeanParents", sep="\t", row.names = F, quote = F)
# tmp$genotype <- 0
# write.table(tmp[,1:6], file="~/Documents/Rprojects/Soybeans/snp/soybeanYield.tfam", sep="\t", row.names = F, col.names=F, quote = F)

ndummies <- length(which(is.na(tmp$fatID)))
tmp$fatID[which(is.na(tmp$fatID))] <- paste0("DUMMY", 1:ndummies)
tmp$matID[which(is.na(tmp$matID))] <- paste0("DUMMY", (ndummies + 1):(ndummies + length(which(is.na(tmp$matID)))))
tmp$famID <- tmp$indID

write.table(tmp, file="~/Documents/Rprojects/Soybeans/snp/soybeanYieldInfo.ped", sep="\t", row.names = F, col.names=F, quote = F)

