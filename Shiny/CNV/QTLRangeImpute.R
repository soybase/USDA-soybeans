qtl <- read.table("./QTL_MAG.txt", sep="\t", header=T)

library(reshape2)
library(dplyr)

qtlrange <- qtl %>% group_by(Chromosome, QTLID, QTLName) %>% summarize(start=min(Start_bp), end=max(End_bp), LocusName=paste(unique(LocusName), collapse=";"), num=length(QTLID)) %>% arrange(Chromosome, start, end)

impute.start <-  function(start, end, num){
  tmp <- subset(data.frame(start=start, end=end, num=num), num==1)
  if(nrow(tmp)==0) return(start)
  startidx <- sapply(tmp$start, function(i) which(sort(unique(tmp$start))==i))
  tmp$prev.end <- c(min(tmp$start), sort(unique(tmp$end)))[startidx]
  imp.start <- start
  imp.start[which(num==1)] <- 3*tmp$start/4 + tmp$prev.end/4
  imp.start
}

impute.end <- function(start, end, num){
  tmp <- subset(data.frame(start=start, end=end, num=num), num==1)
  if(nrow(tmp)==0) return(end)
  endidx <- sapply(tmp$end, function(i) which(sort(unique(tmp$end))==i))
  tmp$last.start <- c(sort(unique(tmp$start)), max(tmp$end))[endidx+1]
  imp.end <- end
  imp.end[which(num==1)] <- 3*tmp$end/4+tmp$last.start/4
  imp.end
}

qtlrange.imputed <- qtlrange %>% transform (imp.start = impute.start(start, end, num), imp.end = impute.end(start, end, num), tmp=(1:length(QTLID)))

qtlrange.imputed$Chromosome <- gsub("Gm", "Chr", qtlrange.imputed$Chromosome)

glymaIDmatch <- function(Start, End, chromosome){
  tmp <- unique(filter(qtlrange.imputed, Start>end & End>start & Chromosome==chromosome)$QTLID)
  if(length(tmp)==0){
    tmp <- NA
  } else {
    tmp <- paste(tmp, collapse="; ")
  }
  
  tmp
}

glymaIDs2 <- glymaIDs %>% group_by(Chromosome, ID, Start, End) %>% mutate(qtlMatch = glymaIDmatch(Start, End, Chromosome))



# 
# qplot(data=subset(one.point.qtl, num==1), geom="errorbar", x=tmp, ymin=start, ymax=end, alpha=I(.5)) + 
#   geom_errorbar(data=subset(one.point.qtl, num==1), aes(ymax=imp.end, ymin=imp.start, x=tmp), alpha=I(.5))
