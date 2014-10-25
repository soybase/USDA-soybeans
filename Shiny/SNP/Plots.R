load("SNPCombinedData.RData")

library(ggplot2)
library(dplyr)
library(reshape2)


tmp.df <- filter(tmp.df, !grepl("scaffold", chr))
tmp.df$Alt.Allele.Count <- round(tmp.df$Alt.Allele.Freq)


range <- c(0, 50000)

qplot(data=subset(tmp.df, chr=="Chr02" & pos>=range[1] & pos<=range[2]), x=factor(pos), y=..count.., position="stack", fill=alt, geom="histogram", stat="bin", binwidth=1)


chr.pos <- group_by(tmp.df, chr, pos)
snps <- summarize(chr.pos, 
                 ref=unique(ref), 
                 A=sum(Alt.Allele.Count*(alt=="A")), 
                 G=sum(Alt.Allele.Count*(alt=="G")), 
                 C=sum(Alt.Allele.Count*(alt=="C")), 
                 T=sum(Alt.Allele.Count*(alt=="T")),
                 total=sum(Alt.Allele.Count))
ref <- cbind(A=(snps$ref=="A")*(20-snps$total), 
             G=(snps$ref=="G")*(20-snps$total),
             C=(snps$ref=="C")*(20-snps$total),
             T=(snps$ref=="T")*(20-snps$total))

snps[,4:7] <- snps[,4:7] + ref

snps$total <- 20
snps <- snps[,1:7]


snp.counts <- melt(snps, id.vars=1:3, value.name="Count", variable.name="Nucleotide")

library(RColorBrewer)
pal <- brewer.pal(8, "Paired")[c(1, 2, 7, 8)]
qplot(data=subset(snp.counts, chr=="Chr18" & pos>=57810000 & pos < 57830000), x=factor(pos), position="stack", fill=Nucleotide, geom="histogram", weight=Count) + 
  scale_fill_manual("Nucleotide", values=c("A"=pal[1], "G"=pal[2], "T"=pal[3], "C"=pal[4])) + 
  xlab("Position on Chr 18") + 
  ylab("Total number of SNPs in 10 lines") + 
  theme(axis.text.x=element_text(angle=90))



library(entropy)

tmp.df$bin <- floor(tmp.df$pos/10000)*10000

chr.pos.var <- group_by(tmp.df, chr, pos, Variety)
snps.var <- summarize(chr.pos.var, 
                  ref=unique(ref), 
                  A=sum(Alt.Allele.Count*(alt=="A")), 
                  G=sum(Alt.Allele.Count*(alt=="G")), 
                  C=sum(Alt.Allele.Count*(alt=="C")), 
                  T=sum(Alt.Allele.Count*(alt=="T")),
                  total=sum(Alt.Allele.Count))
ref.var <- with(snps.var, 
                cbind(A=(ref=="A")*(2-total), 
                      G=(ref=="G")*(2-total),
                      C=(ref=="C")*(2-total),
                      T=(ref=="T")*(2-total)))

snps.var[,5:8] <- snps.var[,5:8] + ref.var



varsnps$bin <- floor(snps$pos/1000)*1000
snps <- snps %>% group_by(chr, bin) 
snps.summary <- snps %>% summarize(mi.empirical(.[4:7]))
