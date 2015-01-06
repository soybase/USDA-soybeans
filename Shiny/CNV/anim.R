library(ggplot2)
library(animint)

genome <- chr.summary
genome$end <- cumsum(genome$end)
genome$start <- 1 + c(0, genome$end[1:19])
genome$mid <- (genome$start + genome$end)/2

chr.summary$mid <- (chr.summary$start  + chr.summary$end)/2

segranges.df.sum <- ddply(segranges.df, .(seqnames), transform, 
                          start=start + genome$start[as.numeric(gsub("Chr", "", seqnames))]-1, 
                          end=end + genome$start[as.numeric(gsub("Chr", "", seqnames))]-1)
segranges.df.sum$midpoint <- (segranges.df.sum$start + segranges.df.sum$end)/2

overview.dens <- ddply(segranges.df.sum, .(seqnames), function(df) {
  res <- density(df$midpoint, adjust=.15, cut=T)
  res2 <- data.frame(seqnames=unique(df$seqnames), x=res$x, y=res$y)
  return(res2)
})
overview.dens$y <- overview.dens$y/max(overview.dens$y)


overview <- 
  ggplot() + 
  geom_tallrect(data=genome, aes(xmin=start, xmax=end, clickSelects=seqnames), alpha=.25, fill="grey80") + 
  geom_area(data=overview.dens, aes(x=x, y=y, group=seqnames, colour=seqnames, fill=seqnames), inherit.aes=F) + 
  scale_x_continuous("", breaks=genome$mid, labels=seqnames, minor_breaks=genome$end) + 
  scale_colour_discrete(guide="none") + 
  scale_fill_discrete(guide="none") + 
  ylab("Scaled Density") + 
  ggtitle("Distribution of CNVs (all varieties)") + 
  theme_animint(width=800, height=400) 

chromosomeplot <- 
  ggplot() + 
  geom_density(data=segranges.df, aes(x=midpoint, y=..scaled.., showSelected=seqnames), adjust=.05, trim=T) + 
  geom_text(data=chr.summary, 
            aes(label=sprintf("Density of CNVs along chromosome %s", gsub("Chr0?", "", seqnames)), 
                x=mid, y=1, showSelected=seqnames)) + 
  scale_x_continuous("", limits=c(0, max(chr.summary$end))) + 
  scale_y_continuous("Scaled Density", limits=c(0, 1.05)) + 
  theme_animint(width=800, height=400)

animint2dir(list(overview=overview, chromosomeplot=chromosomeplot))
