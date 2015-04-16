# This function plots the family tree of one of more varieties, with overlapping branches highlighted

# # Filter family tree by varieties for which we have genetic data
# 
# load("FullTree.rda")
# 
# library(doMC)
# 
# vars <- mclapply(varieties, function(i){
#   buildAncDesTotalDF(i, tree, 8, 0)[,c("variety", "label", "gen", "type")]
# })
# vars <- do.call("rbind", vars)
# 
# subtree <- subset(tree, child%in%vars$label | child%in%vars$variety)
# tree <- subtree
# save(tree, file="SubTree.rda")

library(ggenealogy)
load("SubTree.rda")
load("ShinyStart.rda")

# treeGraph <- treeToIG(tree, isDirected=FALSE)

plotFamilyTree <- function(variety, gens, ncol=2){
  gen.vars2 <- variety
  vals <- list()
  vals$gen.vars <- gen.vars2
  
  if(length(variety)>0){
    temp2 <- lapply(gen.vars2, function(i){
      if(i %in% tree$child | i %in% tree$parent){
        temp <- buildAncDesTotalDF(i, tree, mAnc=gens, mDes=gens)
        temp <- subset(temp, !(leaf & !label%in%varieties))
        # Subset temp to get accurate number of generations according to the slider
        #           temp <- subset(temp, abs(gen)<=input$gens)
        # get plot coords of subsetted data frame
        temp$label2 <- temp$label
      } else {
        temp <- data.frame(variety=i, label=i, root=NA, root.gen=0, gen=0, type=0, branch=0, x=0, y=0, xstart=0, ystart=0, xend=0, yend=0, branchx=0, branchy=0, id=0, par.id=0, size=6, label2=paste(i, "-", "Query not found in database", sep=""))
      }
      unique(temp)
    }) %>% rbind_all()
    
    temp2$color <- c("grey", "#000000")[temp2$label2%in%varieties + 1]
    
    vals$df <- temp2
    
    temp <- merge(data.frame(variety=variety,NewName=vals$gen.vars), summarize(vals$df %>% group_by(label), gen=mean(gen*c(-1,1)[(type=="descendant")+1])), by.x=2, by.y=1)
    vals$match <- temp[order(temp$gen, temp$variety),]
  } else {
    vals$df <- data.frame()
    vals$match <- data.frame()
    vals$gen.vars <- variety
  }
  
  if(nrow(vals$df)>0){
    plot <- ggplot(data=vals$df) + 
      geom_text(aes(x=x, y=y, label=label2, colour=color, size=size), vjust=-.25, hjust=.5) +
      geom_segment(data=vals$df, aes(x=xstart, xend=xend, y=ystart, yend=yend),inherit.aes=F) + 
      geom_segment(data=vals$df, aes(x=xend, xend=branchx, y=yend, yend=branchy),inherit.aes=F) +
      facet_wrap(~variety, scales="free", ncol=ncol) +
      scale_size_continuous(range=c(4,6),guide="none") +
      scale_colour_identity() +
      theme_bw() +
      theme(axis.title=element_blank(), 
            axis.text=element_blank(), 
            axis.ticks=element_blank()) + 
      scale_x_continuous(expand = c(.1, 1.075)) + 
      scale_y_continuous(expand = c(.1, 1.075))
  } else {
    plot <- ggplot() + 
      geom_text(aes(x=0, y=0, label="Please select varieties\n\n Note: It may take a minute to process the input")) +         
      theme_bw() + 
      theme(axis.text=element_blank(), 
            axis.ticks=element_blank(), 
            axis.title=element_blank())
  }
  return(plot)
}

