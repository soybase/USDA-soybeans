
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)


load("tree.rda")
# tree is a data frame with columns child, parent, year, year.imputed (not necessary...yet), min.repro.year (not necessary...yet). It describes edges of the tree connecting parent and child nodes.

# Returns the parents of a particular variety (if they exist)
getparent <- function(bean){
  subset(tree, child==bean)$parent
}

# Returns children of a particular variety (if they exist)
getchild <- function(bean){
  subset(tree, parent==bean)$child
}

# Recursive function to return a list-style tree with all ancestors of a variety
getancestors <- function(bean, gen=0, maxgen=Inf){
  if(is.na(bean)) return()
  
  temp <- getparent(bean)
  if(length(temp)==0) return()
  
  res <- lapply(temp[!is.na(temp)], function(i){
    # print(i)
    if(gen+1<maxgen) {
      temp2 <- getancestors(i, gen=gen+1, maxgen=maxgen)
    } else {
      temp2 <- NULL
    }
    if(length(temp2)<1) return(list(label=i, root=bean, root.gen=gen, gen=gen+1, type="ancestor"))
    return(c(label=i, root=bean, root.gen=gen, gen=gen+1, type="ancestor", temp2))
  })
  
  if(gen==0){
    return(c(label=bean, root=bean, root.gen=gen, gen=gen, type="ancestor", res))
  } else {
    return(res)
  } 
}

# recursive function to return a list-style tree of all descendants of a particular variety
getdescendants <- function(bean, gen=0, maxgen=Inf){
  if(is.na(bean)) return()
  
  temp <- getchild(bean)
  if(length(temp)==0) return()
  
  res <- lapply(temp[!is.na(temp)], function(i){
    # print(i)
    if(gen+1<maxgen){
      temp2 <- getdescendants(i, gen=gen+1, maxgen=maxgen)
    } else {
      temp2 <- NULL
    }
    
    if(length(temp2)<1) return(list(label=i, root=bean, root.gen=gen, gen=gen+1, type="descendant"))
    return(c(label=i, root=bean, root.gen=gen, gen=gen+1, type="descendant", temp2))
  })
  
  if(gen==0){
    return(c(label=bean, root=bean, root.gen=gen, gen=gen, type="descendant", res))
  } else{
    return(res)
  } 
}

# getparent("Adams")
# getchild("Adams")
# 
# getancestors("Amsoy")
# getancestors("Lawrence")
# getdescendants("Kent")

id.offset <<- 1

# converts the list-style-tree to a data frame, where each variety has an id value 
# and references its' parent's id value. ID value ranges correspond to generation. 
# It is possible that with more complex trees the range of id values may need to 
# expand to reduce the probability of two varieties being assigned the same id value. 
node.to.data.frame <- function(tlist, branch=0, par.id = NA){
  if(is.null(tlist)) return(data.frame())
  
  listidx <- which(sapply(tlist, mode)=="list")

  if(length(listidx)==0){
    temp <- as.data.frame(tlist)
    if(nrow(temp)==0) return(data.frame())
    if(!"gen"%in%names(temp)){
      id.offset <<- id.offset+1
      return(cbind(as.data.frame(tlist), branch=branch, par.id=par.id, id=sample(0:99, 1)*10+id.offset/10))
    }
    id.offset <<- id.offset+1
    return(cbind(as.data.frame(tlist), branch=branch, par.id=par.id, id=sign(temp$gen)*sample((abs(temp$gen)*100):((abs(temp$gen)+1)*100-1), 1)*10+id.offset/10))
  } else {
    temp <- as.data.frame(tlist[-listidx])
    branchidx <- listidx-min(listidx)+1
     if(length(branchidx)>1){
      branchidx <- seq(-.5, .5, length.out=length(branchidx)) 
     } else branchidx <- c(-.5, .5)[temp$gen%%2+1]
    id.offset <<- id.offset+1
    id <- sign(temp$gen)*sample((abs(temp$gen)*100):((abs(temp$gen)+1)*100-1), 1)*10+id.offset/10
    return(rbind.fill(cbind(temp, branch=branch, id=id, par.id=par.id), 
                      ldply(1:length(listidx), function(i) 
                        node.to.data.frame(tlist[[listidx[i]]], branch=branchidx[i], par.id=id))))
  }
}

# node.to.data.frame(getdescendants("Kent"))
# node.to.data.frame(getancestors("Kent"))


# Calculates coordinates to plot each node on a tree (with corresponding edges). 
# x, y describe the coordinates of the label
# the line from (xstart,  ystart) to (xend, yend) is the node (the horiz. line by the label)
# the line from (xend, yend) to (branchx, branchy) is the edge
plotcoords <- function(df){
  if(nrow(df)==0) return(data.frame())

  if(nrow(subset(df, root.gen==0 & gen==0))>1){
    temp <- subset(df, root.gen==0 & gen==0)
    temp$type <- "center"
    old.ids <- temp$id[-1]
    df$par.id[which(df$par.id%in%old.ids)] <- temp$id[1]
    temp$id[-1] <- temp$id[1]
    temp <- unique(temp)
    df <- rbind.fill(subset(df, !(root.gen==0 & gen==0)), temp)
  }
  
  df$leaf <- !df$id%in%df$par.id
  
  
  # initialize x, y coords
  df$genside <- df$gen*(df$type=="descendant")-df$gen*(df$type=="ancestor")
  df$x <- 0
  df$y <- 0
  
  # ensure y coords are spread among all branches
  df$y[df$leaf & df$genside<0] <- seq(1, sum(df$leaf), length.out=sum(df$leaf & df$genside<0))
  df$y[df$leaf & df$genside>0] <- seq(1, sum(df$leaf), length.out=sum(df$leaf & df$genside>0))
  
  # ensure a single branch on one side is centrally placed
  if(sum(df$leaf & df$genside<0)==1) df$y[df$leaf & df$genside<0] <- sum(df$leaf)/2
  if(sum(df$leaf & df$genside>0)==1) df$y[df$leaf & df$genside>0] <- sum(df$leaf)/2
  
  df$xstart <- 0
  df$ystart <- 0
  df$xend <- 0
  df$yend <- 0
  df$branchx <- 0
  df$branchy <- 0
  
  # sort data frame
  df <- df[rev(order(df$type, df$root.gen, df$gen, df$par.id, df$branch)),]
  
  # set y coordinates
  for(i in unique(df$genside)[rev(order(abs(unique(df$genside))))]){
    if(i!=0){
      kids <- subset(df, df$genside==i)
      for(j in unique(kids$par.id)){
        df$y[df$id==j] <- mean(subset(kids, par.id==j)$y)
      }
    } else { # center coordinate needs to account for both sides of the tree
      df$y[df$gen==0] <- mean(subset(df, is.na(par.id))$y)
    }
    
  }
  
  yfac <- diff(range(df$y))
  
  df$ystart <- df$y
  df$yend <- df$y

  
  widths <- ddply(df, .(genside), summarize, len=max(nchar(as.character(label))))
  
  widths$len <- widths$len*(1+yfac/(1+yfac)) # padding for label length
  
  gap <- mean(widths$len)/2 # set x gap for connecting branch
  
  for(i in unique(df$genside)[order(abs(unique(df$genside)))]){
    j <- which(df$genside==i)
    if(i == 0){ # Center of the two trees    
      df$x[j] <- 0
      df$xstart[j] <- -widths$len[widths$genside==i]/2
      df$xend[j] <- widths$len[widths$genside==i]/2
      df$branchx[j] <- df$xend[j] # make the "branch" for this node nonexistent
      df$branchy[j] <- df$yend[j]
    } else if(i==1){ # positive side of tree - first "reversed" branch
      for(k in j){
        par <- df[which(df$id==df$par.id[k]),] # find parent
        # only the center node has no "parents" 
        df$branchx[k] <- par$xend
        df$branchy[k] <- par$yend
        df$xend[k] <- df$branchx[k]+gap*sign(i)
        df$xstart[k] <- df$xend[k]+widths$len[widths$genside==i]*sign(i)
        df$x[k] <- (df$xend[k]+df$xstart[k])/2
      } 
    } else if(i>1){ # positive side of tree - "reversed" branch
      for(k in j){
        par <- df[which(df$id==df$par.id[k]),] # find parent
        # only the center node has no "parents" 
        df$branchx[k] <- par$xstart
        df$branchy[k] <- par$ystart
        df$xend[k] <- df$branchx[k]+gap*sign(i)
        df$xstart[k] <- df$xend[k]+widths$len[widths$genside==i]*sign(i)
        df$x[k] <- (df$xend[k]+df$xstart[k])/2
      } 
    } else { # negative side of tree - branch is more positive than node
      for(k in j){
        par <- df[which(df$id==df$par.id[k]),] # find parent
        # only the center node has no "parents" 
        df$branchx[k] <- par$xstart
        df$branchy[k] <- par$ystart
        df$xend[k] <- df$branchx[k]+gap*sign(i)
        df$xstart[k] <- df$xend[k]+widths$len[widths$genside==i]*sign(i)
        df$x[k] <- (df$xend[k]+df$xstart[k])/2
      } 
    }
  }
 
  df$size <- (1-(nrow(df))/412)*4+2
  return(df)
}

# 
# library(ggplot2)
# test <- plotcoords(node.to.data.frame(getancestors("Skylla")))
# qplot(data=test, x=x, y=y, label=label, geom="text", vjust=0, hjust=.5) + 
#   geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend)) + 
#   geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy))
# 
# test <- plotcoords(node.to.data.frame(getdescendants("Skylla")))
# qplot(data=test, x=x, y=y, label=label, geom="text", vjust=0, hjust=.5) + 
#   geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend)) + 
#   geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy))
# 
# test <- plotcoords(node.to.data.frame(getancestors("Williams")))
# qplot(data=test, x=x, y=y, label=label, geom="text", vjust=0, hjust=.5) + 
#   geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend)) + 
#   geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy)) +
#   theme_bw()
# 
# 
# test <- plotcoords(rbind(node.to.data.frame(getancestors("Skylla")), node.to.data.frame(getdescendants("Skylla"))))
# qplot(data=test, x=x, y=y, label=label, geom="text", vjust=0, hjust=.5) + 
#   geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend)) + 
#   geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy)) +
#   theme_bw()
# 
# test <- plotcoords(rbind(node.to.data.frame(getancestors("Blackhawk")), node.to.data.frame(getdescendants("Blackhawk"))))
# qplot(data=test, x=x, y=y, label=label, geom="text", vjust=0, hjust=.5) + 
#   geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend)) + 
#   geom_segment(aes(x=xend, y=yend, xend=branchx, yend=branchy)) +
#   theme_bw()