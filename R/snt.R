###Sackin's Index###
#Originally written by S. Frost
#Adapted for use with heterochronous trees by B. Dearlove

#Define thresh:
#Utility function that returns the length of a vector below a certain threshold.
thresh <- function(v,th){
  length(v[v<=th])
}
 
#Define mtips:
#Utility function returning the number of tips in a tree, if tree is of class 'phylo'
mytips <- function(phy){
  if (!inherits(phy, "phylo")){
    return(NA)
  }
  else{
    n <- length(phy$tip.label)
    return(n) #returns number of tips in tree (if not of type phylo, returns NA)
  }
}

#' snt
#' 
#' Gets the Sackin index over the history of a phylogenetic tree.
#' @param tree an object of class \code{phylo}.
#' @return A list containing the node times, and Sackin's index at those times.
#' @export
#' @examples
#' tree<-rtree(8)
#' sackint<-snt(tree)

#Get Sackin's index over time
snt <- function(phy,tm=NULL,nsteps=NULL){
  if (!inherits(phy, "phylo")){
    return(NA)
  }
  phy.ntips <- mytips(phy)    #Get number of tips
  phy4 <- as(phy,"phylo4")    #Convert to phy4 format (allows root to be found easily)
  phy4.root <- rootNode(phy4) #Get root of tree
  phy4.nodelist <- list()     #Initialise matrix to hold nodes from tip i to root
  for(i in 1:phy.ntips){
    phy4.nodelist[[i]] <- .tipToRoot(phy4, i, phy4.root, include.root=TRUE) #Lists all nodes between root and tip i.
  }
  #phy.bt <- branching.times(phy)  #Get branching times - only works for homochronous trees
  phy.bt<-max(dist.nodes(phy)[phy4.root,1:phy.ntips])-dist.nodes(phy)[phy4.root,(phy.ntips+1):(2*phy.ntips-1)]
  #if(length(phy.bt)!=length(unique(phy.bt))){
  #  phy.bt[which(duplicated(phy.bt))]<-phy.bt[which(duplicated(phy.bt))]+1/10000000
  #}	
  phy4.nodetimelist <- list()     #Initialise node times
  for(i in 1:phy.ntips){
    nodelabels <- as.character(phy4.nodelist[[i]])  #Get node labels for nodes above tip i
    idx <- match(nodelabels,names(phy.bt))          #Get branching times for nodes above tip i
    bt <- as.double(phy.bt[idx])  #Branching times
    phy4.nodetimelist[[i]] <- bt  #save branching times
  }
  tmrca <- max(unlist(phy4.nodetimelist))   #Find tmrca 
  if(is.null(tm)){
    if(!is.null(nsteps)){
      tm <- seq(0,tmrca,by=tmrca/(nsteps-1))
    }
    else{
      tm <- c(0,sort(unique(unlist(phy.bt))))  #Get branching times- either given, or from tree
    }
  }
  nst <- length(tm)   #Time intervals, starting at zero
  snt <- rep(0,nst)   #Initialise Sackin through time - all zero
  for(i in 1:nst){
    lens <- lapply(phy4.nodetimelist,thresh,tm[i])
    snt[i] <- sum(unlist(lens))
  }
  list(tm,snt)
}