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
snt <- function(tree){
  if (!inherits(tree, "phylo")){
    return(NA)
  }
  tree.ntips <-tree$Nnode+1 #Get number of tips
  tree4 <- as(tree,"phylo4")    #Convert to phy4 format (allows root to be found easily)
  tree4.root <- rootNode(tree4) #Get root of tree
  tree.bt<-max(dist.nodes(tree)[tree4.root,1:tree.ntips])-dist.nodes(tree)[tree4.root,(tree.ntips+1):(2*tree.ntips-1)]
  tm<-c(0,sort(tree.bt))
  tree.nodes<-as.numeric(names(sort(tree.bt)))
  tree.snt<-array(0,(tree$Nnode+1))
  tree.delta<-array(0,(tree$Nnode+1))
  for(i in 1:(tree$Nnode)){
    tree.delta[i+1]<-length(extract.clade(tree,tree.nodes[i])$tip.label)
    tree.snt[i+1]<-tree.delta[i+1]+tree.snt[i]
  }
  list(tm,tree.snt,tree.delta)
}

