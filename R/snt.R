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
snt <- function(phy){
  if (!inherits(phy, "phylo")){
    return(NA)
  }
  phy.ntips <-phy$Nnode+1 #Get number of tips
  phy4 <- as(phy,"phylo4")    #Convert to phy4 format (allows root to be found easily)
  phy4.root <- rootNode(phy4) #Get root of tree
  phy.bt<-max(dist.nodes(phy)[phy4.root,1:phy.ntips])-dist.nodes(phy)[phy4.root,(phy.ntips+1):(2*phy.ntips-1)]
  tm<-c(0,sort(phy.bt))
  phy.nodes<-as.numeric(names(sort(phy.bt)))
  phy.snt<-array(0,(phy$Nnode+1))
  phy.delta<-array(0,(phy$Nnode+1))
  for(i in 1:(phy$Nnode)){
    phy.delta[i+1]<-length(extract.clade(phy,phy.nodes[i])$tip.label)
    phy.snt[i+1]<-phy.delta[i+1]+phy.snt[i]
  }
  list(tm,phy.snt,phy.delta)
}

