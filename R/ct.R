#' ct
#' 
#' Gets the number of cherries over the history of a phylogenetic tree.
#' @param tree an object of class \code{phylo} or \code{treeshape}.
#' @return A list containing the node times, and number of cherries at those times.
#' @export
#' @examples
#' tree<-rtree(8)
#' cherries<-ct(tree)
 
#Cherries through time
ct<-function(tree){
  if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
    return("Input needs to be of class phylo or treeshape.")
  }
  else{
    if(inherits(tree,"treeshape")){
      tree<-as.phylo(tree)
    }
    N<-tree$Nnode+1
    n.node<-tree$Nnode
    tree4 <- as(tree,"phylo4")    #Convert to phy4 format (allows root to be found easily)
    tree4.root <- rootNode(tree4) 
    
    if(n.node!=N-1){
      stop("\"phy\" is not fully dichotomous")
    }
    
    if(N<4){ 
      return(NA)
    }
    tree.bt<-max(dist.nodes(tree)[(N+1),1:N])-dist.nodes(tree)[(N+1),(N+1):(2*N-1)]
    tree.times<-tree.bt[order(tree.bt)]
    node.names.ordered<-as.numeric(names(tree.times))
    cherries<-0
    
    for(i in 1:length(node.names.ordered)){
      edges.sofar<-node.names.ordered[1:i]
      cherries[(i+1)]<-sum(tabulate(tree$edge[which(tree$edge[,1]%in%edges.sofar),1][tree$edge[which(tree$edge[,1]%in%edges.sofar), 2]<=N])==2)      
      
    }
    tm<-c(0,tree.times)
    list(tm,cherries)
  }
}