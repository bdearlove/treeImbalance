#' getTipTimes
#' 
#' Gets the times of the tips of a heterochronously sampled tree with most recently sampled tip set to 0.
#' @param tree A tree of class \code{phylo}.
#' @return A object class \code{numeric} giving the tree tip sampling times.
#' @export

#Get tip times (most recently sampled tip at time 0)
#In order of tree$tip.label
getTipTimes<-function(tree){
  n<-tree$Nnode+1
  times<-max(dist.nodes(tree)[(n+1),1:n])-dist.nodes(tree)[(n+1),1:n]
  names(times)<-tree$tip.label
  return(times)
}

#' getCoalescenceEvents
#' 
#' Gets the times of coalescence events (internal nodes), backwards in time from most recently sampled tip.
#' @param tree A tree of class \code{phylo}.
#' @return A object class \code{numeric} giving the coalescence times.
#' @export


#Get coalescence times, backwards in time from most recently sampled tip.
getCoalescenceEvents<-function(tree){
  n<-tree$Nnode+1
  times<-max(dist.nodes(tree)[(n+1),1:n])-dist.nodes(tree)[(n+1),(n+1):(2*n-1)]
  times<-times[order(times)] #Reorder so earliest event listed first
  names(times)<-1:(n-1)
  return(times)
}

#' getSimTree
#' 
#' A function to permute a tree with the same tip and internal node times.
#' @param tree an object of class \code{phylo}
#' @return a tree of class \code{phylo} with the same tip and internal node times as the input tree.
#' @export
#' @examples
#' tree<-rtree(8)
#' tree2<-getSimTree(tree)

#Simulate new topology using same times for sampled tips and coalescence events
getSimTree<-function(tree){
  n<-tree$Nnode+1	#number of taxa at tips
  node.remaining<-1:(2*n-1)	#Vector containing an index for each node
  int.node.times<-dist.nodes(tree)[(n+1),(n+2):(2*n-1)]		#Get internal node times
  int.node.times.ordered<-int.node.times[order(int.node.times)] #Order internal nodes (loop will go from earliest to latest)
  node.times<-c(dist.nodes(tree)[(n+1),1:(n+1)],int.node.times.ordered) 	#Get time from nodes to MRCA (MRCA given at node n+1)
  names(node.times)<-1:(2*n-1)	#Rename nodes where they have been re-ordered (easier to keep track in tree later)
  ancestor<-array(,c(2*n-1))	#Initialise vector to store ancestor information
  new.tree.nodes<-array(,c(n-1,2))	#Initialise array to contain daughter nodes for each internal node
  row.names(new.tree.nodes)<-(2*n-1):(n+1)
  
  #Simulate which nodes coalesce at each event, starting at earliest ancestral node backwards in time(only nodes below current one are viable)
  for(i in (2*n-1):(n+1)){
    this.node<-i	#Define ancestral node
    this.node.time<-node.times[i]		#Get time for current node
    if(this.node>(n+1)){
      valid.nodes.coal<-setdiff(node.remaining[which(node.times[c(node.remaining)]>=this.node.time)],this.node:(n+1))  #Only allow coalescence of nodes below this one
    }else{
      valid.nodes.coal<-setdiff(node.remaining[which(node.times[c(node.remaining)]>=this.node.time)],this.node)  #Only allow coalescence of nodes below this one
    }
    nodes.coal<-sample(valid.nodes.coal,2,replace=F)	#Sample which two valid nodes should coalesce
    ancestor[nodes.coal]<-i		#Set current node to be ancestor of nodes just sampled
    new.tree.nodes[(2*n-i),]<-nodes.coal	#Save daughter nodes of this node
    node.remaining<-setdiff(node.remaining,nodes.coal)	#Remove nodes that have coalesced from those available.
    
  }
  
  #Get branch lengths for new tree
  branches<-node.times-node.times[ancestor]
  
  #Get new tree, based on original	
  tree.sim<-tree
  tree.sim$edge<-cbind(ancestor,1:(2*n-1))[-(n+1),]
  tree.sim$edge.length<-branches[-(n+1)]
  
  #This fixes plotting issues of overlapping branches
  tree.sim<-read.tree(text=write.tree(tree.sim))
  
  #Add node labelling so can compare simulated trees to observed easily 
  #Fix to problem trying to match by node height which fails due to branches calculation above
  #1 = root, n-1 = node closest to present
  tree.sim$node.label<-order(sapply((n+1):(2*n-1),function(x) nodeheight(tree.sim,x)))

  return(tree.sim)
}