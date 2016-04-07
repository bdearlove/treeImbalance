#' matchSimNode
#' 
#' A function that takes a node number in an observed tree, and finds the matching node in a tree simulated with \code{getSimTree}.
#' @param obstree an object of class \code{phylo}
#' @param simtree an object of class \code{phylo}
#' @param node a number
#' @return the number of the node of height nheight in the simulated tree
#' @export
#' @examples
#' tree<-rtree(8)
#' tree2<-getSimTree(tree)

matchSimNode<-function(obstree,simtree,node){
	n<-tree$Nnode+1	#number of taxa at tips
	
	#Order node heights in observed tree 
	nodeheights.ordered<-order(sapply((ntips+1):(2*ntips-1),function(x) nodeheight(tree,x)))
	
	#Where observed node height falls in order
	obsnode.location<-which(nodeheights.ordered==(node-n))

	#Match in simulated tree using that order
	simnode.location<-simtree$node.label[cluster.location]+ntips
}
