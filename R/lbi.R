#' LBI
#' 
#' This calculates the local branching index (LBI), defined by Neher, Russell and Shraiman (2014) as the the tree length surrounding a given node or tip exponentially discounted with increasing distance from that node. 
#' @param tree a tree of class \code{phylo}.
#' @param tau a number giving the scaling factor to use, suggested by Neher et al. to be 0.0625 times the average pairwise distance in the sample.
#' @param transform a function by which to transform the LBI. Default is no transformation.
#' @return an vector of class \code{numeric} giving the LBI for the tips and internal nodes in the same order as given in the tree \code{phylo} object.
#' @keywords local branching index, LBI
#' @references Neher RA, Russell CA, Shraiman BI (2014) Predicting evolution from the shape of genealogical trees. eLife 3:e03568. doi:10.7554/eLife.03568
#' @export
#' @examples
#' tree<-rtree(50)
#' tree.lbi<-lbi(tree)
#' tree$node.label<-tree.lbi[(tree$Nnode+2):length(tree.lbi)]
#' plot(tree, show.node.label = TRUE)

lbi<-function(tree, tau=0.0005,transform=function(x){x}){
  nnodes<-tree$Nnode
  ntips<-tree$Nnode+1

  node.times<-dist.nodes(tree)[(ntips+1),1:(2*ntips-1)]
    
  node.postorder<-as.numeric(names(sort(node.times,decreasing=T)))
  node.preorder<-as.numeric(names(sort(node.times)))
  
  #Initialise the arrays
  uppolarizer<-array(0,(2*ntips-1))
  downpolarizer<-array(0,(2*ntips-1)) 
  lbi<-array(0,(2*ntips-1))
  
  #Traverse tree postorder to calculate message to parents
  for(i in node.postorder){
    node.uppolarizer<-0
    
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    node.uppolarizer<-node.uppolarizer+sum(uppolarizer[node.children])

    node.edge<-which(tree$edge[,2]==i)
    if(length(node.edge>0)){
      bl<-tree$edge.length[node.edge]/tau
      node.uppolarizer<-node.uppolarizer*exp(-bl)
      node.uppolarizer<-node.uppolarizer+tau*(1-exp(-bl))
    }
    uppolarizer[i]<-node.uppolarizer
  }
  
  #Traverse tree preorder to calculate message to children
  for(i in node.preorder){
    node.downpolarizer<-downpolarizer[i]
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    
    for(j in node.children){
      downpolarizer[j]<-node.downpolarizer+sum(uppolarizer[setdiff(node.children,j)])
      
      child.edge<-which(tree$edge[,2]==j)
      bl<-tree$edge.length[child.edge]/tau
      
      downpolarizer[j]<-downpolarizer[j]*exp(-bl)+tau*(1-exp(-bl))
    }
  }
  
  #Calculate lbi
  for(i in node.postorder){
    node.lbi<-downpolarizer[i]
    
    node.children<-tree$edge[which(tree$edge[,1]==i),2]
    node.lbi<-node.lbi+sum(uppolarizer[node.children])
    
    lbi[i]<-transform(node.lbi)
  }
  return(lbi)
}

