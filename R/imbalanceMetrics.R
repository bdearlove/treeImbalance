#' I1
#' 
#' This calculates I1, a weighted average of the balance of internal nodes, where \eqn{N} is the number of tips, \eqn{j} is the number of nodes, and \eqn{r_j} and \eqn{l_j} represent the number of tips in the left and right subtrees respectively. Then,
#' \deqn{ I_1 = \frac{2}{(N-1)(N-2)} \sum_{j \in \mathcal{I} } {|r_j-l_j|}.}
#' I1 is closely related to the Colless index, which can be found using \code{\link[apTreeshape]{colless}}, or by multiplying I1 by \deqn{\frac{(N-1)(N-2)}{2}.}
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry, Colless Index
#' @references
#' \itemize{
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Purvis A, Agapow PM (2002) Phylogeny imbalance: Taxonomic level matters. Systematic Biology 51: 844-854.
#'  \item Fusco G, Cronk Q (1995) A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology 175: 235-243. doi: 10.1006/jtbi.1995.0136
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' I1(tree)

I1<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"phylo")){
			tree<-as.treeshape(tree)
		}

		N<-dim(tree$merge)[1]+1
		rl<-abs(smaller.clade.spectrum(tree)[,1]-2*smaller.clade.spectrum(tree)[,2])
		result<-2*sum(rl)/((N-1)*(N-2))
		return(result)
	}
}

#' I2
#' 
#' This calculates I2, a weighted average of the balance of internal nodes, where \eqn{N} is the number of tips, \eqn{j} is the number of nodes, and \eqn{r_j} and \eqn{l_j} represent the number of tips in the left and right subtrees respectively. Then,
#' \deqn{ I_2 = \frac{1}{(N-2)} \sum_{j \in \mathcal{I},r_j+l_j>2 } {|r_j-l_j|}/{|r_j+l_j-2|}. }
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references 
#' \itemize{
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Purvis A, Agapow PM (2002) Phylogeny imbalance: Taxonomic level matters. Systematic Biology 51: 844-854.
#'  \item Fusco G, Cronk Q (1995) A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology 175: 235-243. doi: 10.1006/jtbi.1995.0136
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' I2(tree)

I2<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"phylo")){
			tree<-as.treeshape(tree)
		}
	
		N<-dim(tree$merge)[1]+1
		rl.num<-abs(smaller.clade.spectrum(tree)[,1]-2*smaller.clade.spectrum(tree)[,2])
		rl.denom<-abs(smaller.clade.spectrum(tree)[,1]-2)
		rl.sum<-smaller.clade.spectrum(tree)[,1]
		result<-sum(rl.num[which(rl.sum>2)]/rl.denom[which(rl.sum>2)])/(N-2)
		return(result)
	}
}

#' Ic
#' 
#' This calculates Ic (as defined in Pompei et al., 2012) for a phylogenetic tree. This asymmetry metric uses a weighted average of the balance of internal nodes:
#' \deqn{ I_c = \frac{1}{(N-1)} \sum_{j \in \mathcal{I}} w_j \frac{\max(r_j,l_j)-m_j}{r_j+l_j-m_j-1}. }
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references
#' \itemize{ 
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Purvis A, Agapow PM (2002) Phylogeny imbalance: Taxonomic level matters. Systematic Biology 51: 844?854.
#'  \item Fusco G, Cronk Q (1995) A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology 175: 235?243. doi: 10.1006/jtbi.1995.0136
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' Ic(tree)


Ic<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"phylo")){
			tree<-as.treeshape(tree)
		}
	
		N<-dim(tree$merge)[1]+1
		w<-array(,c(N-1))
		m<-array(,c(N-1))

		rj<-smaller.clade.spectrum(tree)[,1]-smaller.clade.spectrum(tree)[,2]
		lj<-smaller.clade.spectrum(tree)[,2]
		rl.max<-apply(cbind(rj,lj),1,max)
		m[which(abs(rj+lj)%%2==0)]<-(rj[which(abs(rj+lj)%%2==0)]+lj[which(abs(rj+lj)%%2==0)])/2
		m[which(abs(rj+lj)%%2==1)]<-(rj[which(abs(rj+lj)%%2==1)]+lj[which(abs(rj+lj)%%2==1)]+1)/2
		w[which(abs(rj+lj)%%2==0)]<-1
		w[which(abs(rj+lj)%%2==1)]<-(rj[which(abs(rj-lj)%%2==1)]+lj[which(abs(rj-lj)%%2==1)]-1)/(rj[which(abs(rj-lj)%%2==1)]+lj[which(abs(rj-lj)%%2==1)])
		rl.sum<-smaller.clade.spectrum(tree)[,1]
		result<-sum(w[which(rl.sum>3)]*(rl.max[which(rl.sum>3)]-m[which(rl.sum>3)])/(rj[which(rl.sum>3)]+lj[which(rl.sum>3)]-m[which(rl.sum>3)]-1))/(N-1)
		return(result)
	}
}

#' mean.Iprime
#' 
#' This calculates the mean of I' (Purvis and Agapow, 2012) for a phylogenetic tree. This asymmetry metric uses a weighted average of the balance of internal nodes:
#' \deqn{ mean.Iprime = \frac{1}{(N-1)} \sum_{j \in \mathcal{I}} w_j \frac{\max(r_j,l_j)-m_j}{r_j+l_j-m_j-1}, } where \eqn{w_j = \frac{r_j+l_j-1}{r_j+l_j}} if \eqn{n} is even and \eqn{w_j = 1} if \eqn{n} is odd.
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references
#' \itemize{ 
#'  \item Purvis A, Katzourakis A, Agapow PM (2002) Evaluating phylogenetic tree shape: Two modifications to Fusco and Cronk's method. Journal of Theoretical Biology 214:99-103.
#'  \item Fusco G, Cronk Q (1995) A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology 175: 235?243. doi: 10.1006/jtbi.1995.0136
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' mean.Iprime(tree)


meanIprime<-function(tree){
  if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
    return("Input needs to be of class phylo or treeshape.")
  }
  else{
    if(inherits(tree,"phylo")){
      tree<-as.treeshape(tree)
    }
    
    N<-dim(tree$merge)[1]+1
    w<-array(,c(N-1))
    m<-array(,c(N-1))
    
    rj<-smaller.clade.spectrum(tree)[,1]-smaller.clade.spectrum(tree)[,2]
    lj<-smaller.clade.spectrum(tree)[,2]
    rl.sum<-smaller.clade.spectrum(tree)[,1]
    rl.max<-apply(cbind(rj,lj),1,max)
    m[which(rl.sum%%2==0)]<-(rj[which(rl.sum%%2==0)]+lj[which(rl.sum%%2==0)])/2
    m[which(rl.sum%%2==1)]<-(rj[which(rl.sum%%2==1)]+lj[which(rl.sum%%2==1)]+1)/2
    w[which(rl.sum%%2==0)]<-(rl.sum[which(rl.sum%%2==0)]-1)/rl.sum[which(rl.sum%%2==0)]
    w[which(rl.sum%%2==1)]<-1
    result<-mean(w[which(rl.sum>3)]*(rl.max[which(rl.sum>3)]-m[which(rl.sum>3)])/(rj[which(rl.sum>3)]+lj[which(rl.sum>3)]-m[which(rl.sum>3)]-1))
    return(result)
  }
}

#' M
#' 
#' This calculates M (as defined in Pompei et al., 2012), a measure of asymmetry based on the topological distance, \eqn{M_i}, between tip \eqn{i} and the root:
#' \deqn{ M = \frac{1}{N} \sum_{i \in \mathcal{L}} {M_i}.}
#' This is equivalent to Sackin's index divided by the number of tips (also known as the normalised Sackin index). See also \code{\link[apTreeshape]{sackin}}
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry, Sackin's index
#' @references 
#' \itemize{ 
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Sackin MJ (1972) Good and bad phenograms. Systematic Biology 21: 225-226.
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' M(tree)

M<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"phylo")){
			tree<-as.treeshape(tree)
		}

		N<-dim(tree$merge)[1]+1
		clades = smaller.clade.spectrum(tree)
		result<-sum(clades[,1]/N)
		return(result)
	}
}

#' varM
#' 
#' This calculates \eqn{\sigma^2_M}, a measure of asymmetry based on the topological distance, $M_i$, between tip $i$ and the root:
#' \deqn{ \sigma^2_M = \frac{1}{N} \sum_{i \in \mathcal{L}} {(M_i-M)^2} .}
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references
#' \itemize{ 
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Sackin MJ (1972) Good and bad phenograms. Systematic Biology 21: 225-226.
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' varM(tree)

varM<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"treeshape")){
			tree<-as.phylo(tree)
		}
	
		N<-tree$Nnode+1
		tree<-as(tree, "phylo4")
		root <- rootNode(tree)
		Mi<-array(,c(N))
		for(i in 1:N){
			Mi[i]<-length(.tipToRoot(tree, i, root,include.root=T))
		}
		M.tree<-sum(Mi)/N
		result<-sum((Mi-M.tree)^2)/N
		return(result)
	}
}

#' B1
#' 
#' This calculates B1 for a phylogenetic tree. For each internal node \eqn{j}, we consider the leaf with the largest topological distance \eqn{Z_j}. Then,
#' \deqn{ B_1 = \sum_{j \in \mathcal{I}-root} {Z_j^{-1} }.}
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references 
#' \itemize{ 
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Shao K, Sokal R (1990) Tree balance. Systematic Zooology 39: 266-276
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' B1(tree)

B1<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"treeshape")){
			tree<-as.phylo(tree)
		}
	
		subtree<-subtrees(tree, wait=FALSE)
		Zj<-array(,c(length(subtree)-1))
		
		for(j in 1:length(subtree)){
			N<-subtree[[j]]$Nnode
			tmp.tree<-as(subtree[[j]], "phylo4")
			root<-rootNode(tmp.tree)
			Mi<-array(,c(N))
		
			for(i in 1:N){
				Mi[i]<-length(.tipToRoot(tmp.tree, i, root,include.root=T))
			}

			Zj[j]<-max(Mi)
		}
		
		result<-sum(1/Zj[-1])
		return(result)
	}	
}

#' B2
#' 
#' This calculates B2 for a phylogenetic tree. If we consider the topological distance, \eqn{M_i}, between each leaf \eqn{i} and the root, then
#' \deqn{ B_2 = \sum_{i \in \mathcal{L}} {\frac{M_i}{2^{M_i}} }.}
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references 
#' \itemize{ 
#'  \item Pompei S, Loreto V, Tria F (2012) Phylogenetic Properties of RNA Viruses. PLoS ONE 7(9): e44849. doi:10.1371/journal.pone.0044849
#'  \item Shao K, Sokal R (1990) Tree balance. Systematic Zooology 39: 266-276
#' }
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' B2(tree)

B2<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"treeshape")){
			tree<-as.phylo(tree)
		}
		N<-tree$Nnode+1
		tree<-as(tree, "phylo4")
		root <- rootNode(tree)
		Mi<-array(,c(N))
		for(i in 1:N){
			Mi[i]<-length(.tipToRoot(tree, i, root,include.root=T))
		}
		result<-sum(Mi/(2^Mi))
		return(result)
	}
}

#' Ncherries
#' 
#' This calculates the number of cherries in a phylogenetic tree, where a cherry represents two tips that share a direct ancestor.
#' @param tree A tree of class \code{phylo} or \code{treeshape}.
#' @return An object of class \code{numeric}.
#' @keywords imbalance, asymmetry
#' @references 
#' McKenzie A, Steel M (2000) Distributions of cherries for two models of trees. Math Biosci 164: 81-92.
#' @export
#' @examples
#' N=30
#' tree<-rtreeshape(1,tip.number=N,model="pda")[[1]]
#' Ncherries(tree)

Ncherries<-function(tree){
	if (!(inherits(tree, "phylo")||inherits(tree,"treeshape"))){
		return("Input needs to be of class phylo or treeshape.")
	}
	else{
		if(inherits(tree,"treeshape")){
			tree<-as.phylo(tree)
		}
    		N<-tree$Nnode+1
		n.node<-tree$Nnode
    		
		if(n.node!=N-1){
        		stop("\"phy\" is not fully dichotomous")
		}
    		
		if(N<4){ 
        		return(NA)
		}
    		
		result <- sum(tabulate(tree$edge[, 1][tree$edge[, 2]<=N])==2)
    		return(result)
	}
}