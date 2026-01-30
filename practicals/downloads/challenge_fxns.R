library(ape)

##determines whether a node is a tip or not
isTip<-function(phy,nd){
  
  ##ape is pretty clever about how they organize trees
  ##you may have noticed that for n nodes, the node numbers of the tips are always numbered 1 to n.
  
  ##This means that in order to check whether a node is a tip we only have to ask whether the number is less than or equal to n
  ntips<-length(tree$tip.label) ##The number of tips
  val<- nd<=ntips ##checks if our node number is smaller or equal to the number of tips
  return(val)
  
  ##For a more convoluted but equally valid way of determining whether a node is a trip you can see the code commented out below.
  
  # ##The two columns of the edge matrix actually have some meaning, the represent the directionality of the edges
  #   ##The first column represents the parent node of the edge
  #   ##the second column represents the child node of the edge
  # 
  # ##To find the tips, all we need to do is find the node numbers that never appear in the first column, they are never parents and thus must be tips
  # tips<-setdiff(tree$edge[,2],tree$edge[,1]) ##Finds nodes only in the second column
  # val<- nd %in% tips ##Checks to see if our node number is amongst the tip numbers
  # return(val)
  
}

##Gets the node numbers of the children of a node
getChildren<-function(phy,nd){
  
  ##to get the children we find all the rows in the edge matrix where our nd appears in the first column
  rows<-which(tree$edge[,1]==nd)
  
  ##next we get all of the values of the second column for those rows, those are the children
  children<-tree$edge[rows,2]
  
  if(length(children)==0){ ##if no children were found we want to return NULL
    children<-NULL
  }
  
  return(children)
}

##gets the branch length of the branch that leads to the given node
getBranchLength<-function(phy,nd){
  
  ##We want to find the row index in the edge matrix where the nd is the child
  rw<-which(tree$edge[,2]==nd)
  
  ##we return that edge length 
  val<-tree$edge.length[rw]
  
  if(length(val)==0){##If there are no branches where the node is the child then we are at the root with no branch. Return NULL
    val<-NULL
  }
  return(val)
}

#computes the probability of a nucleotide transitioning to another nucleotide along a branch given a Juke-Cantor substitution model
subst_probsJC<-function(i,j,v){
  ##There are two cases we need to consider:
    ##The nucleotide is the same at the beginning and at the end. i.e. i=j
    ##The nucleotide is the different at the beginning and the end. i.e. i=/=j
  
  ##The math here comes from the JC69 nucleotide substitution model
  if(i==j){##the nucleotide is the same from start to end of the branch
    val<-0.25+0.75*exp((-4*v)/3)
  }else{ ##The nucleotide changes from the start to the end of the branch
    val<-0.25-0.25*exp((-4*v)/3)
  }
  return(val)
}

##This computes the conditional probability of each nucleotide of the parent node given:
  ##the likelihood of each nucleotide at the two child nodes
  ##the branch lengths that lead to the child nodes
  ##a model for nucleotide substitution
nodeLikelihood<-function(l1,l2,v1,v2,sub_model=subst_probsJC){
  
  L<-c() ##This will store the node likelihood for each nucleotide
  
  nucleotides<-names(l1) ##This gets us the name of each base. In most cases this will be "A","C","G",and "T"
  
  for(i in nucleotides){ ##compute for each base..
    
    ##..the conditional probability of that base for one of the two children
    
    running_sum1<-0
    for(j in nucleotides){##go thru the likelihood of each base for that child
      running_sum1<-running_sum1 + sub_model(i,j,v1)*l1[j]  
    }
    
    ##We then do the same for the other child
    running_sum2<-0
    for(j in nucleotides){
      running_sum2<-running_sum2 + sub_model(i,j,v2)*l2[j]
    }
    
    l_i<-running_sum1*running_sum2 ##We then use the AND rule to compute the conditional probability of the observed states of the children given that the parent was in state i
    L[i]<-l_i ##store the value
  }
  return(L)
}



##Gets the nucleotide for every species at a given site
##Then turns the site sequence into a likelihood form that the recursive function can use
siteSeqs2Likelihood<-function(alignment, site_no){
  
  ##First get the nucleotide for every species at a site
  
  site_seqs<-c() ##This will store the nucleotides for each species at a given site in the allignment
  species<-names(seq) ##get all the names of the species
  for(s in species){ ##go thru each species and get the nucleotide for that site
    
    species_seq<-alignment[[s]] ##Get the sequence for the species we are interested in
    site_seqs[s]<-species_seq[site_no] ##store the nucleotide at the site of interest
  }
  
  ##Now turn those site sequences into a likelihood format
  
  likes<-list() ##we want the nucleotide for each species in a likelihood form. We will store a likelihood form for each species here
  
  ##This will serve as a template for each the likelihood form of each species
  lik_template<-c(0,0,0,0)
  names(lik_template)<-c("a","c","g","t")
  for(sp in names(site_seqs)){ ##For every species..
    ##..make a likelihood form for the sequence data
    
    sp_lik<-lik_template ##copy the template
    sp_lik[site_seqs[[sp]]]<-1 ##put a 1 at the observed nucleotide
  
    likes[[sp]]<-sp_lik ##Store the likelihood formatted sequence data  
  }
  return(likes)
}

##Compute the likelihood for each nucleotide at a given node
getLikelihoodRecursive<-function(nd,phy,tip_lik,sub_model){
  
  ##First we want to check if we are at a tip. This is our Base Case.
  if(isTip(phy,nd)){
    ##If we are at a tip then getting the likelihood is as simple as getting the the values from tip_lik
    lik<-tip_lik[[phy$tip.label[nd]]]
  }
  else{ ##we are at an internal node 
    ##Computing the likelihood for an internal node requires us to use nodeLikelihood()
    ##however in order to use nodeLikelihood we need the likelihood of the descendent nodes
    
    ##first let's find the descendents
    descs<-getChildren(phy,nd)
    
    ##here's where the recursion is, we call getLikelihoodRecursive on each of the descendents and just assume that it works fine.
    lik1<-getLikelihoodRecursive(descs[1],phy,tip_lik,sub_model) ##get the likelihood of the first descendent
    lik2<-getLikelihoodRecursive(descs[2],phy,tip_lik,sub_model) ##get the likelihood of the first descendent
    
    ##Great! so now we have the likelihood at each base for both descendents. Now all we need are the branch lengths
    v1<-getBranchLength(phy,descs[1])
    v2<-getBranchLength(phy,descs[2])
    
    ##Now compute the likelihood as we normally would with nodeLikelihood()
    lik<-nodeLikelihood(lik1,lik2,v1,v2,sub_model)
    
    
  }
  return(lik)
}

##Compute the site likelihood
siteLikelihood<-function(phy,tip_lik,sub_model,stat_freqs){
  
  ##first we want to find the root of the tree
  rt<-length(phy$tip.label)+1 ##for a tree with n tips, the root node number is n+1
  
  rt_lik<-getLikelihoodRecursive(rt,phy,tip_lik,sub_model) ## This is the recursive function. We call it at the root.
  ##We will just assume right now that it works and that it gets us the likelihood for each base at the root
  
  
  ##Now we just need to multiply the likelihoods by their base frequencies and add them together
  site_lik<-rt_lik*stat_freqs ##multiply by stationary frequencies
  site_lik<-sum(site_lik) ##add all values together
  
  return(site_lik)
}



