---
title: "Practical: Likelihood of a Tree"
nav_exclude: true
math: katex
---

## Download the following files

1. [Nexus file with data](downloads/challenge_data.nex)
2. [Newick phylogenetic tree](downloads/challenge_tree.tree)
3. [Functions to navigate tree](downloads/challenge_fxns.R)

Open your BIO508-practicals project and deposit the files you just downloaded there. Open a new R script file and name it something you remember (likelihood_for_trees.R for example) 


## Reading And Working With The Data

First, let's read in the data and see how it is represented.

```
##Load in relevant packages
library(ape) ##This package has a lot of tree-related functions and has the 'phylo' class

##Load in prewritten functions
source("challenge_fxns.R")


##read in our files 
tree<-read.tree("challenge_tree.tree")
seq<-read.nexus.data("challenge_data.nex")

```


### Basic information about your tree and data

```


print(tree)##We can see some basic information by just printing our tree

str(tree)##We can see the attributes and class of our tree by using the str() function


```

Describe these four attributes of your tree 

* **edge:** 

* **edge.length:** 

* **Nnode:** 

* **tip.label:**


We can access the individual attributes of the tree by typing our variable followed by *$* and then the name of the attribute. For example:

```
tree$edge
```

Plotting your tree

```
{
  plot(tree)
  edgelabels(tree$edge.length)
}
```


Alternatively, we can plot the node and edge numberings on the tree. These node and edge numberings correspond to what is in the *edge* matrix. 

```
{
  plot(tree,show.tip.label = FALSE)
  edgelabels()
  nodelabels()
  tiplabels()  
}
```

Which color are each of the attributes
* **Edge Numbers** 
* **Internal Node Numbers** 
* **Tip Node Numbers** 

This plotting is also a good place to make better sense of the *edge* matrix attribute.
```
tree$edge
```

We can see that the second row, denoted by $$[2,]$$, has $$5$$ and $$1$$ in the columns. We can look at our plot and find the branch labeled with a $$2$$ and see that it is connected to nodes labeled $$5$$ and $$1$$. 

We can also notice that the ordering the columns matters. Node numbers that appear in the first column of the edge matrix are considered the ancestral node of the branch while node numbers in the second column are considered the descendant node of the branch. We can easily see that $$5$$ is the ancestor, or parent, node for edge $$2$$ while node $$1$$ is the descendent, or child, of the branch. Pretty slick!

### The Sequence Data

Similar things can be done for the sequence data to get a feel for how it is represented.
```
str(seq)
```
Here we can see that we have a list where each attribute in the list corresponds to one of our tip names. Then, within each list is our sequence information. Usually each attribute contains a vector of characters, however, we only read in an allignment with one character. 

## Preliminary functions

Before we get into the inner workings of the pruning algorithm, we will need a few functions to make our lives easier. Examples of these functions will refer to the node and edge numbers that were plotted above, so either be ready to scroll up to that tree or it may be handy to jot the tree down on a piece of scratch paper. Lastly, for brevity's sake I won't bother showing the code for the functions here. The interested reader can find the code along with explainations of the code in the *challenge_fxns.R*.

* `isTip(phy,nd)`: this function returns `TRUE` if the given node is a tip and `FALSE` otherwise.This function has two inputs:
    + $$phy$$: the phylogenetic tree  
    + $$nd$$: the number of the node that we are interested in

```
###Example###
isTip(tree,2) ##Returns?
isTip(tree,4) ##Returns?
```
   
* `getChildren(phy,nd)`: This function finds the child node numbers of a given node on a tree. This function has two inputs:
    + $$phy$$: the phylogenetic tree  
    + $$nd$$: the number of the node that we are interested in
```
###Example###
getChildren(tree,4) ##Returns ?
getChildren(tree,2) ##Returns ?
```    


* `getBranchLength(phy,nd)`: This function returns the length of the branch between a given node and its parent. This function has two inputs:
    + $$phy$$: the phylogenetic tree  
    + $$nd$$: the number of the node that we are interested in
```
###Example###
getBranchLength(tree,3) ##Returns ?
getBranchLength(tree,5) ##Returns ?
getBranchLength(tree,4) ##The root has no branch that leads to it. Returns NULL
```  

* `subst_probsJC(i,j,v)`: This function computes the probability that a site transitions from nucleotide $$i$$ to nucleotide $$j$$ along a branch with length $$\nu$$ when the rate of evolution is $$\lambda=1/3$$ The probabilities are calculated according to the Jukes Cantor model of sequence evolution:
$$
p_{ij}(v) = \left\{
        \begin{array}{ll}
            \frac{1}{4} + \frac{3}{4}e^{\frac{-4\nu}{3}} & \quad i = j \\
            \frac{1}{4} - \frac{1}{4}e^{\frac{-4\nu}{3}} & \quad i \neq j
        \end{array}
    \right.
$$

This function has three inputs:
    + $$i$$: A string containing the starting nucleotide. Either `"a"`, `"c"`, `"g"`, or `"t"`
    + $$j$$: A string containing the ending nucleotide. Either `"a"`, `"c"`, `"g"`, or `"t"`
    + $$\nu$$: a non-negative floating point number that represents the expected number of substitutions along a branch.
    
```
###Example###
subst_probsJC("a","g",1) ##Returns what value?
subst_probsJC("a","a",1) ##Returns  what value

##Compute the probability of transitioning from an A to any other base 
A_A<-subst_probsJC("a","a",0.1)  #A to A
A_C<-subst_probsJC("a","c",0.1)  #A to C
A_G<-subst_probsJC("a","g",0.1)  #A to G
A_T<-subst_probsJC("a","t",0.1)  #A to T
A_A+A_C+A_G+A_T ##Returns 1. This is a check for correctness. The probability of transitioning from A to any other base should be 1.
```  

* `siteSeqs2Likelihood(alignment, site_no)`: This function takes the observed nucleotide data at the tips for a given site and converts it into a likelihood format that our other functions can use to compute likelihoods for internal nodes. This function has two inputs:
    + $$alignment$$: a list of aligned sequences  
    + $$site\_no$$: the position of the alignment that we want to put into a likelihood format 

```

###Example###
lik_seq<-siteSeqs2Likelihood(seq,1) ##put the first site of our alignment into the correct format

##We made a list where each tip name is an attribute that contains a vector with a 1 for the observed base and 0s everywhere else
lik_seq 


```

* `nodeLikelihood(l1,l2,v1,v2,sub_model=subst_probsJC)`: This function computes the conditional probability of each nucleotide for the ancestral node given the likelihood of the nucleotides at the descendent nodes and given the branch lengths that lead to each ancestral node. The likelihood of the ancestor node, denoted $$anc$$, for a given base $$i$$ is computed as follows:
$$
\ell_{anc}(i) = \left(\sum_{j} p_{ij} (\nu_{d1})\ell_{d1}(j) \right)*\left(\sum_{j} p_{ij} (\nu_{d2})\ell_{d2}(j) \right)
$$
 Where $$d1$$ and $$d2$$ represent the two descendent nodes of the ancestral node. 
 The computation is repeated for each of the four nucleotides. The function has five inputs:
    + $$l_{1}$$, $$l_{2}$$ These inputs contain the likelihood of each base for the two descendent nodes
    + $$\nu_{1}$$, $$\nu_{2}$$: The expected number of substitutions along the two branches that lead to the descendent nodes
    + $$sub\_model$$: The model for sequence evolution. The default value is the Jukes Cantor model of sequence evolution, `subst_probsJC`

```
###Example###

##compute the likelihood of each base for the node numbered 5 in our tree

##We can see that the nodes numbered 1 and 2 are the descendents of 5, so we need the sequence data in a likelihood form
lik_seq<-siteSeqs2Likelihood(seq,1)

##the nucleotide likelihood for each child
node1_lik<-lik_seq$t1
node2_lik<-lik_seq$t2

##the branch lengths for each child
v1<-getBranchLength(tree,1) ##branch length for the branch between nodes 5 and 1
v2<-getBranchLength(tree,2) ##branch length for the branch between nodes 5 and 2

node5_lik<-nodeLikelihood(node1_lik,node2_lik,v1,v2,sub_model=subst_probsJC) ##compute the likelihood of each base for node 5
node5_lik
##These numbers look like what we saw in lecture! Nice!


```

## The Pruning Algorithm

### Doing the Algorithm Manually

We should now have all the infrastructure we need to compute the site likelihood. All we need to do for this is:
  
  1. Compute the likelihood for each base at the root
  2. Multiply the likelihood for each base by their respective stationary frequency
  3. Add all the elements together.
  
For step 1. we can use the `nodeLikelihood()` function to compute the likelihood for each base at the root, the node numbered $$4$$. However, in order to do this we need the likelihood at each of the descendent nodes, the nodes numbered $$3$$ and $$5$$. 

The likelihood for node $$3$$ is straightforward, it is a tip and we have observed data for that node. As such, we can just convert the sequence data into a likelihood format for that node. 

```
node3_lik<-lik_seq$t3 ##get the likelihood for node 3
##We stored
## a c g t 
## 1 0 0 0

v3<-getBranchLength(tree,3) ##also get the branch length that connects nodes 4 and 3


```


Node $$5$$ is a bit more tricky since this isn't a tip. In order to get the likelihood for the bases at node $$5$$ we need to use the `nodeLikelihood()` function again but this time on node $$5$$ where nodes $$1$$ and $$2$$ are the descendents. Luckily, we already did this computation in the example for `nodeLikelihood()`.


```

###See the example for nodeLikelihood to see the likelihood computation for node 5

##The likelihood is stored in a variable called node5_lik
v5<-getBranchLength(tree,5) ##also get the branch length that connects nodes 4 and 5

v5
```

We only need to plug in values to calculate the base likelihoods at the root
```
root_lik<-nodeLikelihood(node3_lik,node5_lik,v3,v5,subst_probsJC) ##compute the likelihood at the root for each base
root_lik #These numbers look like what we saw in class!

```

So we've calculated the likelihood for each base, now for step 2. We need to multiply these by their stationary frequencies. The Jukes Cantor model assumes stationary frequencies of $$A=C=G=T=0.25$$. After we do this multiplication we only need to sum the values together to get the site likelihood

```
stationary_freqs<-c(0.25,0.25,0.25,0.25) ##The stationary frequencies are all 0.25
site_lik<-root_lik*stationary_freqs ##multiply by stationary frequencies

site_lik<-sum(site_lik) ##add all values together
site_lik  #what a beaut


```
Although we did the calculations by hand, hopefully this example gives somewhat of an idea for how this process can be applied more generally. Because doing these calculations manually quickly becomes cumbersome as the number of taxa and length of sequences grow. We will want to write a recursive algorithm to automate the process of going through the nodes to compute likelihoods and slowly work our way to the root.


## Fitting a phylogeny with real data

We will use a common example dataset of 12 species of primates plus two outgroups (cow and mouse). This alignment was assembled by Dr. Masami Hasegawa of the Institute of Statistical Mathematics in Tokyo, from sequencing done by Kenji Hayasaka, Takashi Gojobori, and Satoshi Horai (Molecular Biology and Evolution 5: 626-644, 1988). For more details, see [J. Felsenstein's website](http://evolution.gs.washington.edu/gs570/2016/#data).

Download the data file using R:

```
download.file("https://raw.githubusercontent.com/EEOB-Macroevolution/EEOB565X-Spring2020/master/practicals/01-intro-to-phylo/primates.dna","primates.dna")
```

This will download the file `primates.dna` into your current working directory. Now we can read in the file using the function `read.phyDat()`.

```
#install.packages("phangorn") # Remember to remove the hashtag for installation but once installed you can put it back in this line because it will be already in RStudio
library("phangorn")
primates <- read.phyDat("primates.dna", format = "interleaved")
primates
```

### Parsimony 

To begin, we will first estimate a tree using parsimony. To perform a heuristic search under parsimony, we can start with a tree built using a random-addition algorithm.

```
primates_start_tree <- random.addition(primates)
```

First let's check the parsimony score of our starting tree using the `parsimony()` function:

```
parsimony(primates_start_tree, primates)
```


Now let's plot the tree.

```
plot(primates_start_tree)
```

You will notice that the tree does not have meaningful branch lengths. We can represent the branch lengths in units of the number of changes (under parsimony) using the `acctran()` function. 

```
primates_start_tree <- acctran(primates_start_tree, primates)
```

The tree should now be depicted as a phylogram where the branch lengths correspond to the number of state changes along a branch.

```
plot(primates_start_tree)
```

Now we can optimize the tree topology under parsimony using the `optim.parsimony()` function. This performs a series of tree rearrangements under an algorithm called the nearest-neighbor-joining. 

```
primates_opt <- optim.parsimony(primates_start_tree, primates)
primates_opt <- acctran(primates_opt, primates)
```

What is the parsimony score of the optimized tree?

```
parsimony(primates_opt, primates)
```


Plot the optimized tree:

```
plot(primates_opt)
```

What is different between the starting tree and the optimized tree?


### Maximum Likelihood

Compute the likelihood of the parsimony tree and parsimony branch lengths. For this we use the `pml()` function of _phangorn_. By default, this approach computes the likelihood under the Jukes-Cantor (1969) model:

```
fitJC <- pml(primates_opt, data=primates)
```

The object `fitJC` contains the parameter estimates and likelihood of the JC model given the data and the tree topology and branch lengths we estimated under parsimony. We can view these values by calling the variable:

```
fitJC
```

We do not want to rely on parsimony for the tree topology and branch lengths, but this tree provides a reasonable starting place. We can use the `optim.pml()` function to optimize the tree and model parameters under JC69. This will perform NNI tree rearrangements to identify the topology that maximizes the likelihood.

```
fitJC <- optim.pml(fitJC, optNni=TRUE)
plot(fitJC$tree)
```

If we want to optimize the tree using a different model like GTR, we can use the `pml` object as an argument and specify the `GTR` model.

```
fitGTR <- optim.pml(fitJC, model="GTR", optNni=TRUE)
plot(fitGTR$tree)
```

We can try to reroot the tree using the two outgroup taxa. However, this will return an error if the outgroup is not monophyletic.

```
GTR_tree_rooted <- root(fitGTR$tree,c("Bovine","Mouse"))
```

View the details of the `fitGTR` object. What is the log-likelihood?

```
fitGTR
```

How are the two trees different? What are possible reasons for this ?

One thing that may be bothering you is how do you decide if you should analyze your data under JC69 or under GTR. One option is to compare the models statistically. You can use a likelihood ratio test to compare the two models.

```
anova(fitJC, fitGTR)
```


The results of this comparison show that GTR (model 2) is best supported. 

