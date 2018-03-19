# Marina Golivets
# 11/10/2016
# Function to calculate phylogenetic correlation matrix

library(ape)

phylo_corr <- function (tree) {
  
    tree <- read.tree(tree)
    tree <- root(tree, "Dennstaedtia_punctilobula", resolve.root = TRUE)
    #is.rooted(tree)
    tree <- chronos(tree)
    vcv_matrix <- vcv.phylo(tree)
    return(vcv_matrix)
}


