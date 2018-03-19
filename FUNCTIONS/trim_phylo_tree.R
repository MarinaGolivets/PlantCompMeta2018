# Marina Golivets
# Last modified: 23/10/2017

# The function trims phylogeny and main data set

trimTree.fn <- function (Tree, Neighbor = TRUE, Target = TRUE) {
  
  if (Neighbor == TRUE & Target == FALSE) {
    
    # neighbor species phylogeny
    sp.exclude.nb <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$neighbor.name.phylo)[which(unique(dat$neighbor.name.phylo) %in% Tree$tip.label)])]
  
    tree.trim.nb <- drop.tip(Tree, sp.exclude.nb, trim.internal = TRUE, subtree = FALSE,
                           root.edge = 0, rooted = is.rooted(Tree))
  
    # subset data (based on species included in the two phylogenies)
    dat <- subset(dat, dat$neighbor.name.phylo %in% Tree$tip.label)
    
    return(list(Data = dat, TreeNeighbor = tree.trim.nb))

  }
  
  if (Neighbor == TRUE & Neighbor == FALSE) {
    
    # target species phylogeny
    sp.exclude.tg <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$target.name.phylo)[which(unique(dat$target.name.phylo) %in% Tree$tip.label)])]
    
    tree.trim.tg <- drop.tip(Tree, sp.exclude.tg, trim.internal = TRUE, subtree = FALSE,
                           root.edge = 0, rooted = is.rooted(Tree))
    
    # subset data (based on species included in the two phylogenies)
    dat <- subset(dat, dat$target.name.phylo %in% Tree$tip.label)
    
    return(list(Data = dat, TreeTarget = tree.trim.tg))

  }
  
  if (Neighbor == TRUE & Target == TRUE) {
    
    # neighbor species phylogeny
    sp.exclude.nb <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$neighbor.name.phylo)[which(unique(dat$neighbor.name.phylo) %in% Tree$tip.label)])]
    
    tree.trim.nb <- drop.tip(Tree, sp.exclude.nb, trim.internal = TRUE, subtree = FALSE,
                             root.edge = 0, rooted = is.rooted(Tree))
    
    # target species phylogeny
    sp.exclude.tg <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$target.name.phylo)[which(unique(dat$target.name.phylo) %in% Tree$tip.label)])]
    
    tree.trim.tg <- drop.tip(Tree, sp.exclude.tg, trim.internal = TRUE, subtree = FALSE,
                             root.edge = 0, rooted = is.rooted(Tree))
    
    dat <- subset(dat, dat$neighbor.name.phylo %in% Tree$tip.label)
    dat <- subset(dat, dat$target.name.phylo %in% Tree$tip.label)
    
    
    # trim phylogenies again
    
    # neighbor species phylogeny
    sp.exclude.nb <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$neighbor.name.phylo)[which(unique(dat$neighbor.name.phylo) %in% Tree$tip.label)])]
    
    tree.trim.nb <- drop.tip(Tree, sp.exclude.nb, trim.internal = TRUE, subtree = FALSE,
                             root.edge = 0, rooted = is.rooted(Tree))
    
    # target species phylogeny
    sp.exclude.tg <- Tree$tip.label[-which(Tree$tip.label %in% unique(dat$target.name.phylo)[which(unique(dat$target.name.phylo) %in% Tree$tip.label)])]
    
    tree.trim.tg <- drop.tip(Tree, sp.exclude.tg, trim.internal = TRUE, subtree = FALSE,
                             root.edge = 0, rooted = is.rooted(Tree))
    
    return(list(Data = dat, TreeNeighbor = tree.trim.nb,  TreeTarget = tree.trim.tg))
    
  }
  
}

#### END OF CODE ####