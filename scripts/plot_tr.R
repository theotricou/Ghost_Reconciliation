#!/usr/bin/env Rscript
# by Theo Tricou

library(ape)
library(tidyr)
library(phangorn)
trcol <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))}

tree <- read.tree('T/CompleteTree.nwk')

tree_samp <- read.tree('Sample_Cluster/SampledSpeciesTree.nwk')
spnd<-c(tree$tip.label, tree$node.label)

ll = unique(unlist(sapply(tree_samp$tip.label, function(x){
  node = c(which(spnd == x), Ancestors(tree, which(spnd == x), type='all'))
  which(tree$edge[, 2] %in% node)
})))


which(spnd %in% tr$to) %in% which(list_color == "grey")



list_color = rep("lightgrey", length(spnd))
list_color[ll] <- "red"
tr <- read.table("G/Gene_families/58_events.tsv", h = T)
root_L = tree$root.edge
end = tail(tr[, 1], 1)
tr = tr[tr$EVENT == "T", ]
edge<-tree$edge
plot(tree, direction = "downwards", no.margin = TRUE,
  show.tip.label = F, show.node.label = T, root.edge = T)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
x<-lastPP$xx
y<-lastPP$yy
tr = separate(data = tr, col = "NODES", into = c("from", "from_up","from_rep","from_down","to","to_down"), sep = ";")
par(new=TRUE)
plot(tree, direction = "downwards", no.margin = TRUE,
  show.tip.label = F, show.node.label = T, root.edge = T,
  edge.col=list_color, edge.width = 4)
mat <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
ll_to_ghost <- which(which(tree$edge[,2] %in% which(spnd %in% tr$to)) %in% ll)
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,7]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  if ( from == tree$node.label[1] ) {
    fromdad = wherefrom
  }else if ( to == tree$node.label[1] ) {
    todad = whereto
  }
  ywherefrom <- end -(tr[i,1])
  ywhereto <- end - (tr[i,1])
  col_tr = "orange"
  if (which(tree$edge[,2] == which(spnd == tr[i,7])) %in% ll){
    col_tr = "blue"
  }
  mat[i,] <- c(x[wherefrom], ywherefrom, x[whereto], ywhereto, col_tr)
}
arrows(as.numeric(mat[s,1]),as.numeric(mat[s,2]),as.numeric(mat[s,3]),as.numeric(mat[s,4]), col = mat[s,5],lwd=1.3, length = 0.1)
tiplabels(tree$tip.label, adj = c(0, 0.5), srt=-90, bg = "white", frame = "none")
