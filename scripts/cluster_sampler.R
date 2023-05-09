#!/usr/bin/env Rscript
# by Theo Tricou

require(ape)

args <- commandArgs(TRUE)

t <- read.tree("T/ExtantTree.nwk")
dist <- cophenetic.phylo(t)

name <- sample(t$tip.label, 1)
br_len <- dist[name,]
br_len <- br_len[-which(names(br_len) == name)]
proba <- log((1/br_len+1)) / sum(log((1/br_len+1)))
proba2 <- (log(1/br_len+1)**2) / sum(log(1/br_len+1)**2)
picked <- sample(names(br_len), args[1], replace=FALSE, prob=proba2)
to_keep <- c(name, picked)
# colors = rep("black",length(t$tip.label))
# colors[which(t$tip.label %in% to_keep)] <- "red"
# plot(t, tip.color = colors)

write.table(to_keep, file = "cluster_sample", quote = F, col.names = F, row.names = F)
