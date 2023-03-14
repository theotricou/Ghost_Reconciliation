#!/usr/bin/env Rscript
# by Theo Tricou

library(reshape)
library(ggplot2)
d = read.table('RESULTS.txt', h = T)



res <- as.data.frame(do.call("rbind", by(d, list(d$sim, d$tool), function(x) {

  sTP <- sum(x[x$status == "VP", "score"])
  sFP <- sum(x[x$status == "FP", "score"])
  sFN <- sum(x[x$status == "FN", "score"])
  precis <- sTP / (sTP + sFP)
  recall <- sTP / (sTP + sFN)

  return(list(tool = x[1, 4], sim = x[1, 7], precision = precis, recall = recall))

})))

RES = as.data.frame(lapply(res, unlist))

RES.melt <- melt(RES, id.vars = c("sim", "tool"), measure.vars = c("precision", "recall"))

ggplot(RES.melt, aes(x = tool, y = value, color = variable)) + geom_boxplot()



# GNU Ghost
