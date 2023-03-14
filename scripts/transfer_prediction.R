#!/usr/bin/env Rscript
# by Theo Tricou

REC_ROB <- function(sim){
  # read reconciliation results file
  REC_tr <- read.table(paste(sim, "Transfers_rec.txt", sep = "/"))
  colnames(REC_tr) <- c("gene", "tool", "from", "to", "score")
  # read simulated transfers file
  SIM_tr <- read.table(paste(sim, "Transfers.txt", sep = "/"))
  colnames(SIM_tr) <- c("gene", "from", "to")
  # read node mapping file
  MAP_tr <- read.table(paste(sim, "Nodes_mapping.txt", sep = "/"), h=T)
  # read aletree. Used to determine if transfer are predict between sister
  # branches
  tree = ape::read.tree(paste(sim, "aletree", sep = "/"))
  nodes <- c(tree$tip.label, tree$node.label)
  # all inferred transfers are FP until proven otherwise
  REC_tr$status <- "FP"
  REC_tr$paste <- paste(REC_tr$gene, REC_tr$from, REC_tr$to, sep="@")

  # map simulated transfers
  SIM_tr$map_from <- sapply(SIM_tr$from, function(x) MAP_tr[which(MAP_tr$Node==x), "ALE_ID"])
  SIM_tr$map_to <- sapply(SIM_tr$to, function(x) {
    if (MAP_tr[which(MAP_tr$Node==x), "Cor_Recip"] != "None"){
      MAP_tr[which(MAP_tr$Node==x), "ALE_ID"]
    }else{return("None")}
  })

  # check if transfer is predicted between sister branches or sthe same branch
  # these transfers are undetectable and remove from subsequent analysis
  SIM_tr_red <- SIM_tr[SIM_tr$map_to != "None",]
  SIM_tr_red$map_to <- unlist(apply(SIM_tr_red, 1, function(x){
    if (as.character(as.numeric(x[5])) == as.character(as.numeric(x[4]))){
      return("Self")
    }else{
      if (phangorn::Ancestors(tree, which(nodes == as.character(as.numeric(x[4]))), "parent") == phangorn::Ancestors(tree, which(nodes == as.character(as.numeric(x[5]))), "parent")){
        return("Sister")
      }else{
        return(as.character(as.numeric(x[5])))
      }
    }
  }))

  # remove transfers impossible to detecte, between sister branches or predicte
  # from/to the same branch
  SIM_tr <- SIM_tr_red[!SIM_tr_red$map_to %in% c("Self", "Sister"),]

  # creates an event string for a simpler comparaison between data frames
  SIM_tr$paste <- paste(SIM_tr$gene, SIM_tr$map_from, SIM_tr$map_to, sep="@")
  # a table to account for duplicate events after prediction
  SIM_tr_score <- as.data.frame(table(SIM_tr$paste))
  colnames(SIM_tr_score) <- c("paste", "score")
  # all simulated transfers are FN until proven otherwise
  SIM_tr_score$status <- "FN"

  # creates a dataset of transfer events status by reconciliation tools
  RES <-  as.data.frame(do.call("rbind", by(REC_tr, REC_tr$tool, function(x){
    # transfers that matched both datasets were annotated as TP.
    x[which(x$paste %in% SIM_tr_score$paste), "status"] <- "TP"
    # other transfers stay FN
    rest = SIM_tr_score[-which(SIM_tr_score$paste %in% x$paste),]
    # change "tool" to match current one
    if (nrow(rest) > 0){
      rest$tool = x[1, "tool"]
      res = rbind(x[, c("paste", "tool", "score", "status")], rest)
    }else{
      res = x[, c("paste", "tool", "score", "status")]
    }
    return(res)
  })))
  split_gene <- as.data.frame(stringr::str_split_fixed(RES$paste, "@", 3))
  colnames(split_gene) <- c("gene", "from", "to")
  REC_RES <- cbind(split_gene, RES[, c("tool", "score", "status")])
  REC_RES$sim <- sim
  return(REC_RES)
}


simulation_replicates = Sys.glob("sim*")
results <- as.data.frame(do.call("rbind",lapply(simulation_replicates, function(x) REC_ROB(x))))


write.table(results, file = "RESULTS.txt", quote = F, row.names = F)

# GNU Ghost
