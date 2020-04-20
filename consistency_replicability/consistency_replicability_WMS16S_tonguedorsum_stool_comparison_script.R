library(phyloseq)
source("../eval_functions.R")
load(file = "../data/WMS16Ssubsets_replicability_tonguedorsum_stool.RData")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

comparison <- tonguedorsum_stoolWMS16S$tonguedorsum_stoolWMS[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/WMSsubsets_replicability_DA_tonguedorsum_stool",task_id+100,".RDS"))

comparison <- tonguedorsum_stoolWMS16S$tonguedorsum_stool16S[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/16Ssubsets_replicability_DA_tonguedorsum_stool",task_id+100,".RDS"))
