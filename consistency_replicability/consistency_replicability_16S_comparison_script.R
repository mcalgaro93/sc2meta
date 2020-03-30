library(phyloseq)
source("../eval_functions.R")
load(file = "../data/16Ssubsets_replicability.RData")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

comparison <- subsets_consistency_replicability_16S$subgingival_supragingival[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/16Ssubsets_replicability_DA_subgingival_supragingival",task_id,".RDS"))

comparison <- subsets_consistency_replicability_16S$gingiva_mucosa[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/16Ssubsets_replicability_DA_gingiva_mucosa",task_id,".RDS"))

comparison <- subsets_consistency_replicability_16S$tonguedorsum_stool[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/16Ssubsets_replicability_DA_tonguedorsum_stool",task_id,".RDS"))
