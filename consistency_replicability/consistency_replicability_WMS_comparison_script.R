library(phyloseq)
source("../eval_functions.R")
load(file = "../data/WMSsubsets_replicability.RData")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

comparison <- subsets_consistency_replicability_WMS$CRC_control[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/WMSsubsets_replicability_DA_CRC_control",task_id,".RDS"))

comparison <- subsets_consistency_replicability_WMS$schizophrenia_control[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/WMSsubsets_replicability_DA_schizophrenia_control",task_id,".RDS"))

comparison <- subsets_consistency_replicability_WMS$tonguedorsum_stool[[task_id]]
comparisonDA <- lapply(comparison,oneSimRunGSOwn,epsilon = 1e14)
saveRDS(comparisonDA,file = paste0("../data/WMSsubsets_replicability_DA_tonguedorsum_stool",task_id,".RDS"))
