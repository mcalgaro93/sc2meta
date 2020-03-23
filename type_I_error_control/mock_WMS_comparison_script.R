library(phyloseq)
source("../eval_functions.R")
ps_list <- readRDS(file = "../data/Stool_16S_WMS_mock.RDS")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

ps <- ps_list$Stool_WMS[[task_id]]
Stool_WMS_mockDA <- oneSimRunGSOwn(physeq = ps,epsilon = 1e10)

saveRDS(object = Stool_WMS_mockDA,file = paste0("../mockDA/Stool_WMS_mockDA_",task_id,".RDS"))