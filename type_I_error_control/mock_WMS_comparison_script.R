library(phyloseq)
source("../eval_functions.R")
load(file = "../data/Stool_16S_WMS.RData")
mock_df <- read.table(file = "../data/Stool_16S_WMS_mock_df.tsv",sep = "\t")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

ps <- Stool_16S_WMS$Stool_WMS
ps@sam_data$grp <- unlist(mock_df[task_id,])
Stool_WMS_mockDA <- oneSimRunGSOwn(physeq = ps,epsilon = 1e10)

saveRDS(object = Stool_WMS_mockDA,file = paste0("../mockDA/Stool_WMS_mockDA_",task_id,".RDS"))