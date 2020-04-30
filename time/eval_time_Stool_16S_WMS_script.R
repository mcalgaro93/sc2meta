source("../eval_functions.R")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

task_id <- as.numeric(args[1])

# Not present in repository. Too heavy files
load(file = paste0("../data/Stool_16S_wMG_split/sim",task_id,".RData"))

times <- lapply(sim,function(s) oneSimRunGSOwn_time(physeq = s$ps,epsilon = 1e14))
saveRDS(times,file = paste0("../data/Stool_16S_wMG_split/times",task_id,".RDS"))