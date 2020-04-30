# file_data_in = "./parametric_simulations/data/BritoIL_Stool_Oral_simulations.RData"
# out_dir = "./parametric_simulations/data/BritoIL_Stool_Oral_split"
# nsplits = 30

data_splitter <- function(file_data_in, out_dir, nsplits){
  cat("Loading simulations","\n")
  load(file = file_data_in)
  cat("Splitting",length(sims),"simulations in",nsplits,"blocks","\n")
  indexes <- seq(1,length(sims),length(sims)/nsplits)
  indexes <- c(indexes,length(sims)+1)
  ind <- NULL
  for(i in 1:(length(indexes)-1)){
    ind <- cbind(ind,(indexes[i]):(indexes[i+1]-1))
  } 
  for(i in 1:ncol(ind))
  {
    cat("Splitting simulations... Saving block", i, "of", ncol(ind),"\n")
    sim <- sims[ind[,i]]
    save(sim,file = paste0(out_dir,"/sim",i,".RData")[1])
  }
}

data_merger <- function(file_dir,nsplits){
  cat("Merging simulated data results:","\n")
  evals <- list()
  # ind <- length(list.files(path = file_dir,pattern = ".RDS"))
  for(i in 1:nsplits)
  {
    cat("loading evals",i,"\n")
    sim <- readRDS(file = paste0(file_dir,i,".RDS"))
    evals <- c(evals,sim)
  }
  cat("Saving merged data","\n")
  saveRDS(evals,file = paste0(file_dir,".RDS"))
}

# Example:
# data_merger(file_dir="./parametric_simulations/data/TD_16S_wMG_split/evals_sim_epsilon1e14_",nsplits=30)



