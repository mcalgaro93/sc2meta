### This code was supplied to sbatch, e.g.:
#                                       to be generated     output filedir     log                  errors
#                                             |                       |         |                     |
# Rscript ./power/eval_functions_call.R ./data/sims.RData ./data/evals.RDS > ./data/sims.log 2> ./data/sims.err

source("./eval_functions.R")
register(SerialParam())
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("At least two argument must be supplied (input_file and output_file)", call.=FALSE)
}  
### simulations are required in input
input_file <- args[1]
### name for the output is the second input to supply
output_file <- args[2]

load(file = input_file)
start_time <- Sys.time()
evals <- lapply(X = sim, FUN = function(sim){
  cat(green("Simulation",sim$number,"\n"))
  cat(green("Simulation name:",sim$name,"\n"))
  cat("Overall time:",Sys.time()-start_time,"\n")
  start <- Sys.time()
  eval <- oneSimRunGSOwn(physeq = sim$ps)
  end <- Sys.time()
  cat(magenta("dataset evaluation time:",end-start,"\n"))
  return(eval)
})
saveRDS(evals,file = output_file)