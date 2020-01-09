library(plyr)
library(AUC)

evalPVals <- function(resi, alpha = 0.05, pvalsType = "adjP", rawPvalsType = "rawP") {
  # Rarely a fit has failed, then we return 0 for sens, 1 for specificity and
  # NA for FDR, AUC and the lib/cons areas
  if (!is.matrix(resi)) {
    cat("Resi is not a matrix! \n")
    return(c(Sensitivity = 0, Specificity = 1, FDR = 0, AUC = 0.5))
  }
  # Some DESeq2 results (for example) had NA adjusted p-values Replace NA
  # $adjP values to highest possible value (1.0)
  # resi[is.na(resi[, pvalsType]), pvalsType] <- 1
  # Or just count them
  NA_prop <- sum(is.na(resi[, pvalsType]))/nrow(resi)
  
  # Evaluate detection performance.
  wh.pred = (resi[, pvalsType] < alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("TP", rownames(resi))
  # Calc number of differentially abundant taxa
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  # Sensitivity: True Positives divided by all positives (sum of true
  # positives and false negatives)
  Sensitivity = TPs/(TPs + FNs)
  # Specificity: True negatives divided by all negatives (sum of true
  # negatives and false positives)
  Specificity = TNs/(TNs + FPs)
  # false discovery rate: false positives divided by all detected positives
  FDR = if ((TPs + FPs) == 0) 
    0 else FPs/(TPs + FPs)
  
  # If no true positives, return NA's for irrelevant measures
  wh.truth = (1:nrow(resi) %in% wh.TP)
  
  # IF AUC cannot be calculated, return NA
  rocObj = try(AUC::roc(1 - resi[, pvalsType], factor(as.numeric(wh.truth))))
  return(c("NA_Proportion" = NA_prop, 
           Sensitivity = Sensitivity, 
           Specificity = Specificity, 
           FDR = FDR, 
           AUC = ifelse(class(rocObj)[1] == "try-error", NA, AUC::auc(rocObj))))
}  # END - function: evalPVals

df_creator <- function(evals_file, sim_flow_file, out_dir){
  cat("Reading evals","\n")
  evals <- readRDS(file = evals_file)
  load(file = sim_flow_file)
  cat("Creating data.frame from evals","\n")
  eval_stats <- ldply(.data = evals,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      evalPVals(resi = m$pValMat,alpha = 0.05,pvalsType = "adjP",rawPvalsType = "rawP")
    })
  })
  colnames(eval_stats) <- c("method",colnames(eval_stats)[-1])
  nmethods <- length(unique(eval_stats$method))
  simulation_flow_df <- apply(simulation_flow[1:length(evals),], 2, function(col) sapply(col,function(cell) rep(cell,each = nmethods)))
  evals_stats_df <- data.frame(eval_stats,simulation_flow_df) 
  
  cat("Computing ROC from pVals","\n")
  eval_ROC <- ldply(.data = evals,.fun = function(methods){
    ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP",x = rownames(m$pValMat))))
      cbind(fpr = ROC$fpr, tpr = ROC$tpr)
    })
  })
  colnames(eval_ROC) <- c("method",colnames(eval_ROC)[-1])
  lengths_ROC <- ldply(.data = evals,.fun = function(methods){
    sum(ldply(.data = methods,.fun = function(m){
      ROC <- AUC::roc(predictions = 1-m$pValMat[,"adjP"], labels = as.factor(grepl(pattern = "TP",x = rownames(m$pValMat))))
      return(length(ROC$tpr))
    })$V1)
  })
  simulation_flow_ROC_df <- apply(simulation_flow[1:length(evals),], 2, function(col) unlist(mapply(col,lengths_ROC$V1,FUN = function(cell,times) rep(x = cell,times)),use.names = FALSE))
  evals_ROC_df <- cbind(eval_ROC,simulation_flow_ROC_df) 
  cat("Summarizing ROC values","\n")
  evals_ROC_summary_df <- ddply(.data = evals_ROC_df[,-ncol(evals_ROC_df)],.variables = ~ 
                                  method + 
                                  dataset + 
                                  distribution + 
                                  sampleSize +
                                  simulation +
                                  TPR + 
                                  foldEffect + 
                                  compensation + 
                                  sparsityEffect, .fun = function(x){
                                    support <- seq(0,1,length.out = 101)
                                    fpr_tpr <- data.frame(fpr = 0, tpr = 0)
                                    for(i in 2:length(support)){
                                      fpr_s <- support[i-1]
                                      if(sum(x$fpr>=support[i-1] & x$fpr<support[i]) > 0)
                                        tpr_s <- mean(x$tpr[x$fpr>=support[i-1] & x$fpr<support[i]])
                                      else tpr_s <- fpr_tpr[i-2,2]
                                      fpr_tpr[i-1,] <- c(fpr_s,tpr_s)
                                    }
                                    fpr_tpr[1,] <- c(0,0)
                                    fpr_tpr[length(support),] <- c(1,1)
                                    return(fpr_tpr)
                                  })
  evals_ROC_summary_mean_df <- ddply(.data = evals_ROC_summary_df,.variables = ~ 
                                  method + 
                                  dataset + 
                                  distribution + 
                                  sampleSize +
                                  TPR + 
                                  foldEffect + 
                                  compensation + 
                                  sparsityEffect +
                                  fpr, .fun = function(x){
                                    tpr = mean(x$tpr)
                                    se = sqrt(var(x$tpr))
                                    return(data.frame(tpr = tpr, se = se))
                                  })
  cat("Saving data","\n")
  saveRDS(evals_stats_df,file = paste0(out_dir,"evals_stats_df.RDS"))
  saveRDS(evals_ROC_summary_df,file = paste0(out_dir,"evals_ROC_summary_df.RDS"))
  saveRDS(evals_ROC_summary_mean_df,file = paste0(out_dir,"evals_ROC_summary_mean_df.RDS"))
}

### Example code to generate power data.frames 
### The simulation files are heavy, for this reason they are not saved in github
### However the final data.frame is available.
df_creator(evals_file="./data/Stool_16S_WMS_evals_simulations.RDS",
           sim_flow_file="./data/Stool_16S_WMS_simulation_flow.RData",
           out_dir=".data/Stool_16S_WMS/")

df_creator(evals_file="./data/TD_16S_WMS_evals_simulations.RDS",
           sim_flow_file="./data/TD_16S_WMS_simulation_flow.RData",
           out_dir=".data/TD_16S_WMS/")

df_creator(evals_file="./data/BritoIL_Stool_Oral_simulation_evals_sim.RDS",
           sim_flow_file="./data/BritoIL_Stool_Oral_simulation_flow.RData",
           out_dir=".data/BritoIL_Stool_Oral/")
