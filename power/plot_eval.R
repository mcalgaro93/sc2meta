library(plyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggpubr)

# Since files containing simulations' evaluations are very memory expensive,
# here we provide the ready to use data.frame with partial AUROC values for each
# set of variables of the simulation framework: evals_AUC_ROC_tot.RDS
# The code to compute it, is reported in this RScript and in the other RScripts
# in the ./power/ directory.

read_data <- function(main_dir = "Stool_16S_WMS"){
  evals_ROC_summary_mean_df <- readRDS(file = paste0("./data/",main_dir,"/evals_ROC_summary_mean_df.RDS"))
  evals_ROC_summary_df <- readRDS(file = paste0("./data/",main_dir,"/evals_ROC_summary_df.RDS"))
  ### Partial Area Under the ROC curve, from 0 to 0.1
  evals_AUC_ROC_tot <- ddply(evals_ROC_summary_df,.variables = ~ method + dataset + distribution + sampleSize + simulation + TPR + foldEffect + compensation + sparsityEffect,
                         function(x) sum(x$tpr[x$fpr<=0.1]/10))

  evals_AUC_ROC_tot$method <- factor(evals_AUC_ROC_tot$method,levels = unique(evals_AUC_ROC_tot$method))
  evals_AUC_ROC_tot <- evals_AUC_ROC_tot[!is.na(evals_AUC_ROC_tot$method),]
  return(evals_AUC_ROC_tot)
}
# 
# Stool <- read_data(main_dir = "Stool_16S_WMG")
# TongueDorsum <- read_data(main_dir = "TD_16S_WMG")
# BritoIL <- read_data(main_dir = "BritoIL_Stool_Oral")
# 
# ### rbind all evals_AUC_ROC_tot from each dataset evaluated
# evals_AUC_ROC_tot <- data.frame(rbind(Stool,TongueDorsum,BritoIL))
# saveRDS(evals_AUC_ROC_tot,file = "./data/evals_AUC_ROC_tot.RDS")

evals_AUC_ROC_tot <- readRDS(file = "./data/evals_AUC_ROC_tot.RDS")
evals_AUC_ROC_tot_corncob <- readRDS(file = "./data/evals_AUC_ROC_tot_corncob.RDS")
evals_AUC_ROC_tot <- rbind(evals_AUC_ROC_tot,evals_AUC_ROC_tot_corncob)
# evals_AUC_ROC_tot$distribution_sparsityEffect <- paste0(evals_AUC_ROC_tot$distribution,"\n",evals_AUC_ROC_tot$sparsityEffect)
# evals_AUC_ROC_tot$foldEffect_sampleSize <- paste0(evals_AUC_ROC_tot$foldEffect,"\n",evals_AUC_ROC_tot$sampleSize)
# evals_AUC_ROC_tot$TPR_sampleSize <- paste0(evals_AUC_ROC_tot$TPR,"\n",evals_AUC_ROC_tot$sampleSize)

# Text formatting for some variables
evals_AUC_ROC_tot$distribution_sparsityEffect <- paste0(evals_AUC_ROC_tot$distribution,"-",evals_AUC_ROC_tot$sparsityEffect)
evals_AUC_ROC_tot$foldEffect_sampleSize <- paste0(evals_AUC_ROC_tot$foldEffect,"-",evals_AUC_ROC_tot$sampleSize)
evals_AUC_ROC_tot$TPR_sampleSize <- paste0(evals_AUC_ROC_tot$TPR,"-",evals_AUC_ROC_tot$sampleSize)

# mean(1- pAUROC) computation in order to rank methods
evals_AUC_ROC_mean <- ddply(evals_AUC_ROC_tot, ~ dataset + distribution + sampleSize + TPR + foldEffect + compensation + sparsityEffect + distribution_sparsityEffect + foldEffect_sampleSize + TPR_sampleSize + method, function(x) mean(1-x[,"V1"]))
rank_all <- ddply(evals_AUC_ROC_mean, ~ dataset + distribution + sampleSize + TPR + foldEffect + compensation + sparsityEffect + distribution_sparsityEffect + foldEffect_sampleSize + TPR_sampleSize, function(x) rank(x[,"V1"]))
mean_all <- ddply(evals_AUC_ROC_mean, ~ dataset + distribution + sampleSize + TPR + foldEffect + compensation + sparsityEffect + distribution_sparsityEffect + foldEffect_sampleSize + TPR_sampleSize, function(x) x[,"V1"])

# Function to properly aggregate results for plotting

aggregation <- function(rank_all,evals_AUC_ROC_tot){
  colnames(rank_all)[11:26] <- levels(evals_AUC_ROC_tot$method)
  rank_all$technology <- factor(ifelse(grepl(x = rank_all$dataset,pattern = "16S"),"16S","WMS"))
  
  # rank_agg_distribution <- ddply(rank_all, ~ distribution, function(x) colMeans(x[,9:22]))
  # rank_agg_dataset <- ddply(rank_all, ~ distribution + technology + dataset, function(x) colMeans(x[,9:22]) )
  # rank_agg_sampleSize <- ddply(rank_all, ~ distribution + technology + sampleSize, function(x) colMeans(x[,9:22]) )
  # rank_agg_foldEffect <- ddply(rank_all, ~ distribution + technology + foldEffect, function(x) colMeans(x[,9:22]) )
  # rank_agg_TPR <- ddply(rank_all, ~ distribution + technology + TPR, function(x) colMeans(x[,9:22]) )
  # rank_agg_sparsityEffect <- ddply(rank_all, ~ distribution + technology + sparsityEffect, function(x) colMeans(x[,9:22]) )
  # rank_agg_compensation <- ddply(rank_all, ~ distribution + technology + compensation, function(x) colMeans(x[,9:22]) )
  
  rank_agg_technology <- ddply(rank_all, ~ distribution + technology, function(x) colMeans(x[,11:26]) )
  rank_agg_technology <- ddply(rank_agg_technology, ~ technology, function(x) colMeans(x[,3:18]) )
  
  rank_agg_distribution <- ddply(rank_all, ~ distribution + technology, function(x) colMeans(x[,11:26]) )
  rank_agg_distribution <- ddply(rank_agg_distribution, ~ distribution, function(x) colMeans(x[,3:18]) )
  
  rank_agg_dataset <- ddply(rank_all, ~ distribution + dataset, function(x) colMeans(x[,11:26]) )
  rank_agg_dataset <- ddply(rank_agg_dataset, ~ dataset, function(x) colMeans(x[,3:18]) )
  
  rank_agg_sampleSize <- ddply(rank_all, ~ distribution + technology + sampleSize, function(x) colMeans(x[,11:26]) )
  rank_agg_sampleSize <- ddply(rank_agg_sampleSize, ~ sampleSize, function(x) colMeans(x[,4:19]) )
  
  rank_agg_foldEffect <- ddply(rank_all, ~ distribution + technology + foldEffect, function(x) colMeans(x[,11:26]) )
  rank_agg_foldEffect <- ddply(rank_agg_foldEffect, ~ foldEffect, function(x) colMeans(x[,4:19]) )
  
  rank_agg_TPR <- ddply(rank_all, ~ distribution + technology + TPR, function(x) colMeans(x[,11:26]) )
  rank_agg_TPR <- ddply(rank_agg_TPR, ~ TPR, function(x) colMeans(x[,4:19]) )
  
  rank_agg_sparsityEffect <- ddply(rank_all, ~ distribution + technology + sparsityEffect, function(x) colMeans(x[,11:26]) )
  rank_agg_sparsityEffect <- ddply(rank_agg_sparsityEffect, ~ sparsityEffect, function(x) colMeans(x[,4:19]) )
  
  rank_agg_compensation <- ddply(rank_all, ~ distribution + technology + compensation, function(x) colMeans(x[,11:26]) )
  rank_agg_compensation <- ddply(rank_agg_compensation, ~ compensation, function(x) colMeans(x[,4:19]) )
  
  rank_agg_distribution_sparsityEffect <- ddply(rank_all, ~ distribution + technology + sparsityEffect + distribution_sparsityEffect, function(x) colMeans(x[,11:26]) )
  rank_agg_distribution_sparsityEffect <- ddply(rank_agg_distribution_sparsityEffect, ~ distribution_sparsityEffect, function(x) colMeans(x[,5:20]) )
  
  rank_agg_foldEffect_sampleSize <- ddply(rank_all, ~ distribution + technology + foldEffect + sampleSize + foldEffect_sampleSize, function(x) colMeans(x[,11:26]) )
  rank_agg_foldEffect_sampleSize <- ddply(rank_agg_foldEffect_sampleSize, ~ foldEffect_sampleSize, function(x) colMeans(x[,6:21]) )
  
  rank_agg_TPR_sampleSize <- ddply(rank_all, ~ distribution + technology + TPR + sampleSize + TPR_sampleSize, function(x) colMeans(x[,11:26]) )
  rank_agg_TPR_sampleSize <- ddply(rank_agg_TPR_sampleSize, ~ TPR_sampleSize, function(x) colMeans(x[,6:21]) )
  
  rank_agg_tot <- list(technology = rank_agg_technology,
                       distribution = rank_agg_distribution,
                       dataset = rank_agg_dataset,
                       sampleSize = rank_agg_sampleSize,
                       TPR = rank_agg_TPR,
                       foldEffect = rank_agg_foldEffect,
                       compensation = rank_agg_compensation,
                       sparsityEffect = rank_agg_sparsityEffect,
                       distribution_sparsityEffect = rank_agg_distribution_sparsityEffect,
                       foldEffect_sampleSize = rank_agg_foldEffect_sampleSize,
                       TPR_sampleSize = rank_agg_TPR_sampleSize)
  
  rank_agg_tot_melted_list <- lapply(rank_agg_tot,melt)
  rank_agg_tot_melted_list <- lapply(rank_agg_tot,function(x){
    colnames(x)[1] <- "sim_variable"
    return(x)
  })
  rank_agg_tot_df <- ldply(rank_agg_tot_melted_list)
  rank_agg_tot_df_melted <- melt(rank_agg_tot_df)
  colnames(rank_agg_tot_df_melted) <- c("sim_variable","sim_value","method","value")
  rank_agg_tot_df_melted$sim_variable <- factor(rank_agg_tot_df_melted$sim_variable,
                                                levels = unique(rank_agg_tot_df_melted$sim_variable),
                                                labels = c(unique(rank_agg_tot_df_melted$sim_variable)[1:8],
                                                           "distribution\nsparsityEffect",
                                                           "foldEffect\nsampleSize",
                                                           "TPR\nsampleSize"),
                                                ordered = TRUE)
  rank_agg_tot_df_melted
}

rank_agg_tot_df_melted <- aggregation(rank_all = rank_all,evals_AUC_ROC_tot = evals_AUC_ROC_tot)
mean_agg_tot_df_melted <- aggregation(rank_all = mean_all,evals_AUC_ROC_tot = evals_AUC_ROC_tot)

rank_agg_tot_df_melted$mean <- 1-mean_agg_tot_df_melted$value

library(ggpubr)
univariate <- c("technology",
                "distribution",
                "dataset",
                "sampleSize",
                "TPR",
                "foldEffect",
                "compensation",
                "sparsityEffect")
bivariate <- c("distribution\nsparsityEffect",
               "foldEffect\nsampleSize",
               "TPR\nsampleSize")

simulations_summary <- ddply(rank_agg_tot_df_melted,.variables = ~ method,function(x) mean(x[,"value"]))
colnames(simulations_summary) <- c("method","value")
saveRDS(simulations_summary, file = "./data/summary/simulations_summary.RDS")

ord <- order(ddply(rank_agg_tot_df_melted,.variables = ~ method,function(x) mean(x[,"value"]))$V1)

a1 <- ggplot(data = rank_agg_tot_df_melted[rank_agg_tot_df_melted$sim_variable %in% univariate,], mapping = aes(x = sim_value, y = method, fill = value)) + 
  facet_grid(~ sim_variable, scales = "free_x", space = "free_x") + 
  geom_tile(width = 0.8, height = 0.8) +
  geom_text(aes(label = round(mean*100,digits = 0))) +
  scale_y_discrete(limits = levels(rank_agg_tot_df_melted$method)[ord]) +
  scale_fill_distiller(palette = "RdYlBu",guide = "colorbar",limits = c(1,16)) +
  xlab("Variable") + ylab("Method") + labs(fill = "Mean rank") +
  ggtitle(label = "Simulation framework", subtitle = "Ranked methods by partial-AUROC curve") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

a2 <- ggplot(data = rank_agg_tot_df_melted[rank_agg_tot_df_melted$sim_variable %in% bivariate,], mapping = aes(x = sim_value, y = method, fill = value)) + 
  facet_grid(~ sim_variable, scales = "free_x", space = "free_x") + 
  geom_tile(width = 0.8, height = 0.8) +
  geom_text(aes(label = round(mean*100,digits = 0))) +
  scale_fill_distiller(palette = "RdYlBu",guide = "colorbar",limits = c(1,16)) +
  scale_y_discrete(limits = levels(rank_agg_tot_df_melted$method)[ord]) +
  xlab("Variable") + ylab("Method") + labs(fill = "Mean rank") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank()
        # strip.text = element_text(color = "white")
        ) +
  cowplot::panel_border(colour = "black",size = 1,linetype = 2)

a3 <- ggplot(data = rank_agg_tot_df_melted[rank_agg_tot_df_melted$sim_variable %in% bivariate,], mapping = aes(x = sim_value, y = method, fill = value)) + 
  facet_grid(~ sim_variable, scales = "free_x", space = "free_x") + 
  geom_tile(width = 0.8, height = 0.8) +
  geom_text(aes(label = round(mean*100,digits = 0))) +
  scale_fill_distiller(palette = "RdYlBu",guide = "colorbar",limits = c(1,16)) +
  scale_y_discrete(limits = levels(rank_agg_tot_df_melted$method)[ord]) +
  xlab("Variable") + ylab("Method") + labs(fill = "Mean rank") +
  ggtitle(label = "Simulation framework - Couple of variables", subtitle = "Ranked methods by partial-AUROC curve") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank()
        # strip.text = element_text(color = "white")
  ) +
  panel_border(colour = "black",size = 1,linetype = 2)

fig <- plot_grid(a1,a2,align = "h",rel_widths = c(2,1.1))
fig
