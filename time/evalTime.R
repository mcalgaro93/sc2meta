library(reshape2)
source("../colors_17.R")

if(!file.exists("../data/Stool_16S_WMS_times.RDS")) {
  # Times
  all_times <- readRDS("../data/times.RDS") # Not present, too heavy
  # Simulation framework
  load("../data/Stool_16S_wMG_simulation_flow.RData")
  
  library(plyr)
  all_times_method <- lapply(all_times,function(x) x[names(x)!="physeq"]) 
  all_times_method_df <- lapply(all_times_method,function(x) ldply(x,function(y) y["elapsed"]))
  #all_times_method_rank_df <- lapply(all_times_method_df,function(x) return(cbind(x[,".id"],rank(x[,"elapsed"]))))
  all_times_df <- ldply(all_times_method_df)
  colnames(all_times_df) <- c("method","elapsed")
  
  all_times_df <- data.frame(apply(simulation_flow,2,function(col) rep(col,each = 16)),all_times_df)
  all_times_df <- all_times_df[,!(colnames(all_times_df) %in% c("seed"))]
  rownames(all_times_df) <- NULL
  all_times_df <- data.frame(apply(all_times_df,2,unname))
  all_times_df$elapsed <- as.numeric(as.character(all_times_df$elapsed))
  
  rank_all <- ddply(all_times_df, ~ simulation + dataset + distribution + sampleSize + foldEffect + TPR + sparsityEffect + compensation + distribution, function(x) rank(x[,"elapsed"]))
  colnames(rank_all)[9:24] <- as.character(unique(all_times_df$method))
  rank_all <- ddply(rank_all, ~ dataset + distribution + sampleSize + foldEffect + TPR + sparsityEffect + compensation + distribution, function(x) colMeans(x[,9:24]))
  
  rank_agg_dataset <- ddply(rank_all,.variables = ~ dataset, function(x) colMeans(x[,8:23]))
  rank_agg_distribution <- ddply(rank_all,.variables = ~ distribution , function(x) colMeans(x[,8:23]))
  rank_agg_sampleSize <- ddply(rank_all,.variables = ~ sampleSize , function(x) colMeans(x[,8:23]))
  rank_agg_foldEffect <- ddply(rank_all,.variables = ~ foldEffect , function(x) colMeans(x[,8:23]))
  rank_agg_TPR <- ddply(rank_all,.variables = ~ TPR , function(x) colMeans(x[,8:23]))
  rank_agg_sparsityEffect <- ddply(rank_all,.variables = ~ sparsityEffect , function(x) colMeans(x[,8:23]))
  rank_agg_compensation <- ddply(rank_all,.variables = ~ compensation , function(x) colMeans(x[,8:23]))
  rank_agg_distribution <- ddply(rank_all,.variables = ~ distribution , function(x) colMeans(x[,8:23]))
  
  mean_agg_dataset <- ddply(all_times_df,.variables = ~ method + dataset, function(x) mean(x[,"elapsed"]))
  mean_agg_distribution <- ddply(all_times_df,.variables = ~ method + distribution , function(x) mean(x[,"elapsed"]))
  mean_agg_sampleSize <- ddply(all_times_df,.variables = ~ method + sampleSize , function(x) mean(x[,"elapsed"]))
  mean_agg_foldEffect <- ddply(all_times_df,.variables = ~ method + foldEffect , function(x) mean(x[,"elapsed"]))
  mean_agg_TPR <- ddply(all_times_df,.variables = ~ method + TPR , function(x) mean(x[,"elapsed"]))
  mean_agg_sparsityEffect <- ddply(all_times_df,.variables = ~ method + sparsityEffect , function(x) mean(x[,"elapsed"]))
  mean_agg_compensation <- ddply(all_times_df,.variables = ~ method + compensation , function(x) mean(x[,"elapsed"]))
  mean_agg_distribution <- ddply(all_times_df,.variables = ~ method + distribution , function(x) mean(x[,"elapsed"]))
  
  rank_agg_tot <- list(distribution = rank_agg_distribution,
                       sampleSize = rank_agg_sampleSize,
                       TPR = rank_agg_TPR,
                       foldEffect = rank_agg_foldEffect,
                       compensation = rank_agg_compensation,
                       sparsityEffect = rank_agg_sparsityEffect)
  
  rank_agg_tot_melted_list <- lapply(rank_agg_tot,melt)
  rank_agg_tot_melted_list <- lapply(rank_agg_tot_melted_list,function(x){
    colnames(x)[1] <- "sim_variable"
    return(x)
  })
  rank_agg_tot_df <- ldply(rank_agg_tot_melted_list)
  colnames(rank_agg_tot_df) <- c("sim_variable","sim_value","method","value")
  rank_agg_tot_df$sim_variable <- factor(rank_agg_tot_df$sim_variable,levels = unique(rank_agg_tot_df$sim_variable),ordered = TRUE)
  
  rank_agg_tot_df$method <- factor(rank_agg_tot_df$method,levels = c("DESeq2_poscounts",
                                                                     "DESeq2_TMM",
                                                                     "DESeq2_poscounts_zinbwave",
                                                                     "edgeR_TMM_standard",
                                                                     "edgeR_poscounts_standard",
                                                                     "edgeR_TMM_robustDisp",
                                                                     "edgeR_TMM_zinbwave",
                                                                     "limma_voom_TMM",
                                                                     "limma_voom_TMM_zinbwave",
                                                                     "ALDEx2",
                                                                     "mgsZig_CSS",
                                                                     "corncob_LRT",
                                                                     "corncob_wald",
                                                                     "MAST",
                                                                     "seurat_wilcoxon",
                                                                     "scde"), ordered = TRUE)
  
  saveRDS(list("rank_agg_tot_df" = rank_agg_tot_df,"all_times_df" = all_times_df),file = "./data/Stool_16S_WMS_times.RDS")
} else times_list <- readRDS(file = "../data/Stool_16S_WMS_times.RDS")

rank_agg_tot_df <- times_list$rank_agg_tot_df
all_times_df <- times_list$all_times_df

# devtools::install_github("teunbrand/ggnomics")
library(ggnomics)
library(ggpubr)
# svg("../data/heat_all.svg",width = 13,height = 7)
# ggplot(data = rank_agg_tot_df, mapping = aes(x = sim_value, y = method, fill = value)) + 
#   facet_grid(~ sim_variable, scales = "free_x", space = "free_x") + 
#   geom_tile(width = 0.8, height = 0.8) +
#   scale_fill_distiller(palette = "RdYlBu",guide = "colorbar",limits = c(1,16)) +
#   xlab("Simulation framework") + ylab("Method") + labs(fill = "Mean rank") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank())
# dev.off()

library(plyr)
rank_summary <- ddply(rank_agg_tot_df, ~ sim_variable + method, function(x) mean(x$value))
rank_summary <- ddply(rank_summary, ~ method, function(x) mean(x$V1))
colnames(rank_summary) <- c("method","value")

saveRDS(object = rank_summary,file = "../data/summary/simulations_time_summary.RDS")

# all_times_df$method <- relevel(all_times_df$method,ref = "limma_voom_TMM")
# mod1 <- step(glm(elapsed ~ method*sampleSize + dataset + distribution + sampleSize + TPR + foldEffect + compensation + sparsityEffect, data = all_times_df, family = Gamma),direction = "both")
# summary(mod1)
# plot(mod1)

library(ggplot2)
library(ggforce)
ord <- order(ddply(all_times_df,.variables = ~ method,.fun = function(x) mean(x$elapsed))$V1)

all_times_df$method <- factor(all_times_df$method,levels = levels(all_times_df$method)[ord], ordered = TRUE)
library(cowplot)
g_legend <- get_legend(ggplot(all_times_df,mapping = aes(x = sampleSize, y = elapsed, color = method)) + 
  facet_wrap(~ distribution,nrow = 2) +
  geom_point(size = 3) + 
  theme_minimal() +
  labs(x = "Sample Size", y = "Computational time (s)") + 
  scale_color_manual(values = cols, guide = guide_legend(nrow = 3)) + 
  theme(legend.position = "bottom"))

plot_time_zoom <- function(distribution){
  df <- all_times_df[all_times_df$distribution == distribution,]
  ord <- order(ddply(df,~ method, function(x) median(x[,"elapsed"]))[,"V1"])
  df$method <- factor(df$method, levels = levels(df$method)[ord],ordered = TRUE)
  ggplot(all_times_df[all_times_df$distribution == distribution,],mapping = aes(x = sampleSize, y = elapsed, color = method)) + 
    facet_zoom(ylim = c(0,20)) + 
    #facet_wrap(~ distribution,nrow = 2) +
    geom_boxplot() + 
    #theme_minimal() +
    labs(x = "Sample Size", y = "Computational time (s)") + 
    ggtitle(label = paste("Computational times for",distribution,"distribution"),subtitle = "Zoomed-in from 0 to 20 seconds") +
    scale_color_manual(values = cols) + 
    theme(legend.position = "none")
}


svg("../fig_times.svg",width = 15,height = 10) 
plot_grid(
  plot_time_zoom("NB"),
  plot_time_zoom("ZINB"),
  g_legend,nrow = 3,rel_heights = c(1,1,0.3),labels = c("a","b"))
dev.off()

aggr_times_means <- aggregate(elapsed ~ method + sampleSize + distribution, data = all_times_df, FUN = function(x) round(mean(x),2))
aggr_times_sd <- aggregate(elapsed ~ method + sampleSize + distribution, data = all_times_df, FUN = function(x) round(sd(x),2))
dc_means <- dcast(aggr_times_means,formula =  method ~ distribution + sampleSize)
dc_sd <- dcast(aggr_times_sd,formula =  method ~ distribution + sampleSize)
dc_all <- mapply(FUN = function(x,y){
  paste0(x," (",y,")")
},dc_means,dc_sd)
write.csv(dc_all,file = "../data/aggr_times_sd.csv")

dc <- dcast(aggr_times,formula =  method ~ distribution + sampleSize)
dc[match(levels(all_times_df$method),dc$method),]

