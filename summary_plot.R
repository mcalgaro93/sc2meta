library(ggplot2)
library(plyr)

### Heatmap all
FDR_summary <- readRDS(file = "./data/summary/FDR_divergence_summary.RDS")
KS_summary <- readRDS(file = "./data/summary/KS_summary.RDS")
WMC_summary <- readRDS(file = "./data/summary/WMC_summary.RDS")
simulations_summary <- readRDS(file = "./data/summary/simulations_summary.RDS")
simulations_time_summary <- readRDS(file = "./data/summary/simulations_time_summary.RDS")
real_data_time_summary <- readRDS(file = "./data/summary/real_data_time_summary.RDS")
enrichment_summary <- readRDS(file = "./data/summary/enrichment_summary.RDS")
  
FDR_summary$value <- FDR_summary$value/nrow(FDR_summary)
KS_summary$value <- KS_summary$value/nrow(KS_summary)
WMC_summary$value <- WMC_summary$value/nrow(WMC_summary)
simulations_summary$value <- simulations_summary$value/nrow(simulations_summary)
real_data_time_summary$value <- real_data_time_summary$value/nrow(real_data_time_summary)
enrichment_summary$value <- enrichment_summary$value/nrow(enrichment_summary)

summary_list <- list("Type I Error - Nominal vs Observed" = FDR_summary,
                     "Type I Error - KS" = KS_summary,
                     "Concordance analysis - WMC" = WMC_summary,
                     "Power - Enrichement analysis" = enrichment_summary,
                     "Computational time" = real_data_time_summary)

summary_df <- ldply(summary_list,.id = "Parameter")

methodlist = c("DESeq2_poscounts",
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
               "songbird",
               "mixMC",
               "MAST",
               "seurat_wilcoxon",
               "scde")

summary_df$method <- factor(summary_df$method, levels = methodlist, ordered = TRUE)

ord <- order(ddply(summary_df[!summary_df$Parameter %in% c("Computational time") & !summary_df$method %in% c("songbird","mixMC"),],.variables = ~ method, .fun = function(x){
  mean(x$value,na.action = "na.omit")
})[,2])

index <- which(levels(summary_df$method) %in% c("mixMC","songbird"))
methodlist_ordered <- methodlist[!methodlist %in% c("songbird","mixMC")]
methodlist_ordered <- c(methodlist_ordered[ord],c("mixMC","songbird"))

fig <- ggplot(summary_df,mapping = aes(x = method, y = Parameter, fill = value)) +
  geom_tile(width = 0.8,height = 0.8) +
  #geom_text(aes(label = round(mean,digits = 2))) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,1)) + 
  coord_equal() +
  theme_minimal() +
  ylab("Measure") + xlab("Method") + labs(fill = "Mean rank") +
  ggtitle(label = "Overall normalized ranking", subtitle = "Ranked methods") +
  scale_y_discrete(limits = rev(levels(summary_df$Parameter))) +
  scale_x_discrete(limits = methodlist_ordered) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        axis.ticks = element_blank())
fig

svg("../fig7.svg",height = 5, width = 8)
fig
dev.off()

fig_vertical <- ggplot(summary_df,mapping = aes(x = Parameter, y = method, fill = value)) +
  geom_tile(width = 0.8,height = 0.8) +
  #geom_text(aes(label = round(mean,digits = 2))) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0,1)) + 
  coord_equal() +
  theme_minimal() +
  xlab("Measure") + ylab("Method") + labs(fill = "Mean rank") +
  ggtitle(label = "Overall normalized ranking", subtitle = "Ranked methods") +
  # scale_x_discrete(limits = rev(levels(summary_df$Parameter))) +
  scale_y_discrete(limits = methodlist_ordered) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        axis.ticks = element_blank())
fig_vertical

svg("../fig7_vertical.svg",height = 8, width = 5)
fig_vertical
dev.off()