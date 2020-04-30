library(ggplot2)
library(plyr)

### Heatmap all
FDR_summary <- readRDS(file = "./data/summary/FDR_divergence_summary.RDS")
KS_summary <- readRDS(file = "./data/summary/KS_summary.RDS")
WMC_summary <- readRDS(file = "./data/summary/WMC_summary.RDS")
simulations_summary <- readRDS(file = "./data/summary/simulations_summary.RDS")
simulations_time_summary <- readRDS(file = "./data/summary/simulations_time_summary.RDS")

summary_list <- list("Nominal vs Observed Type I Error" = FDR_summary,
                     "Type I Error - KS" = KS_summary,
                     "Concordance analysis - WMC" = WMC_summary,
                     "Parametric simulations - pAUROC" = simulations_summary,
                     "Computational time" = simulations_time_summary)

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

ord <- order(ddply(summary_df[summary_df$Parameter != "Computational time",],.variables = ~ method, .fun = function(x){
  mean(x$value,na.action = "na.omit")
})[,2])

summary_df$method <- factor(summary_df$method, levels = methodlist[ord], ordered = TRUE)

fig <- ggplot(summary_df,mapping = aes(x = Parameter, y = method, fill = value)) +
  geom_tile(width = 0.8,height = 0.8) +
  #geom_text(aes(label = round(mean,digits = 2))) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(1,18)) + 
  coord_equal() +
  theme_minimal() +
  xlab("Measure") + ylab("Method") + labs(fill = "Mean rank:") +
  ggtitle(label = "Overall ranking", subtitle = "Ranked methods") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks = element_blank())
fig

# svg("../fig6_incomplete.svg",height = 11, width = 7)
# fig
# dev.off()
