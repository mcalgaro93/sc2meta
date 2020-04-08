library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(phyloseq)
library(vegan)
library(reshape2)
library(ffpe)
library(cowplot)
library(scales)
library(ggdendro)


# 16S Data retrieval
generate_dataset_16S <- function(A,B,seed,prop){
  V35_all <<- V35() %>% subset(select = HMP_BODY_SUBSITE %in% c(A,B)) %>% as_phyloseq()
  # group B
  psB <<- subset_samples(V35_all, HMP_BODY_SUBSITE == B)
  # Run Center WUGC to lower numerosity
  psB <<- subset_samples(psB, RUN_CENTER == "WUGC")
  # Keep one sample per individual
  psB <<- subset_samples(psB, !duplicated(RSID))
  # Pruning
  psB <<- prune_samples(sample_sums(psB) >= 10^3, psB)
  # Group A
  psA <<- subset_samples(V35_all, HMP_BODY_SUBSITE == A)
  psA <<- subset_samples(psA, RUN_CENTER == "WUGC")
  psA <<- subset_samples(psA, !duplicated(RSID))
  psA <<- prune_samples(sample_sums(psA) >= 10^3, psA)
  
  grp1 <<- grp2 <<- c()
  set.seed(seed)
  for(i in 1:length(psB@sam_data$RSID)){
    if(psB@sam_data$RSID[i] %in% psA@sam_data$RSID){
      # if(sample(c(1,2),1,prob = c(0.37,0.63)) == 1){
      if(sample(c(1,2),1,prob = prop) == 1){
        grp1 <<- c(grp1,psB@sam_data$RSID[i])
      } else grp2 <<- c(grp2,psB@sam_data$RSID[i])
    } else grp1 <<- c(grp1,psB@sam_data$RSID[i])
  }
  ps <<- merge_phyloseq(subset_samples(psB,RSID %in% grp1),subset_samples(psA,RSID %in% grp2))
  ps <<- filter_taxa(ps,function(x) sum(x>0)>0,1)
  return(ps)
}

# WMS data retrieval
generate_dataset_WMS <- function(psA,psB,seed,prop){
  grp1 <<- grp2 <<- c()
  set.seed(seed)
  for(i in 1:length(psA@sam_data$subjectID)){
    if(psA@sam_data$subjectID[i] %in% psB@sam_data$subjectID){
      if(sample(c(1,2),1,prob = prop) == 1){
        grp1 <<- c(grp1,psA@sam_data$subjectID[i])
      } else grp2 <<- c(grp2,psA@sam_data$subjectID[i])
    } else grp1 <<- c(grp1,psA@sam_data$subjectID[i])
  }
  ps <<- merge_phyloseq(subset_samples(psA,subjectID %in% grp1),subset_samples(psB,subjectID %in% grp2))
  ps <<- filter_taxa(ps,function(x) sum(x>0)>0,1)
  return(ps)
}

# Alpha diversity
alpha_16S <- function(ps){
  alpha_shannon <- vegan::diversity(x = ps@otu_table@.Data, MARGIN = 2, index = "shannon")
  ps@sam_data$Shannon <- alpha_shannon
  HMP_BODY_SUBSITE_medians <- ddply(.data = data.frame(ps@sam_data),.variables = "HMP_BODY_SUBSITE",.fun = function(x) median(x$Shannon))
  ord <- HMP_BODY_SUBSITE_medians$HMP_BODY_SUBSITE[order(HMP_BODY_SUBSITE_medians$V1)]
  cols <- hue_pal()(2)
  names(cols) <- unique(ps@sam_data$HMP_BODY_SUBSITE)
  num <- as.vector(table(ps@sam_data$HMP_BODY_SUBSITE)[ord])
  g_alpha <- ggboxplot(data.frame(ps@sam_data),x = "HMP_BODY_SUBSITE", y = "Shannon",order = ord, varwidth = TRUE, color = "HMP_BODY_SUBSITE",ggtheme = theme_minimal()) + 
    geom_hline(yintercept = sum(HMP_BODY_SUBSITE_medians$V1*num/sum(num)), color = "black", lty = 2) +
    theme(axis.text.x = element_text(hjust = 1,angle = 45),
          legend.position = "none",
          axis.title.x = element_blank()) +
    ggtitle("Alpha diversity",subtitle = "Shannon") +
    stat_compare_means() +
    scale_color_manual(values = cols)
  g_alpha
}

alpha_WMS <- function(ps, variable){
  alpha_shannon <- vegan::diversity(x = ps@otu_table@.Data, MARGIN = 2, index = "shannon")
  ps@sam_data$Shannon <- alpha_shannon
  study_condition_medians <- ddply(.data = data.frame(ps@sam_data),.variables = variable,.fun = function(x) median(x$Shannon))
  ord <- study_condition_medians[order(study_condition_medians$V1),variable]
  cols <- hue_pal()(2)
  names(cols) <- unlist(unname(c(unique(ps@sam_data[,variable]))))
  num <- as.vector(table(ps@sam_data[,variable])[ord])
  g_alpha <- ggboxplot(data.frame(ps@sam_data),x = variable, y = "Shannon",order = ord, varwidth = TRUE, color = variable,ggtheme = theme_minimal()) + 
    geom_hline(yintercept = sum(study_condition_medians$V1*num/sum(num)), color = "black", lty = 2) +
    theme(axis.text.x = element_text(hjust = 1,angle = 45),
          legend.position = "none",
          axis.title.x = element_blank()) +
    ggtitle("Alpha diversity",subtitle = "Shannon") +
    stat_compare_means() +
    scale_color_manual(values = cols)
  g_alpha
}

#MDS
compute_MDS <- function(ps, normalization = "none" , method = "MDS", distance = "bray", color = "RUN_CENTER", names = "HMP_BODY_SUBSITE", ellipse = TRUE){
  if(normalization == "rarefy"){
    ps <- rarefy_even_depth(ps,rngseed = 123)
  } else if(normalization == "TSS"){
    ps@otu_table@.Data <- apply(ps@otu_table@.Data,2,function(col) col/sum(col)*median(colSums(ps@otu_table@.Data))) 
  } else {}

  ordination <- ordinate(ps, method = method, distance = distance) 
  ps@sam_data$A1 <- ordination$vectors[,1]
  ps@sam_data$A2 <- ordination$vectors[,2]
  p_o <- ggplot(data.frame(ps@sam_data), aes(x = A1, y = A2, color = eval(parse(text = color)))) +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle(label = paste0(method," on ",paste(unlist(unique(ps@sam_data[,names])),collapse = " vs ")),subtitle = paste0("Normalization: ",normalization,", Distance: ", distance)) +
    labs(color = color) + xlab(paste0("PCoA.1 (",round(ordination$values$Relative_eig[1]*100,digits = 2),"%)")) +
    ylab(paste0("PCoA.2 (",round(ordination$values$Relative_eig[2]*100,digits = 2),"%)"))
  if(ellipse) p_o <- p_o + stat_ellipse()
  
  p_l <- get_legend(ggplot(data.frame(ps@sam_data), aes(x = A1, y = A2, color = eval(parse(text = color)))) +
                      geom_point(size = 5) +
                      theme_minimal() +
                      theme(legend.position = "bottom") +
                      guides(color=guide_legend(nrow=2,byrow=TRUE)) +
                      labs(color = color))
  
  return(list(plot = p_o,legend = p_l))
}

compute_concordance <- function(ps_fitted_list){
  conc_df <- NULL
  for(i in 1:length(ps_fitted_list)){ # i in 1:100 comparisons
    # adjusted pval extraction
    cat("Comparison",i,"\n")
    adjP_df1 <- llply(.data = ps_fitted_list[[i]]$Subset1,.fun = function(method) method$pValMat[!is.na(method$pValMat[,2]) & method$pValMat[,2]<0.1,2]) # first half data
    adjP_df2 <- llply(.data = ps_fitted_list[[i]]$Subset2,.fun = function(method) method$pValMat[!is.na(method$pValMat[,2]) & method$pValMat[,2]<0.1,2]) # second half data
    # pval extraction
    P_df1 <- llply(.data = ps_fitted_list[[i]]$Subset1,.fun = function(method){
      if("pValMat" %in% names(method)){
        out <- method$pValMat[!is.na(method$pValMat[,1]) & method$pValMat[,1]<1,1]
      } else {
        out <- abs(method[,2])
        names(out) <- rownames(method)
      }
      return(out)
    }) # first half data
    P_df2 <- llply(.data = ps_fitted_list[[i]]$Subset2,.fun = function(method){
      if("pValMat" %in% names(method)){
        out <- method$pValMat[!is.na(method$pValMat[,1]) & method$pValMat[,1]<1,1]
      } else {
        out <- abs(method[,2])
        names(out) <- rownames(method)
      }
      return(out)
    })  # second half data
    
    nmethods <- length(names(ps_fitted_list[[i]]$Subset1))
    for(j in 1:nmethods){ # j in method names
      cat("Mehod",names(ps_fitted_list[[i]]$Subset1)[j],"with \n")
      for(k in 1:nmethods){ # k in method names again
        cat("\t",names(ps_fitted_list[[i]]$Subset1)[k],"\n")
        if(j != k){ # BMC computation
          # BMC for Subset1 data
          conc_subset1 <- data.frame(CATplot(vec1 = P_df1[[j]],vec2 = P_df1[[k]],make.plot = FALSE,maxrank = 100), 
                                     method1 = names(ps_fitted_list[[i]]$Subset1)[j], 
                                     method2 = names(ps_fitted_list[[i]]$Subset1)[k],
                                     #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                                     #ndisc_0.1_method2 = length(adjP_df1[[k]]),
                                     nfeatures = ifelse(test = (j < nmethods),yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]$pValMat),no = nrow(ps_fitted_list[[i]]$Subset1[[j]])),
                                     comparison = i,
                                     subset = "1")
          # BMC for Subset2 data
          conc_subset2 <- data.frame(CATplot(vec1 = P_df2[[j]],vec2 = P_df2[[k]],make.plot = FALSE,maxrank = 100), 
                                     method1 = names(ps_fitted_list[[i]]$Subset2)[j], 
                                     method2 = names(ps_fitted_list[[i]]$Subset2)[k],
                                     #ndisc_0.1_method1 = length(adjP_df2[[j]]),
                                     #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                                     nfeatures = ifelse(test = (j < nmethods),yes = nrow(ps_fitted_list[[i]]$Subset2[[j]]$pValMat),no = nrow(ps_fitted_list[[i]]$Subset2[[j]])),
                                     comparison = i,
                                     subset = "2")
          conc <- rbind(conc_subset1,conc_subset2)
        } else {
          # WMC computed between Subset1 and Subset2
          conc <- data.frame(CATplot(vec1 = P_df1[[j]],vec2 = P_df2[[k]],make.plot = FALSE,maxrank = 100), 
                             method1 = names(ps_fitted_list[[i]]$Subset1)[j], 
                             method2 = names(ps_fitted_list[[i]]$Subset2)[k],
                             #ndisc_0.1_method1 = length(adjP_df1[[j]]),
                             #ndisc_0.1_method2 = length(adjP_df2[[k]]),
                             nfeatures = mean(ifelse(test = (j < nmethods),yes = nrow(ps_fitted_list[[i]]$Subset1[[j]]$pValMat),no = nrow(ps_fitted_list[[i]]$Subset1[[j]])),
                                              ifelse(test = (j < nmethods),yes = nrow(ps_fitted_list[[i]]$Subset2[[k]]$pValMat),no = nrow(ps_fitted_list[[i]]$Subset2[[k]]))),
                             comparison = i,
                             subset = "1vs2")
        }
        conc_df <- rbind(conc_df,conc)
      }
    }
  }
  return(conc_df)
}

### Compute AUC and AOC of concordance distribution 
# This method comes from p-value analysis
AucAocFun <- function(cVals, nfeatures, threshold = 0, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line
  MaxArea <- threshold # Total Area over = under the y=x line
  estimated <- cVals # Estimated values
  theoretic <- seq_along(cVals)/unique(nfeatures) # y=x values
  
  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }
  
  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated # Values over the y=x line
  conservArea <- sum(-diffPVals[indConserv])/MaxArea # Area over the y = x
  liberalArea <- sum(diffPVals[!indConserv])/MaxArea # Area under the y = x
  
  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
      conservArea)
  
}

# Concordance plot function

gheat <- function(AUC_AOC_between_methods,concordance_df_summary,tech,comp){
  # Filtering
  AUC_AOC_between_methods <- AUC_AOC_between_methods[AUC_AOC_between_methods$tech == tech & AUC_AOC_between_methods$comp == comp,]
  concordance_df_summary <- concordance_df_summary[concordance_df_summary$tech == tech & concordance_df_summary$comp == comp,]
  forlegend <- AUC_AOC_between_methods
  forlegend$method1 <- factor(forlegend$method1,levels = c("DESeq2_poscounts",
                                                           "DESeq2_TMM",
                                                           "DESeq2_poscounts_zinbwave",
                                                           "edgeR_TMM_standard",
                                                           "edgeR_poscounts_standard",
                                                           "edgeR_TMM_robustDisp",
                                                           "edgeR_TMM_zinbwave",
                                                           "limma_voom_TMM",
                                                           "limma_voom_TMM_zinbwave",
                                                           "ALDEx2",
                                                           "songbird",
                                                           "mgsZig_CSS",
                                                           "corncob_LRT",
                                                           "corncob_wald",
                                                           "MAST",
                                                           "seurat_wilcoxon",
                                                           "scde"),ordered = TRUE)
  g_legend_dendrogram <- get_legend(ggplot() + 
                                      geom_point(data=forlegend, aes(x = method1, y = 1, color = method1),size = 5) +
                                      scale_color_manual(values = cols) +
                                      theme_minimal() +
                                      theme(legend.position = "bottom") +
                                      guides(color = guide_legend(title = "Methods:",title.position = "left",nrow = 3)))
  
  # Clustering
  dist_matrix <- dcast(data = AUC_AOC_between_methods, formula = method1 ~ method2,value.var = "conservArea")
  dist_df <- dist_matrix[,2:ncol(dist_matrix)]
  rownames(dist_df) <- colnames(dist_df)
  distances <- as.dist(1-dist_df)
  hc <- hclust(d = distances)
  # Area extraction
  area <- apply(concordance_df_summary,1,function(x){
    area <- AUC_AOC_between_methods$conservArea[AUC_AOC_between_methods$method1 == x["method1"] & AUC_AOC_between_methods$method2 == x["method2"]]
    return(1-area)
  })
  concordance_df_summary_area <- cbind(concordance_df_summary,area = area)
  # As factor
  concordance_df_summary_area$method1 <- factor(concordance_df_summary_area$method1,levels = levels(concordance_df_summary_area$method1)[hc$order]) 
  concordance_df_summary_area$method2 <- factor(concordance_df_summary_area$method2,levels = levels(concordance_df_summary_area$method2)[hc$order]) 
  # edges
  edges <- data.frame(x = c(0,0,100,100),
                      xend = c(0,100,100,0),
                      y = c(0,1,1,0),
                      yend = c(1,1,0,0))
  # heatmap
  g_heat <- ggplot(concordance_df_summary_area,aes(x = rank, y = concordance)) +
    #geom_line(size = 1) +
    facet_grid(method1 ~ method2,scales = "free_x",switch = "y") +
    xlab("Rank") + # ylab("Concordance") +
    theme_pubr() + 
    theme(axis.text = element_blank(),
          #axis.text.x = element_text(hjust = 1, angle = 45),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.line.y.right = element_blank(),
          # strip.text = element_text(hjust = 100, vjust = 100),
          # strip.background = element_rect(fill = "gray",linetype = 1,color = "white")) + 
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(0,"cm"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
    #geom_abline(mapping = aes(intercept = 0,slope = 1/nfeatures),color = "red",lty = 2) +
    coord_cartesian(xlim = c(0,100), ylim = c(0,1)) +
    geom_ribbon(aes(ymin = rank/nfeatures, ymax = concordance, fill = area)) +
    geom_segment(concordance_df_summary_area[concordance_df_summary_area$method1 == concordance_df_summary_area$method2,],mapping = aes(x = 0, xend = 0, y = 0, yend = 1, color = "red")) +
    geom_segment(concordance_df_summary_area[concordance_df_summary_area$method1 == concordance_df_summary_area$method2,],mapping = aes(x = 100, xend = 100, y = 1, yend = 0, color = "red")) +
    geom_segment(concordance_df_summary_area[concordance_df_summary_area$method1 == concordance_df_summary_area$method2,],mapping = aes(x = 0, xend = 100, y = 1, yend = 1, color = "red")) +
    geom_segment(concordance_df_summary_area[concordance_df_summary_area$method1 == concordance_df_summary_area$method2,],mapping = aes(x = 100, xend = 0, y = 0, yend = 0, color = "red")) +
    #scale_fill_gradientn(colours = c("red","yellow","turquoise"),limits = c(-0.01,1)) +
    scale_fill_distiller(palette = "RdYlBu",limits = c(-0.01,1)) +
    #scale_color_gradientn(colours = c("red","yellow","turquoise"),limits = c(-0.1,1)) +
    scale_y_continuous(breaks = c(0,0.5,1),position = "right") +
    scale_x_continuous(breaks = c(0,50,100))
  
  g_vertical_dendrogram <- ggplot() + 
    geom_segment(data=dendro_data(hc)$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_label(data=dendro_data(hc)$labels, aes(x=x, y=y, label=label, hjust=1,color=label), nudge_y = 0) +
    coord_flip() + scale_y_reverse(expand = c(0,0,0,0)) + scale_x_reverse() +
    scale_color_manual(values = cols) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          legend.position = "none",
          panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
    ggtitle(label = paste(tech,paste(strsplit(comp,split = "_")[[1]],collapse = " vs ")),
            subtitle = "Concordance heatmap")
  
  g_horizontal_dendrogram <- ggplot() + 
    geom_segment(data=dendro_data(hc)$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
    geom_point(data=dendro_data(hc)$labels, aes(x=x, y=y,color=label),size = 5) +
    scale_y_continuous() +
    #scale_y_reverse(expand=c(2,1)) + scale_x_reverse(expand=c(2,1)) +
    scale_color_manual(values = cols) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          legend.position = "none",
          panel.spacing = unit(0, "lines"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  
  addline_format <- function(x,...){
    gsub(':\\s',':\n',x)
  }
  g_heat_w_legend <- get_legend(ggplot(concordance_df_summary_area,aes(x = rank, y = concordance)) +
                                  facet_grid(method1 ~ method2,scales = "free_x",switch = "y") +
                                  labs(fill = addline_format("Rescaled Area from Rank 1 to 100 between: bisector and concordance")) +
                                  theme_minimal() + 
                                  theme(legend.position = "bottom") +
                                  guides(fill = guide_colorbar(title.position = "top",barwidth = 15)) +
                                  geom_ribbon(aes(ymin = rank/nfeatures, ymax = concordance, fill = area),alpha = 0.8) +
                                  scale_fill_distiller(palette = "RdYlBu",limits = c(0,1)) +
                                  scale_y_continuous(breaks = c(0,0.5,1),position = "right") +
                                  scale_x_continuous(breaks = c(0,50,100)))
  
  a <- plot_grid(plotlist = list(g_horizontal_dendrogram,g_horizontal_dendrogram,g_vertical_dendrogram,g_heat),align = "hv",axis = "lrtb")
  b <- g_heat_w_legend
  c <- g_legend_dendrogram
  return(list(plot = a,legend_heat = b, legend_dendro = c))
}


# div = diversity: high, mid, low
# tech = data type: 16S or WMS
g_AUC50 <- function(conc_df,div,tech){
  conc_df_sub <- conc_df[conc_df$rank == 50 & conc_df$subset == "1vs2",]
  case <- conc_df_sub[conc_df_sub$tech == tech & conc_df_sub$diversity == div,]
  ord <- order(ddply(case,.variables = ~ method1, function(x) median(x[,"concordance"]))$V1)
  ggplot(case,aes(x = method1, y = concordance, color = method1)) +
    geom_boxplot() +
    coord_flip() +
    scale_x_discrete(limits = levels(case$method1)[rev(ord)]) +
    xlab("Method") + ylab("Concordance") +
    ggtitle(label = paste(strsplit(as.character(unique(case$comp)),split = "_")[[1]],collapse = " vs "),
            subtitle = paste(unique(case$tech),"-",unique(case$diversity),"diversity")) +
    theme_minimal() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none",
          panel.spacing = unit(1,"lines")) +
    scale_color_manual(values = cols) +
    scale_y_continuous(limits = c(0,1))
}


