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