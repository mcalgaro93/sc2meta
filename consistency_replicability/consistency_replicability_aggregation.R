### To aggregate the single output of the job array we use a list

# 16S
subgingival_supragingival <- list()
for(i in 1:100){
  subgingival_supragingival[[i]] <- readRDS(paste0("../data/16Ssubsets_replicability_DA_subgingival_supragingival",i,".RDS"))
}
names(subgingival_supragingival) <- paste0("Comparison",1:100)

gingiva_mucosa <- list()
for(i in 1:100){
  gingiva_mucosa[[i]] <- readRDS(paste0("../data/16Ssubsets_replicability_DA_gingiva_mucosa",i,".RDS"))
}
names(gingiva_mucosa) <- paste0("Comparison",1:100)

tonguedorsum_stool <- list()
sequence = c(1:31,101,33:39,102,41:59,103,61:67,104,69:72,105,74:100)
for(i in 1:100){
  tonguedorsum_stool[[i]] <- readRDS(paste0("../data/16Ssubsets_replicability_DA_tonguedorsum_stool",sequence[i],".RDS"))
}
names(tonguedorsum_stool) <- paste0("Comparison",1:100)

subsets_replicability_DA_16S <- list(subgingival_supragingival = subgingival_supragingival,
                                     gingiva_mucosa = gingiva_mucosa,
                                     tonguedorsum_stool = tonguedorsum_stool)

saveRDS(subsets_replicability_DA_16S,file = "../data/16Ssubsets_replicability_DA.RDS")

# WMS
CRC_control <- list()
for(i in 1:100){
  CRC_control[[i]] <- readRDS(paste0("../data/WMSsubsets_replicability_DA_CRC_control",i,".RDS"))
}
names(CRC_control) <- paste0("Comparison",1:100)

schizophrenia_control <- list()
for(i in 1:100){
  schizophrenia_control[[i]] <- readRDS(paste0("../data/WMSsubsets_replicability_DA_schizophrenia_control",i,".RDS"))
}
names(schizophrenia_control) <- paste0("Comparison",1:100)

tonguedorsum_stool <- list()
for(i in 1:100){
  tonguedorsum_stool[[i]] <- readRDS(paste0("../data/WMSsubsets_replicability_DA_tonguedorsum_stool",i,".RDS"))
}
names(tonguedorsum_stool) <- paste0("Comparison",1:100)

WMSsubsets_replicability_DA <- list(CRC_control = CRC_control,
                                    schizophrenia_control = schizophrenia_control,
                                    tonguedorsum_stool = tonguedorsum_stool)

saveRDS(WMSsubsets_replicability_DA,file = "../data/WMSsubsets_replicability_DA.RDS")

################ songbird ##############################

songbird_data_collector <- function(filename,sequence = 1:100){
  comparisons <- list()
  for(i in 1:100){
    subsets <- list()
    for(j in 1:2){
      subset = read.table(file = paste0(filename,sequence[i],"_Subset",j,"/differentials.tsv"),
                          header = TRUE,sep = "\t",row.names = 1)
      subsets[[j]] <- subset
      
    }
    names(subsets) <- c("Subset1","Subset2")
    comparisons[[i]] <- subsets
  }
  names(comparisons) <- paste0("Comparison",1:100)
  return(comparisons)
}

# 16S
filename = "../consistency_replicability/songbird/16S_subgingival_supragingival_Comparison"
subgingival_supragingival <- songbird_data_collector(filename)

filename = "../consistency_replicability/songbird/16S_gingiva_mucosa_Comparison"
gingiva_mucosa <- songbird_data_collector(filename)

filename = "../consistency_replicability/songbird/16S_tonguedorsum_stool_Comparison"
# Since some corncob method failed in some sims, this is the matched sequence to obtain 100 sims
sequence = c(1:31,101,33:39,102,41:59,103,61:67,104,69:72,105,74:100)
tonguedorsum_stool <- songbird_data_collector(filename,sequence = sequence)

subsets_replicability_DA_16S_songbird <- list(subgingival_supragingival = subgingival_supragingival,
                                              gingiva_mucosa = gingiva_mucosa,
                                              tonguedorsum_stool = tonguedorsum_stool)

saveRDS(subsets_replicability_DA_16S_songbird,file = "../data/16Ssubsets_replicability_DA_songbird.RDS")

# WMS
filename = "../consistency_replicability/songbird/WMS_CRC_control_Comparison"
CRC_control <- songbird_data_collector(filename)

filename = "../consistency_replicability/songbird/WMS_schizophrenia_control_Comparison"
schizophrenia_control <- songbird_data_collector(filename)

filename = "../consistency_replicability/songbird/WMS_tonguedorsum_stool_Comparison"
tonguedorsum_stool <- songbird_data_collector(filename)

WMSsubsets_replicability_DA_songbird <- list(CRC_control = CRC_control,
                                             schizophrenia_control = schizophrenia_control,
                                             tonguedorsum_stool = tonguedorsum_stool)

saveRDS(WMSsubsets_replicability_DA_songbird,file = "../data/WMSsubsets_replicability_DA_songbird.RDS")

#############################################
#                                           #
## put together songbird and other methods ##
#                                           #
#############################################

join_all_methods <- function(method_list,songbird_list){
  
  for(dataset in names(method_list)){
    # print(paste0("dataset: ",dataset))
    data <- method_list[[dataset]]
    data_songbird <- songbird_list[[dataset]]
    comparisons <- list()
    for(i in 1:100){
      # print(paste0("comparison: ",i))
      comparison <- data[[i]]
      comparison_songbird <- data_songbird[[i]]
      subsets <- list()
      for(j in 1:2){
        # print(paste0("subset: ",j))
        subset <- comparison[[j]]
        subset_songbird <- comparison_songbird[[j]]
        subset$songbird <- subset_songbird
        subsets[[j]] <- subset
        # print(paste(names(subset)))
      }
      names(subsets) <- c("Subset1","Subset2")
      comparisons[[i]] <- subsets
    }
    names(comparisons) <- paste0("Comparison",1:100)
    method_list[[dataset]] <- comparisons
  }
  return(method_list)
}

WMSsubsets_replicability_DA <- readRDS(file = "../data/WMSsubsets_replicability_DA.RDS")
subsets_replicability_DA_16S <- readRDS(file = "../data/16Ssubsets_replicability_DA.RDS")
WMSsubsets_replicability_DA_songbird <- readRDS(file = "../data/WMSsubsets_replicability_DA_songbird.RDS")
subsets_replicability_DA_16S_songbird <- readRDS(file = "../data/16Ssubsets_replicability_DA_songbird.RDS")

WMSsubsets_replicability_DA_total <- join_all_methods(method_list = WMSsubsets_replicability_DA,
                                                      songbird_list = WMSsubsets_replicability_DA_songbird)
subsets_replicability_DA_total_16S <- join_all_methods(method_list = subsets_replicability_DA_16S,
                                                       songbird_list = subsets_replicability_DA_16S_songbird)

saveRDS(WMSsubsets_replicability_DA_total,file = "../data/WMSsubsets_replicability_DA_total.RDS")
saveRDS(subsets_replicability_DA_total_16S,file = "../data/16Ssubsets_replicability_DA_total.RDS")