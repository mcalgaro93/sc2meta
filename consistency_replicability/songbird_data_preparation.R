library(phyloseq)

songbird_data_creator <- function(ps_list, filename, offset = 0){
  for(i in 1:length(ps_list)){
    comparison <- ps_list[[i]]
    for(j in 1:length(comparison)){
      ps <- comparison[[j]]
      sample_names(ps) <- paste0("S_",sample_names(ps))
      write.table(x = "#OTU ID\t",
                  file = paste0(filename,i+offset,"_Subset",j,"_otutable.tsv"),quote = FALSE,row.names = FALSE, col.names = FALSE,eol = "")
      write.table(x = ps@otu_table, append = TRUE,
                  file = paste0(filename,i+offset,"_Subset",j,"_otutable.tsv"),
                  quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      write.table(x = "sampleID\t",
                  file = paste0(filename,i+offset,"_Subset",j,"_samdata.tsv"),quote = FALSE,row.names = FALSE, col.names = FALSE,eol = "")
      write.table(x = ps@sam_data, append = TRUE,
                  file = paste0(filename,i+offset,"_Subset",j,"_samdata.tsv"),
                  quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      
    }
  }
}

load("../data/16Ssubsets_replicability.RData")
ps_list <- subsets_consistency_replicability_16S$subgingival_supragingival
filename <- "./songbird/16S_subgingival_supragingival_Comparison"
songbird_data_creator(ps_list,filename)

ps_list <- subsets_consistency_replicability_16S$gingiva_mucosa
filename <- "./songbird/16S_gingiva_mucosa_Comparison"
songbird_data_creator(ps_list,filename)

ps_list <- subsets_consistency_replicability_16S$tonguedorsum_stool
filename <- "./songbird/16S_tonguedorsum_stool_Comparison"
songbird_data_creator(ps_list,filename)

###########################################################################à
load("../data/WMSsubsets_replicability.RData")
ps_list <- subsets_consistency_replicability_WMS$CRC_control
filename <- "./songbird/WMS_CRC_control_Comparison"
songbird_data_creator(ps_list,filename)

ps_list <- subsets_consistency_replicability_WMS$schizophrenia_control
filename <- "./songbird/WMS_schizophrenia_control_Comparison"
songbird_data_creator(ps_list,filename)

ps_list <- subsets_consistency_replicability_WMS$tonguedorsum_stool
filename <- "./songbird/WMS_tonguedorsum_stool_Comparison"
songbird_data_creator(ps_list,filename)

######### for simulations over 100 (in those cases where corncob failed) ######
load("../data/WMS16Ssubsets_replicability_tonguedorsum_stool.RData")
ps_list <- tonguedorsum_stoolWMS16S$tonguedorsum_stool16S
filename <- "./songbird/16S_tonguedorsum_stool_Comparison"
songbird_data_creator(ps_list,filename,offset = 100)

ps_list <- tonguedorsum_stoolWMS16S$tonguedorsum_stoolWMS
filename <- "./songbird/WMS_tonguedorsum_stool_Comparison"
songbird_data_creator(ps_list,filename,offset = 100)