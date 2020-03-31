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
for(i in 1:100){
  tonguedorsum_stool[[i]] <- readRDS(paste0("../data/16Ssubsets_replicability_DA_tonguedorsum_stool",i,".RDS"))
}
names(tonguedorsum_stool) <- paste0("Comparison",1:100)

subsets_replicability_DA_16S <- list(subgingival_supragingival = subgingival_supragingival,
                                     gingiva_mucosa = gingiva_mucosa,
                                     tonguedorsum_stool = tonguedorsum_stool)

saveRDS(subsets_replicability_DA_16S,file = "../data/16Ssubsets_replicability_DA.RData")

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

saveRDS(WMSsubsets_replicability_DA,file = "../data/WMSsubsets_replicability_DA.RData")