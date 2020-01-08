library(phyloseq)
library(plyr)
library(zinbwave)
library(HMP16SData)
library(curatedMetagenomicData)

### STOOL
if(!file.exists("./data/Stool_16S_WMS_models.RData")){
  load("./data/Stool_16S_WMS.RData")
  Stool_16S_WMS_models <- lapply(Stool_16S_WMS,function(ps){
    
    zinbFit(ps@otu_table@.Data,
            X = model.matrix(~ 1,data = data.frame(ps@sam_data)),
            commondispersion = TRUE,
            epsilon = 1e10,
            verbose = TRUE, 
            K = 0, 
            BPPARAM = SerialParam())
  })
  save(Stool_16S_WMS_models,file = "./data/Stool_16S_WMS_models.RData")
} else load(file = "./data/Stool_16S_WMS_models.RData")

### TONGUE DORSUM
if(!file.exists("./data/TD_16S_WMS_models.RData")){
  # WMS Data
  # You can choose if to use the ready models and extract ZINB ...
  load("./data/HMPWMSData_models.RData")
  load("./data/HMPWMSData_list.RData")
  TD_16S_WMS = list("TD_16S" = HMP16SData_list$`Tongue Dorsum`,"TD_WMS" = HMPWMSData_list$`Tongue Dorsum`)
  TD_WMS <- HMPWMSData_models$`Tongue Dorsum`$ZINB$model
  # 16S Data
  # ... or to estimate ZINB from data.
  load("./data/HMP16SData_list.RData")
  TD_16S <- zinbFit(HMP16SData_list$`Tongue Dorsum`@otu_table@.Data,
              X = model.matrix(~ 1,data = data.frame(HMP16SData_list$`Tongue Dorsum`@sam_data)),
              commondispersion = TRUE,
              epsilon = 1e10,
              verbose = TRUE, 
              K = 0, 
              BPPARAM = SerialParam())
  TD_16S_WMS_models <- list(TD_16S = TD_16S, TD_WMS = TD_WMS)
  save(TD_16S_WMS,file = "./data/TD_16S_WMS.RData")
  save(TD_16S_WMS_models,file = "./data/TD_16S_WMS_models.RData")
} else load("./data/TD_16S_WMS_models.RData")

### Brito Il 2016
# WMS Data
if(!file.exists("./data/BritoIL_Stool_Oral_models.RData")){
  load("./data/AllWMSData_list.RData")
  AllWMSData_list$BritoIL_2016_oralcavity_control
  
  BritoIL_Stool <- subset_samples(AllWMSData_list$BritoIL_2016_stool_control, gender == "female" & age_category == "adult")
  BritoIL_Oral <- subset_samples(AllWMSData_list$BritoIL_2016_oralcavity_control, gender == "female" & age_category == "adult")
  # Pruning
  BritoIL_Stool <- prune_samples(sample_sums(BritoIL_Stool) >= 10^6, BritoIL_Stool)
  BritoIL_Oral <- prune_samples(sample_sums(BritoIL_Oral) >= 10^6, BritoIL_Oral)
  # Filtering rarest taxa
  BritoIL_Stool <- filter_taxa(BritoIL_Stool,function(x) sum(x>10)>1,1)
  BritoIL_Oral <- filter_taxa(BritoIL_Oral,function(x) sum(x>10)>1,1)
  
  BritoIL_Stool_Oral <- list(BritoIL_Stool = BritoIL_Stool,BritoIL_Oral = BritoIL_Oral)
  save(BritoIL_Stool_Oral,file = "./data/BritoIL_Stool_Oral.RData")
  BritoIL_Stool_Oral_models <- lapply(BritoIL_Stool_Oral,function(ps){
    zinbFit(ps@otu_table@.Data,
            X = model.matrix(~ 1,data = data.frame(ps@sam_data)),
            commondispersion = TRUE,
            epsilon = 1e10,
            verbose = TRUE, 
            K = 0, 
            BPPARAM = SerialParam())
  })
  save(BritoIL_Stool_Oral_models,file = "./data/BritoIL_Stool_Oral_models.RData")
} else load(".data/BritoIL_Stool_Oral_models.RData")


