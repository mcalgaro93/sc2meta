## Between Dataset Concordance and Within Dataset Concordance 

source("./additional_functions.R")

### WMS

# high diversity
WMSsubsets_replicability_DA <- readRDS(file = "../data/WMSsubsets_replicability_DA_total.RDS")
ps_fitted_list <- WMSsubsets_replicability_DA$tonguedorsum_stool
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/WMSconcordance_tonguedorsum_stool.RData")

# mid diversity
ps_fitted_list <- WMSsubsets_replicability_DA$schizophrenia_control
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/WMSconcordance_schizophrenia_control.RData")

# low diversity
ps_fitted_list <- WMSsubsets_replicability_DA$CRC_control
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/WMSconcordance_CRC_control.RData")