## Between Dataset Concordance and Within Dataset Concordance 

source("./additional_functions.R")

# high diversity
subsets_replicability_DA_16S <- readRDS(file = "../data/16Ssubsets_replicability_DA_total.RDS")
ps_fitted_list <- subsets_replicability_DA_16S$tonguedorsum_stool
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/16Sconcordance_tonguedorsum_stool.RData")

# mid diversity
ps_fitted_list <- subsets_replicability_DA_16S$gingiva_mucosa
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/16Sconcordance_gingiva_mucosa.RData")

# low diversity
ps_fitted_list <- subsets_replicability_DA_16S$subgingival_supragingival
  conc_df <- compute_concordance(ps_fitted_list = ps_fitted_list)
  save(conc_df,file = "../data/16Sconcordance_subgingival_supragingival.RData")