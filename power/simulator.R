library(phyloseq)
library(zinbwave)
library(edgeR)
library(EDASeq)
library(ggplot2)
library(gridExtra)
library(BiocParallel)
library(digest)

## Load only one of these:

### Stool
load(file = "./data/Stool_16S_WMS.RData")
load(file = "./data/Stool_16S_WMS_models.RData")
data_list = Stool_16S_WMS
models = Stool_16S_WMS_models

### Tongue Dorsum
load(file = "./data/TD_16S_WMS.RData")
load(file = "./data/TD_16S_WMS_models.RData")
data_list = TD_16S_WMS
models = TD_16S_WMS_models

### BritoIL Stool and Oral Cavity
load(file = "./data/BritoIL_Stool_Oral.RData")
load(file = "./data/BritoIL_Stool_Oral_models.RData")
data_list = BritoIL_Stool_Oral
models = BritoIL_Stool_Oral_models


addFoldChange = function(rhos, fc, H1frac = TPR, compensate = FALSE) {
  if (fc == 1) {
    return(rhos)
  }
  nTaxa = length(rhos)
  if (compensate) {
    nOTUsUp = round(nTaxa * H1frac * (1/(fc + 1)))  #Upregulated taxa
    nOTUsDown = round(nTaxa * H1frac - nOTUsUp)  #Downregulated taxa
    cond=TRUE 
    while(cond){
      OTUids = sample(names(rhos), nOTUsUp + nOTUsDown, replace = FALSE)
      OTUidUps = OTUids[1:nOTUsUp]
      OTUidDowns = OTUids[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
      a = sum(rhos[OTUidUps])
      b = sum(rhos[OTUidDowns])
      cond = (b/a) + 1 <= fc
    }
    rhos[OTUidUps] = rhos[OTUidUps] * fc  # Add fold change up
    rhos[OTUidDowns] = rhos[OTUidDowns] * ((a/b) * (1-fc) + 1)  #And compensate the downs. This way the average FC is 5 in both directions and the TN taxa are really left untouched
  } else {
    nOTUsUp = nOTUsDown = round(nTaxa * H1frac/2)  #Upregulated taxa and Downregulated taxa
    OTUids = sample(names(rhos), nOTUsUp + nOTUsDown, replace = FALSE)
    OTUidUps = OTUids[1:nOTUsUp]
    OTUidDowns = OTUids[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
    rhos[OTUidUps] = rhos[OTUidUps] * fc  # Add fold change up
    rhos[OTUidDowns] = rhos[OTUidDowns] / fc # fold change down
  }
  indTPup <- names(rhos) %in% OTUidUps
  newTaxaNamesUp <- paste0(names(rhos)[indTPup], "-TPup")
  indTPdown <- names(rhos) %in% OTUidDowns
  newTaxaNamesDown <- paste0(names(rhos)[indTPdown], "-TPdown")
  names(rhos)[indTPup] <- newTaxaNamesUp
  names(rhos)[indTPdown] <- newTaxaNamesDown
  rhos/sum(rhos)  #Renormalize. 
}

# # Trim by prevalence and total OTU reads
simpleTrimGen <- function(otuTab, minReads = 10, minPrev = 1) {
  # `prevalence` is the fraction of samples in which an OTU is observed at
  # least `minReads` times.
  # prevalence <- rowMeans(otuTab > minReads)
  absolute <- rowSums(otuTab > minReads)
  # cv_index <- apply(otuTab,1,function(row) stats::sd(row)/mean(row))
  ## Will only keep OTUs that appear in more than minPrev samples 
  # indOTUs2Keep <- (prevalence > minPrev)
  indOTUs2Keep <- (absolute > minPrev)
  return(otuTab[indOTUs2Keep, ])
}  # END - function: simpleTrim general

computeExactWeights <- function (model, x) 
{
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  zinbwg
}

dataset <- names(data_list)
distribution <- c("NB","ZINB")
simulation = 1:50
sampleSize = c(10,20,40)
TPR = c(0.1, 0.5)
foldEffect <- c(2, 5)
compensation <- c("yes","no")
sparsityEffect <- c(0, 0.05, 0.15)

simulation_flow <- data.frame(expand.grid(dataset = dataset,
                                          distribution = distribution,
                                          simulation = simulation,
                                          sampleSize = sampleSize,
                                          TPR = TPR,
                                          foldEffect = foldEffect, 
                                          compensation = compensation,
                                          sparsityEffect = sparsityEffect))
# Removing senseless simulations: 
# 1. nb generated dataset with sparsity_effect!=0
simulation_flow <- simulation_flow[-which(simulation_flow$sparsityEffect != 0 & simulation_flow$distribution == "NB"),]
rownames(simulation_flow) <- 1:nrow(simulation_flow)
# Unique seed for each mix of variables
for(i in 1:nrow(simulation_flow))
{
  simulation_flow$seed[i] <- strtoi(paste0("0x",substr(digest::digest(simulation_flow[i,1:ncol(simulation_flow)],algo = "sha256"),start = 1,stop = 7)))
}

### Stool
save(simulation_flow,file = "./data/Stool_16S_WMS_simulation_flow.RData")
### Tongue Dorsum
save(simulation_flow,file = "./data/TD_16S_WMS_simulation_flow.RData")
### BritoIL Stool and Oral Cavity
save(simulation_flow,file = "./data/BritoIL_Stool_Oral_simulation_flow.RData")

# Generating simulations
sims <- apply(simulation_flow,1,function(sim){
  simName <- paste(colnames(simulation_flow)[1:ncol(simulation_flow)],sim[1:ncol(simulation_flow)],sep = ":",collapse = "_")
  cat(simName,"\n")
  set.seed(sim[ncol(simulation_flow)])
  # taxa-specific mean extraction for NB component of the mixture model ZINB
  true_model <- models[[sim[1]]]
  tot_zinb_mu <- sum(exp(true_model@beta_mu))
  # link function is log
  zinb_mu_rel <- as.numeric(exp(true_model@beta_mu))/tot_zinb_mu 
  zinb_mu <- as.numeric(true_model@beta_mu)
  
  # Simulation with original parameters for experimental group 1
  nb_zinb_sim_data1 <- zinbSim(true_model) 
  
  # OTU naming
  OTU_names <- paste("OTU_",1:nrow(nb_zinb_sim_data1$counts),sep = "") 
  names(zinb_mu_rel) <- OTU_names
  
  # Sampling number sample_size samples from generated data
  # Index for group 1 and 2
  gr1 <- sample(x = 1:nsamples(data_list[[sim[1]]]),replace = FALSE,size = as.numeric(sim[4])) 
  gr2 <- sample(x = 1:nsamples(data_list[[sim[1]]]),replace = FALSE,size = as.numeric(sim[4]))
  zinb_mu_rel_fc <- addFoldChange(rhos = zinb_mu_rel, fc = as.numeric(sim[6]), H1frac = as.numeric(sim[5]), compensate = (sim[7]=="yes"))
  zinb_mu_fc <- zinb_mu_rel_fc * tot_zinb_mu
  
  # True weights for glm testing
  true_weights_NB_data1 <- computeExactWeights(model = true_model,x = nb_zinb_sim_data1$dataNB)
  true_weights_ZINB_data1 <- computeExactWeights(model = true_model,x = nb_zinb_sim_data1$counts)
  
  ### Plots mu
  # par(mfrow = c(2,2))
  # plot(zinb_mu,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Original Up",ylim = c(-5,5))
  # plot(zinb_mu,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Original Down",ylim = c(-5,5))
  # plot(log(zinb_mu_fc),col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Added fold effect Up",ylim = c(-5,5))
  # plot(log(zinb_mu_fc),col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Added fold effect Down",ylim = c(-5,5))
  
  true_model@beta_mu <- matrix(log(zinb_mu_fc),nrow = 1)
  
  if(sim[2]=="ZINB"){
    # Pi parameter for zero probability in ZINB
    # logistic link
    zinb_pi <- as.numeric(exp(true_model@beta_pi)/(1+exp(true_model@beta_pi))) 
    ### Plots pi
    # par(mfrow = c(2,2))
    # plot(zinb_pi,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Original Up")
    # plot(zinb_pi,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Original Down")
    
    zinb_pi[grepl(x = names(zinb_mu_rel_fc), pattern = "up")] <- unlist(sapply(zinb_pi[grepl(x = names(zinb_mu_rel_fc), pattern = "up")],function(x) max(x-as.numeric(sim[8]),1e-08)))
    zinb_pi[grepl(x = names(zinb_mu_rel_fc), pattern = "down")] <- unlist(sapply(zinb_pi[grepl(x = names(zinb_mu_rel_fc), pattern = "down")],function(x) min(x+as.numeric(sim[8]),1-1e-8)))
    
    # plot(zinb_pi,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "up"))+1,main = "Up after decreasing sparsity")
    # plot(zinb_pi,col = as.numeric(grepl(x = names(zinb_mu_fc),pattern = "down"))+1,main = "Down after increasing sparsity")
    
    # Back to original scale
    true_model@beta_pi <- matrix(log(zinb_pi/(1-zinb_pi)),nrow = 1)
    # Simulation with changed parameters
    nb_zinb_sim_data2 <- zinbSim(true_model)
    true_weights_data2 <- computeExactWeights(model = true_model,x = nb_zinb_sim_data2$counts)
    true_weights <- cbind(true_weights_ZINB_data1[,gr1],true_weights_data2[,gr2])
    sim_counts <- cbind(nb_zinb_sim_data1$counts[,gr1],
                        nb_zinb_sim_data2$counts[,gr2])
  } else {
    # Simulation with changed parameters for experimental group 2
    nb_zinb_sim_data2 <- zinbSim(true_model) 
    true_weights_data2 <- computeExactWeights(model = true_model,x = nb_zinb_sim_data2$dataNB)
    true_weights <- cbind(true_weights_NB_data1[,gr1],true_weights_data2[,gr2])
    sim_counts <- cbind(nb_zinb_sim_data1$dataNB[,gr1],
                        nb_zinb_sim_data2$dataNB[,gr2])
  }
  
  # Information about True Positives
  OTU_names <- names(zinb_mu_fc) 
  # Sample names
  sample_names <- c(paste("Sample_",1:as.numeric(sim[4]),"_grp1",sep=""),paste("Sample_",(as.numeric(sim[4])+1):(2*as.numeric(sim[4])),"_grp2",sep=""))
  colnames(sim_counts) <- colnames(true_weights) <- sample_names
  rownames(sim_counts) <- rownames(true_weights) <- OTU_names
  # Trim too rare OTUs
  sim_counts_filtered <- simpleTrimGen(sim_counts)
  true_weights_filtered <- true_weights[rownames(sim_counts_filtered),]
  obj <- list(ps = phyloseq(otu_table(sim_counts_filtered,taxa_are_rows = TRUE),sample_data(data.frame(grp = rep(c("grp1","grp2"),each = as.numeric(sim[4])),row.names = sample_names))),
              name = simName,
              true_weights = true_weights_filtered)
  return(obj)
})

for(i in 1:length(sims)){ sims[[i]]$number <- i }
names(sims) <- lapply(sims,function(sim) sim$name)

## Depending on which dataset you chose at the beginning
### Stool
save(sims,file = "./data/Stool_16S_WMS_simulations.RData")
### Tongue Dorsum
save(sims,file = "./data/TD_16S_WMS_simulations.RData")
### BritoIL Stool and Oral Cavity
save(sims,file = "./data/BritoIL_Stool_Oral_simulations.RData")