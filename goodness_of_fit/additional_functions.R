library(zinbwave)
library(MAST)
library(metagenomeSeq)
library(EDASeq)
library(edgeR)
library(MGLM)
library(plyr)

### Extraction of ZINB coefs
computeExp <- function(zinbModel){
  (1 - t(getPi(zinbModel))) * t(getMu(zinbModel))
}
computeVar <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  (1 - pi) * mu * (1 + mu*(phi + pi))
}
computeP0 <- function(zinbModel){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi)
}

# Negative Binomial fitting
fitNB_TMM <- function(counts, design){
  normFacts <- edgeR::calcNormFactors(counts, method = "TMM")
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- normFacts
  disp <- estimateDisp(y = dge,design = design,tagwise = TRUE)
  fit <- glmFit(dge$counts,dispersion = disp$tagwise.dispersion, design = design)
  list(fitted = fit$fitted.values, disp = disp$tagwise.dispersion)
}

# Truncated Gaussian Hurdle model fitting
fitHURDLE <- function(ps,design,scale = c("median","default")){
  counts <- as(ps@otu_table@.Data,"matrix")
  if(scale == "median"){
    tpm <- counts*median(colSums(counts))/colSums(counts)
  } 
  if(scale == "default"){ 
    tpm <- counts*1e6/colSums(counts) 
  }
  tpm <- log2(tpm + 1)
  sca <- FromMatrix(tpm, cData=data.frame(sample_data(ps)))
  # Adaptive thresholding from MAST vignette
  # Non necessary if there isn't bimodal density
  # freq_expressed <- 0.05
  # thres <- thresholdSCRNACountMatrix(assay(sca), conditions = 1)
  # assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
  # expressed_genes <- freq(sca) > freq_expressed
  # sca <- sca[expressed_genes,]
  ngeneson <- apply(counts,2,function(x) mean(x>0))
  CD <- colData(sca)
  CD$ngeneson <- ngeneson
  CD$cngeneson <- CD$ngeneson-mean(ngeneson)
  if(ncol(design) > 1)
    CD$grp <- design[,2]
  colData(sca) <- CD
  # hurdle model estimation
  if(ncol(design) > 1){
    hm <- zlm(~ 1 + grp + cngeneson ,method = "bayesglm", sca = sca)
    betaC <- coef(hm,"C")
    fittedC <- tcrossprod(betaC,model.matrix(~ 1 + grp + cngeneson,data = colData(sca)))
    betaD <- coef(hm,"D")
    fittedD <- tcrossprod(betaD,model.matrix(~ 1 + grp + cngeneson,data = colData(sca)))
  } else{
    hm <- zlm(~ 1 + cngeneson ,method = "bayesglm", sca = sca)
    betaC <- coef(hm,"C")
    fittedC <- tcrossprod(betaC,model.matrix(~ 1 + cngeneson,data = colData(sca)))
    betaD <- coef(hm,"D")
    fittedD <- tcrossprod(betaD,model.matrix(~ 1 + cngeneson,data = colData(sca)))
  } 
  mu <- rowMeans(fittedC)
  pi <- rowMeans(invlogit(fittedD))
  # disp <- hm@dispersion[,1]
  # hmY = log1p(rowMeans(2^assay(sca)-1))
  # hmY0 = apply(assay(sca),1,function(row) mean(row==0))
  hmEY = rowMeans(invlogit(fittedD)*fittedC)*log(2)
  hmEY0 = 1-pi
  list(model = hm, EY = hmEY, EY0 = hmEY0)
}

# Zero-inflated Gaussian mixture model fitting
fitZIG <- function(ps,design,scale = c("median","default")){
  ## sample data converted to Annotated Data Frame
  ADF <- AnnotatedDataFrame(data.frame(sample_data(ps)))
  ## design matrix
  MGS <- newMRexperiment(counts = ps@otu_table@.Data, phenoData = ADF)
  MGSp = cumNormStat(MGS)
  MGS <- cumNorm(MGS,MGSp)
  normFactor = normFactors(MGS)
  if(scale == "median"){
    normFactor = log2(normFactor/median(normFactor) + 1)
  }
  else { normFactor = log2(normFactor/1000 + 1) }
  desMat <- cbind(design,normFactor = normFactor)
  zig <- fitZig(MGS,desMat,control = zigControl(maxit = 1000),useCSSoffset = FALSE)
  # mmCount <- cbind(1L, log2(normFactors(MGS)/1000 + 1)) if useCSSoffset = TRUE
  countMu <- tcrossprod(zig@fit$coefficients, desMat)
  zigEY <- rowMeans(countMu)*log(2)
  zigEY0 <- rowMeans(zig@z)
  list(model = zig, EY = zigEY, EY0 = zigEY0)
}

# Dirichlet-Multinomial fitting
fitDirMult <- function(ps, design){
  # Counts
  Y <- t(ps@otu_table@.Data)
  # design <- model.matrix(~ 1, data.frame(ps@sam_data))

  # N
  ls <- rowSums(Y)
  
  # desMat <- cbind(design,log_g)
  if(length(colnames(design)[-1])==0){
    desFormula = as.formula("Y ~ 1")
  } else desFormula = as.formula(paste0("Y ~ 1 + ",paste0(colnames(design)[-1],collapse = "+")))
  
  dmFit <- MGLMreg(desFormula, dist="DM",display = TRUE)
  fitted_values <- dmFit@fitted * ls
  alpha_i <- t(exp(t(dmFit@coefficients) %*% t(design)))
  alpha_0 <- rowSums(alpha_i)
  
  EY <- log1p(colMeans(ls*alpha_i/alpha_0))
  EY0 <- colMeans(beta(alpha_i,ls+alpha_0-alpha_i)/beta(alpha_i,alpha_0-alpha_i))
  list(EY = EY, EY0 = EY0)
}

# ZINB, NB, HURDLE, ZIG, DM model estimation

fit_models <- function(ps, design, epsilon = 1e10){
  if(!taxa_are_rows(otu_table(ps)))
    otu_table(ps) <- t(otu_table(ps))
  cat(unique(ps@sam_data$HMP_BODY_SUBSITE))
  ### ZINB
  cat("\n Starting ZINB \n")
  zinb <- zinbFit(ps@otu_table@.Data,X = design,commondispersion = TRUE,epsilon = epsilon, verbose = TRUE, K = 0, BPPARAM = SerialParam())
  zinbEY = log1p(rowMeans(computeExp(zinb)))
  zinbEY0 = rowMeans(computeP0(zinb))
  zinb_values = list(model = zinb, EY = zinbEY, EY0 = zinbEY0)
  cat("\n ZINB done \n")
  ### NB
  cat("\n Starting NB \n")
  nb <- fitNB_TMM(ps@otu_table@.Data,design = design)
  nbEY = log1p(rowMeans(nb$fitted))
  nbEY0 = rowMeans((1 + nb$fitted * nb$disp)^(-1/nb$disp))
  nbloglik <- rowSums(dnbinom(x = ps@otu_table@.Data,size=1/nb$disp,mu=rowMeans(nb$fitted),log = TRUE))
  # parameters: 1 dispersion and 1 mean for each feature
  nb_values = list(model = nb, EY = nbEY, EY0 = nbEY0)
  cat("\n NB done \n")
  ### HURDLE MODELS
  cat("\n Starting Hurdle with TPM \n")
  hm_default_values <- fitHURDLE(ps = ps, design = design, scale = "default")
  cat("\n Hurdle with TPM done \n")
  cat("\n Starting Hurdle with TPmedian \n")
  hm_median_values <- fitHURDLE(ps = ps, design = design, scale = "median")
  cat("\n Hurdle with TPmedian done \n")
  ### ZIG
  cat("\n Starting default ZIG \n")
  zig_default_values <- fitZIG(ps = ps, design = design, scale = "default")
  cat("\n Standard ZIG done \n")
  cat("\n Starting ZIG with median scaling factor \n")
  zig_median_values <- fitZIG(ps = ps, design = design, scale = "median")
  cat("\n ZIG with median scaling factor done \n")
  ### DM
  cat("\n Starting Dirichlet-Multinomial \n")
  DM_values <- fitDirMult(ps = ps,design = design)
  cat("\n Dirichlet Multinomial \n")
  
  ### OBSERVED
  Y = log1p(rowMeans(ps@otu_table@.Data))
  Y0 = rowMeans(ps@otu_table@.Data == 0)
  real_values <- list(Y = Y, Y0 = Y0)
  return(list(ZINB = zinb_values, 
              NB = nb_values, 
              HURDLE_default = hm_default_values, 
              HURDLE_median = hm_median_values, 
              ZIG_default = zig_default_values, 
              ZIG_median = zig_median_values,
              DM = DM_values,
              OBSERVED = real_values))
  cat("\n Dataset end \n")
}

extract_values <- function(site,distributions = c("ZINB","NB","HURDLE_default","HURDLE_median","ZIG_default","ZIG_median","DM"), varnames = c("EY","EY0")){
  cbind(ldply(lapply(site[distributions],function(dist){
    as.data.frame(dist[varnames])
  }),.id = "Models"),
  ldply(lapply(site["OBSERVED"],function(dist){
    observed_df <- as.data.frame(dist)
    return(observed_df)
  })))
}

# From model list to data.frame
model_to_data.frame <- function(model_list){
  df_list <- lapply(model_list,extract_values)
  df <- ldply(df_list)
  colnames(df) <- c("Models","EY","EY0","Subsites","Y","Y0")
  df$Subsites <- rep(names(df_list),sapply(df_list,nrow))
  # Mean differences for average and zero probability
  df$mean_diff <- df$EY - df$Y
  df$zero_prob <- df$EY0 - df$Y0
  return(df)
}

# GOF indexes
# Root mean square errors
compute_RMSE <- function(df){
  ddply(df,.variables = ~ Models + Subsites, function(x){
    MD <- sqrt(mean(x$mean_diff^2,na.rm = TRUE))
    ZPD <- sqrt(mean(x$zero_prob^2,na.rm = TRUE))
    return(data.frame("MD" = MD, "ZPD" = ZPD))
  })
}