# Many functions of this script are adapted from the work of:
# A broken promise: microbiome differential abundance methods do not control the false discovery rate.
# The original code is available at https://users.ugent.be/~shawinke/ABrokenPromise/index.html

# It is mandatory to use flexmix v2.3-13 and scde 1.99.1.
# https://cran.r-project.org/src/contrib/Archive/flexmix/
# https://github.com/hms-dbmi/scde/releases
# Seurat v2.3.4
# https://github.com/satijalab/seurat/releases

pkgs <- c("edgeR", 
          "limma", 
          "DESeq2",
          "metagenomeSeq", 
          "phyloseq", 
          "plyr", 
          "reshape2",
          "ROCR",
          "samr",
          "zinbwave",
          "BiocParallel",
          "AUC",
          "genefilter",
          "MAST",
          "scde", # It is important to use flexmix v2.3-13 and scde 1.99.1
          "Seurat",
          "crayon",
          "ALDEx2",
          "corncob",
          "selbal")
for(i in pkgs) { library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE) }

register(SerialParam())
# If you work on a MRAN environment
#setMKLthreads(1)

### Total sum scaling, also known as proportion normalization, dividing by library sizes
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data(physeq) <- aux
  physeq
}# END - function: normNone

### edgeR normalisations: TMM and RLE
normEdgeR <- function(physeq, method = c('TMM', 'RLE', 'upperquartile'))
{
  # require(edgeR)
  otuTab <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq))
  {
    otuTab <- t(otuTab)
  } else {}
  
  if (method == "upperquartile")
  {
    scaledCounts <- t(otuTab) / colSums(otuTab)
    tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
      quantile(x[x != 0], probs = .75))
    normFacts <- tmpNF/exp(mean(log(tmpNF)))
    method <- "UQ"
  } else {
    normFacts <- edgeR:::calcNormFactors(otuTab, method = method)
  }# END - ifelse: upperquartile only of non-zero counts
  #VERY IMPORTANT: multiply by library sizes and renormalize. edgeR calculates scaling factors, which still have to be multiplied by library sizes to get to the size factors of effective sequencing depth, i.e. robust estimates of the library sizes
  #normFacts = normFacts*sample_sums(physeq)
  #normFacts = normFacts/exp(mean(log(normFacts)))
  if (all(is.na(normFacts))) #Resort to proportion normalization in case of failure for all samples
  {
    normFacts = sample_sums(physeq)
  }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normEdgeR

### function that apply different normalisations and build *DESeqDataSet* object
### for DESeq2 analysis
normDESeq2 <- function(physeq, whichOTUs = NULL, method = c("poscounts","ratio"))
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else {}# END - if: whichOTUs
  
  #   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "poscounts")
  {
    obj <- phyloseq_to_deseq2(physeq,design = ~grp)
    normFacts <- sizeFactors(DESeq2::estimateSizeFactors(obj,type = "poscounts"))
  } else {
    otuTab <- as(otu_table(physeq), "matrix")
    if (any(otuTab == 0))
    {
      otuTab <- otuTab + 1L
    } else {}
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  }
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normDESeq2

### Cumulative Sum Scaling from *metagenomeSeq*
normCSS <- function(physeq, geoMean = FALSE, rel = 0.1)
{
  # require(metagenomeSeq)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  otuTab <- as(otu_table(physeq), "matrix")
  
  aux <- newMRexperiment(counts = otuTab)
  normFacts <- metagenomeSeq::calcNormFactors(
    obj = aux, p = cumNormStatFast(aux, rel = rel))
  normFacts <- drop(as.matrix(normFacts))
  normFacts = log2(normFacts/median(normFacts) + 1) # otherwise log2(normfactor/1000 + 1) is computed, remember to useCSSoffset = FALSE
  if (geoMean)
  {
    normFacts <- normFacts / exp(mean(log(normFacts[normFacts > 0]), na.rm = TRUE))
  } else {}
  
  #   physeq@otu_table@.Data <- countsMat * normFac
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- "NF.CSS"
  physeq@sam_data@names <- aux
  return(physeq)
}

### Perform EdgeR, robust version for overdispersion estimation.
### edgeR_QLFTest_robust_3.6.4
#   function (counts, group, design = NULL, prior.df = 10) 
edgeR_robust <- function(physeq, design = as.formula("~ grp"), prior.df = 10, normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dgeW <- estimateGLMRobustDisp(y = dge, design, prior.df = prior.df, maxit = 10)
  glmFit <- glmQLFit(y = dgeW, dispersion = dgeW$tagwise.dispersion, robust = TRUE,
                     design = design)
  glmRes <- glmQLFTest(glmFit, coef = 2)
  pval <- glmRes$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmRes$table$logFC, "logCPM" = glmRes$table$logCPM, "F" = glmRes$table$F)
  list("pValMat" = pValMat, "dispEsts" = dgeW$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeRRobust

edgeR_standard <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dge <- estimateDisp(dge, design)
  glmFit <- glmFit(dge, design)
  glmRes <- glmLRT(glmFit, coef = 2)
  pval <- glmRes$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmRes$table$logFC, "logCPM" = glmRes$table$logCPM, "F" = glmRes$table$F)
  list("pValMat" = pValMat, "dispEsts" = dge$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeR - Standard

edgeR_zinbweights <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights)
{
  # require(edgeR)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  dge$weights <- weights
  dge = estimateDisp(dge,design)
  #plotBCV(dge)
  glmFit <- glmFit(dge, design)
  glmlrt <- glmWeightedF(glmFit,coef = 2)
  pval <- glmlrt$table$PValue
  padj <- p.adjust(pval, "BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = taxa_names(physeq)
  statInfo <- cbind("logFC" = glmlrt$table$logFC, "logCPM" = glmlrt$table$logCPM, "F" = glmlrt$table$F)
  list("pValMat" = pValMat, "dispEsts" = dge$tagwise.dispersion, "statInfo" = statInfo)
}# END: edgeR - ZINBWaVE weights

limma_voom <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  # require(edgeR)
  # require(limma)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }# END: limma-voom
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot=FALSE, lib.size = colSums(counts)*NFs)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none")
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END: limma voom

limma_voom_zinbweights <- function(physeq, design = as.formula("~ grp"), normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"), weights)
{
  # require(edgeR)
  # require(limma)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  counts <- as(otu_table(physeq), "matrix")
  if( normFacts=="TSS"){
    NFs = 1
  } else {
    normFacts <- paste("NF", normFacts, sep = ".")
    NFs = get_variable(physeq, normFacts)
    NFs = NFs/exp(mean(log(NFs)))
  }
  
  dge <- DGEList(counts = counts) #, remove.zeros=TRUE)
  dge$samples$norm.factors <- NFs
  design <- model.matrix(design, data.frame(sample_data(physeq)))
  v <- voom(dge, design, plot = FALSE, weights = weights, lib.size = colSums(counts)*NFs)
  v$weights <- v$weights * weights
  fit <- lmFit(v, design, weights = v$weights)
  fit$df.residual <- rowSums(weights) - ncol(design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="none")
  pval <- tt$P.Value
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(tt)
  statInfo <- cbind("logFC" = tt$logFC, "t" = tt$t)
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END: limma voom - ZINBWaVE weights

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund.
negBinTestDESeq2 <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                             normFacts = c("TMM", "RLE", "poscounts","ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE)
{
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #if (any(otu_table(physeq) == 0))
  #{
  #  otu_table(physeq) <- otu_table(physeq) + 1L
  #} else {}
  dds <- phyloseq_to_deseq2(physeq, design = design)
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### Run DESeq
  ddsRes <- DESeq(object = dds,test = "LRT", reduced = ~1,parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat, "dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund. accounting for zero inflation through zinb weights
negBinTestDESeq2_zinbweights <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                                         normFacts = c("TMM", "RLE", "poscounts", "ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights)
{
  register(SerialParam())
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  # if (any(otu_table(physeq) == 0))
  # {
  #  otu_table(physeq) <- otu_table(physeq) + 0.0001L
  # } else {}
  
  dds <- phyloseq_to_deseq2(physeq, design = design)
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### ZINB-WaVE weights
  counts <- as(otu_table(physeq), "matrix")
  weights[which(weights<1e-6)] <- 1e-06
  assays(dds)[["weights"]] = weights
  ### Run DESeq
  
  ddsRes <- DESeq(object = dds, test = "LRT", reduced = ~1, parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat,"dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2 + ZINBWaVE

### Perform ZIG regression from metagenomeSeq
metagenomeSeqZIG <- function (physeq, design = "~ grp", 
                              normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  # require(metagenomeSeq)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## OTU table
  otuTab <- as(otu_table(physeq), "matrix")
  ## sample data converted to Annotated Data Frame
  ADF <- AnnotatedDataFrame(data.frame(sample_data(physeq)))
  ## extract the chosen Normalisation Factors
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs <- get_variable(physeq, normFacts)
  if(normFacts == "NF.TMM"){
    NFs = NFs * sample_sums(physeq)
  }
  if (all(NFs==1))
  {
    ## needs one normalisation factor, library size in this case
    MGS <- newMRexperiment(counts = otuTab, phenoData = ADF, 
                           normFactors = colSums(otuTab))
  } else{
    MGS <- newMRexperiment(counts = otuTab, phenoData = ADF, normFactors = NFs)
  }# END - ifelse: normalisation factors all 1 or not
  
  ## design matrix
  design <- model.matrix(as.formula(paste0(design,paste0("+",normFacts),collapse = "")), data = data.frame(sample_data(physeq)))
  
  suppressWarnings(fit <- try(fitZig(MGS, design,verbose = FALSE,useCSSoffset = FALSE, control = zigControl(maxit = 1000)), silent = TRUE))
  
  if(class(fit)=="try-error"){
    res=matrix(NA, ncol=2, nrow=ntaxa(physeq))
  } else {
  # You need to specify all OTUs to get the full table from MRfulltable. 
  res <- MRfulltable(fit, number = nrow(get("counts", assayData(MGS))))
  # if any OTUs left out, rm those from x. Detected by NA rownames.
  res <- res[!is.na(rownames(res)), c("pvalues", "adjPvalues"), drop = FALSE]
  }
  colnames(res) <- c("rawP", "adjP")
  pValMat <- as.matrix(res)
  lods <- fit$eb$lods
  Ftest <- fit$eb$F
  names(Ftest) <- taxa_names(physeq)
  colnames(lods) <- c("Intercept","coef_grp2","scalingFactor")
  statInfo <- cbind("lods" = lods, "F" = Ftest)
  return(list("pValMat" = pValMat,"statInfo" = statInfo))
}# END - function: metagenomeSeqZIG

MASTmodel <- function(physeq,design = as.formula("~ grp"))
{
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  counts <- as(otu_table(physeq), "matrix")
  tpm <- counts*median(colSums(counts))/colSums(counts)
  tpm <- log2(tpm + 1)
  sca <- FromMatrix(tpm,cData=sample_data(physeq))
  # here, we keep all OTUs so that we can fairly compare MAST and the other methods. So, no adaptive thresholding or filtering by gene expression
  assays(sca) <- list(tpm=assay(sca))
  ngeneson <- apply(counts,2,function(x) mean(x>0))
  CD <- colData(sca)
  CD$ngeneson <- ngeneson
  CD$cngeneson <- CD$ngeneson-mean(ngeneson)
  colData(sca) <- CD
  ## differential expression
  fit <- zlm(~ grp + cngeneson, sca = sca,
             method = "bayesglm", ebayes = TRUE)
  summaryDt = summary(fit, doLRT='grpgrp2')
  summaryDt = summaryDt$datatable
  fcHurdle <- merge(summaryDt[contrast=='grpgrp2'&component=='H',.(primerid, `Pr(>Chisq)`)],
                    summaryDt[contrast=='grpgrp2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
  fcHurdle[,padj:=p.adjust(`Pr(>Chisq)`, 'BH')]
  df = data.frame(gene = fcHurdle$primerid,
                  pval=fcHurdle[,'Pr(>Chisq)'], padj=fcHurdle$padj,
                  logfc=fcHurdle$coef)
  colnames(df)[2] = 'pval'
  pValMat <- as.matrix(df[, c("pval", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = df$logfc)
  rownames(statInfo) <- rownames(pValMat) <- df$gene
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: MAST

scdemodel <- function(physeq,design = as.formula("~ grp"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  groupVar <- get_variable(physeq, "grp")
  counts <- as(otu_table(physeq), "matrix")
  names(groupVar) <- colnames(counts)
  counts <- apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x})
  o.ifm <- scde.error.models(counts = counts, groups = groupVar, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1,min.size.entries = nrow(counts)*0.2)
  # filter out cells (samples) that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.samples <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.samples, ]
  # estimate gene expression (OTU abundance) prior
  o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
  # Differential abundace analysis
  # define two groups of cells
  groups <- groupVar[valid.samples]
  # run differential expression tests on all genes.
  ediff <- scde.expression.difference(o.ifm, counts, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
  
  rawP <- 1-pnorm(ediff$Z)
  adjP <- p.adjust(rawP,method = "BH")
  pValMat <- as.matrix(cbind(rawP=rawP,adjP = adjP))
  rownames(pValMat) <- rownames(ediff)
  statInfo <- ediff[,1:4]
  list("pValMat" = pValMat, "statInfo" = statInfo)
}# END - function: Single-Cell Differential Expression Analysis scde

Seuratmodel <- function(physeq,design = as.formula("~ grp"),normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  # ## add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  
  groupVar <- get_variable(physeq, "grp")
  counts <- as(otu_table(physeq), "matrix")
  names(groupVar) <- colnames(counts)
  
  sobj <- CreateSeuratObject(raw.data = counts,min.cells = 1,
                             normalization.method = "LogNormalize",
                             do.scale = TRUE,do.center = TRUE)
  sobj <- AddMetaData(sobj,metadata = groupVar,col.name = "grp")
  sobj@ident <- as.factor(groupVar)
  # Gene selection for input to CCA
  sobj <- FindVariableGenes(sobj, do.plot = F)
  response <- FindMarkers(sobj, ident.1 = "grp1", ident.2 = "grp2",print.bar = FALSE)
  rawP <- response$p_val
  adjP <- p.adjust(rawP,method = "BH")
  pValMat <- as.matrix(cbind(rawP=rawP,adjP = adjP))
  rownames(pValMat) <- rownames(response)
  otu_na <- which(is.na(match(rownames(counts),rownames(pValMat))))
  if(!is.null(otu_na))
  {
    pValMat <- rbind(pValMat,matrix(NA,ncol = 2,nrow = length(otu_na),dimnames = list(c(rownames(counts)[otu_na]),c("rawP","adjP"))))
  }
  statInfo <- response[,2:4]
  return(list("pValMat" = pValMat, "statInfo" = statInfo))
}# END - function: Seurat Wilcoxon test

NODESmodel <- function(physeq,design = as.formula("~ grp"),normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  # ## add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  groupVar <- get_variable(physeq, "grp")
  counts <- as(otu_table(physeq), "matrix")
  colnames(counts)=groupVar
  normCounts=pQ(counts,frac = 1,throw_sd = 0,hard_outlier = 0)
  res=NODES::NODES(data=normCounts,group=colnames(normCounts))
  pval=vector(length=nrow(counts))
  names(pval)=rownames(counts)
  pval[rownames(normCounts)]=res$Fisher
  padj <- p.adjust(pval,method="BH")
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(normCounts)
  return(list("pValMat" = pValMat))
}

# Too time consuming
ANCOMmodel <- function(physeq)
{
  ### force orientation samples x OTUs
  if (taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  # ## add 1 to zero counts
  #   if (any(otu_table(physeq) == 0))
  #   {
  #     otu_table(physeq) <- otu_table(physeq) + 1L
  #   } else {}
  counts <- as(otu_table(physeq), "matrix")[,1:100]
  OTUdat <- data.frame(Sample.ID = sample_names(physeq),counts)
  groupVar <- get_variable(physeq, "grp")
  Vardat <- data.frame(Sample.ID = sample_names(physeq),grp = groupVar)
  
  comparison_test=ANCOM.main(OTUdat=OTUdat,
                             Vardat=Vardat,
                             adjusted=FALSE,
                             repeated=F,
                             main.var="grp",
                             adj.formula=NULL,
                             repeat.var=NULL,
                             longitudinal=FALSE,
                             random.formula=NULL,
                             multcorr=3,
                             sig=0.05,
                             prev.cut=0.90)
}# END - function: ANCOM2

ALDEx2model <- function(physeq,design = as.formula("~ grp"),normFacts = c("TMM", "RLE", "poscounts", "CSS", "UQ", "none", "TSS"))
{
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## add 1 to zero counts
  if (any(otu_table(physeq) == 0))
  {
    otu_table(physeq) <- otu_table(physeq) + 1L
  } else {}
  groupVar <- as.character(get_variable(physeq, "grp"))
  counts <- as(otu_table(physeq), "matrix")
  x <- aldex(counts, groupVar, mc.samples=128, test="t", effect=TRUE,
             include.sample.summary = FALSE, denom ="iqlr", verbose=TRUE)
  pval = x$wi.ep
  padj = x$wi.eBH
  pValMat <- cbind("rawP" = pval, "adjP" = padj)
  rownames(pValMat) = rownames(x)
  statInfo <- x[,1:7]
  return(list("pValMat" = pValMat,"statInfo" = statInfo))
}# END - function: ALDEx2

corncobmodel <- function(physeq, design = as.formula("~ grp"), test = c("Wald","LRT"), bootstrap = c("TRUE","FALSE")){
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  da_analysis <- differentialTest(formula = design,
                                  phi.formula = design,
                                  formula_null = ~ 1,
                                  phi.formula_null = design,
                                  test = test, boot = bootstrap,
                                  data = physeq,
                                  fdr_cutoff = 0.05)
  pValMat <- cbind("rawP" = da_analysis$p, "adjP" = da_analysis$p_fdr)
  rownames(pValMat) = names(da_analysis$p)
  statInfo <- ldply(da_analysis$all_models,coef)
  statInfo$coef = c("mu.(Intercept)","mu.grpgrp2","phi.(Intercept)","phi.grpgrp2")
  statInfo$feature = rep(names(da_analysis$p),each = 4)
  return(list("pValMat" = pValMat,"statInfo" = statInfo))
}# END - function: corncob

# Too time consuming
selbalmodel <- function(physeq, variable_name = "grp"){
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  selbal_test <- selbal.cv(x = t(physeq@otu_table),
                           y = physeq@sam_data[,variable_name])
}# END - function: selbal

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

oneSimRunGSOwn <- function(physeq, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  ## all normalisations
  physeq <- normEdgeR(physeq = physeq, method = "TMM")
  # physeq <- normEdgeR(physeq = physeq, method = 'RLE') 
  # physeq <- normEdgeR(physeq = physeq, method = 'upperquartile')
  physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
  physeq <- normCSS(physeq = physeq)
  # physeq <- normTSS(physeq = physeq)
  cat("Normalisations: DONE\n")
  #zinb model estimation
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$grp), K = 0,
                       epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  cat("ZINB model estimation: DONE\n")
  returnList = list()
  #returnList$physeq = physeq
  returnList = within(returnList, {
    ## edgeR Robust
    edgeR_TMM_robustDisp <- edgeR_robust(physeq, normFacts = "TMM")
    cat("EdgeR robust tests: DONE\n")
    ## edgeR Standard
    edgeR_TMM_standard <- edgeR_standard(physeq, normFacts = "TMM")
    cat("EdgeR standard tests: DONE\n")
    ## edgeR with zinbwave weights
    edgeR_TMM_zinbwave <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = weights)
    cat("EdgeR with ZINB-WaVE weights tests: DONE\n")
    ## edgeR with true weights
    if(!is.null(true_weights)){
      edgeR_TMM_trueweights <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = true_weights)
      cat("EdgeR with true weights tests: DONE\n")
    }
    edgeR_poscounts_standard <- edgeR_standard(physeq, normFacts = "poscounts")
    cat("EdgeR with RLE normalisation tests: DONE\n")
    
    ## limma-voom
    limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
    cat("Limma Voom tests: DONE\n")
    ## limma-voom zinbwave weights
    limma_voom_TMM_zinbwave <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = weights)
    cat("Limma Voom ZINB-WaVE weights tests: DONE\n")
    ## limma-voom with true weights
    if(!is.null(true_weights)){
      limma_voom_TMM_trueweights <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = true_weights)
      cat("Limma Voom true weights tests: DONE\n")
    }
    
    ## NB test from DESeq2
    DESeq2_poscounts <- negBinTestDESeq2(physeq, normFacts = "poscounts")
    cat("NB DESeq2 tests: DONE\n")
    ## NB test from DESeq2 with zinbwave weights
    DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = weights)
    cat("NB DESeq2 with ZINB-WaVE weights tests: DONE\n")
    ## NB test from DESeq2 with true weights
    if(!is.null(true_weights)){
      DESeq2_poscounts_trueweights <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = true_weights)
      cat("NB DESeq2 with true weights tests: DONE\n")
    }
    DESeq2_TMM <- negBinTestDESeq2(physeq, normFacts = "TMM")
    cat("NB DESeq2 with TMM normalisation tests: DONE\n")
    
    ## metagenomeSeq Zero-Inflated Gaussian
    mgsZig_CSS <- metagenomeSeqZIG(physeq, normFacts = "CSS")
    cat("MetagenomeSeq ZIG tests: DONE\n")
    ## ALDEx2 Wilcoxon test
    ALDEx2 <- ALDEx2model(physeq)
    cat("ALDEx2 Wilcoxon test: DONE\n")
    ## corncob Wald test
    corcob_wald <- corncobmodel(physeq,test = "Wald",bootstrap = FALSE)
    cat("corncob Wald test: DONE\n")
    ## corncob LRT test
    corcob_LRT <- corncobmodel(physeq,test = "LRT",bootstrap = FALSE)
    cat("corncob LRT test: DONE\n")
      
    ## MAST hurdle models
    MAST <- MASTmodel(physeq)
    cat("MAST LRT tests: DONE\n")
    ## scde single cell differential expression
    scde <- scdemodel(physeq)
    cat("scde single cell differential expression tests: DONE\n")
    ## Seurat 
    seurat_wilcoxon <- Seuratmodel(physeq)
    cat("Seurat Wilcoxon tests: DONE\n")
    ## NODES
    #nodes <- NODESmodel(physeq)
    #cat("NODES Wilcoxon tests: DONE\n")
  })
  return(returnList)
}

oneSimRunGSOwnFastestMethod <- function(physeq, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  physeq <- normEdgeR(physeq = physeq, method = "TMM")
  cat("Normalisations: DONE\n")
  returnList = list()
  #returnList$physeq = physeq
  returnList = within(returnList, {
    ## edgeR Standard
    edgeR_TMM_standard <- edgeR_standard(physeq, normFacts = "TMM")
    cat("EdgeR standard tests: DONE\n")
    ## limma-voom
    limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
    cat("Limma Voom tests: DONE\n")
  })
  return(returnList)
}

oneSimRunGSOwn_time <- function(physeq, true_weights = NULL, epsilon = 1e10) { 
  # Prevent NA when converting to integer due to some outlier generation during simulation
  physeq@otu_table@.Data[which(physeq@otu_table@.Data>.Machine$integer.max)] <- .Machine$integer.max
  ## all normalisations
  
  #zinb model estimation
  zinb_start <- proc.time()
  zinbmodel <- zinbFit(Y = physeq@otu_table@.Data, 
                       X = model.matrix(~ physeq@sam_data$grp), K = 0,
                       epsilon = epsilon, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  weights <- computeExactWeights(model = zinbmodel,x = physeq@otu_table@.Data)
  colnames(weights) <- colnames(physeq@otu_table)
  rownames(weights) <- rownames(physeq@otu_table)
  zinb_end <- proc.time()
  zinb_time <- zinb_end-zinb_start
  
  returnList = list()
  returnList = within(returnList, {
  
    edgeR_TMM_robustDisp <- system.time(
      {## edgeR Robust
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_robustDisp <- edgeR_robust(physeq, normFacts = "TMM")
      cat("EdgeR robust tests: DONE\n")}
    )
    
    edgeR_TMM_standard <- system.time(
      {## edgeR Standard
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_standard <- edgeR_standard(physeq, normFacts = "TMM")
      cat("EdgeR standard tests: DONE\n")}
    )
    
    edgeR_TMM_zinbwave <- system.time(
      {## edgeR with zinbwave weights
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      edgeR_TMM_zinbwave <- edgeR_zinbweights(physeq, normFacts = "TMM", weights = weights)
      cat("EdgeR with ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    edgeR_poscounts_standard <- system.time(
      {physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
      edgeR_poscounts_standard <- edgeR_standard(physeq, normFacts = "poscounts")
      cat("EdgeR with RLE normalisation tests: DONE\n")}
    )
    
    limma_voom_TMM <- system.time(
      {## limma-voom
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      limma_voom_TMM <- limma_voom(physeq, normFacts = "TMM")
      cat("Limma Voom tests: DONE\n")}
    )
    
    limma_voom_TMM_zinbwave <- system.time(
      {## limma-voom zinbwave weights
      physeq <- normEdgeR(physeq = physeq, method = "TMM")
      limma_voom_TMM_zinbwave <- limma_voom_zinbweights(physeq, normFacts = "TMM", weights = weights)
      cat("Limma Voom ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    DESeq2_poscounts <- system.time(
      {## NB test from DESeq2
        physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
        DESeq2_poscounts <- negBinTestDESeq2(physeq, normFacts = "poscounts")
        cat("NB DESeq2 tests: DONE\n")}
    )
    
    DESeq2_TMM <- system.time(
      {## DESeq2 TMM
        physeq <- normEdgeR(physeq = physeq, method = "TMM")
        DESeq2_TMM <- negBinTestDESeq2(physeq, normFacts = "TMM")
        cat("NB DESeq2 with TMM normalisation tests: DONE\n")}
    )
    
    DESeq2_poscounts_zinbwave <- system.time(
      {physeq <- normDESeq2(physeq = physeq)  # poscounts, similar to RLE
        ## NB test from DESeq2 with zinbwave weights
        DESeq2_poscounts_zinbwave <- negBinTestDESeq2_zinbweights(physeq, normFacts = "poscounts",weights = weights)
        cat("NB DESeq2 with ZINB-WaVE weights tests: DONE\n")}
    ) + zinb_time
    
    mgsZig_CSS <- system.time(
      {## metagenomeSeq Zero-Inflated Gaussian
        physeq <- normCSS(physeq = physeq)
        mgsZig_CSS <- metagenomeSeqZIG(physeq, normFacts = "CSS")
        cat("MetagenomeSeq ZIG tests: DONE\n")}
    )
    
    ALDEx2 <- system.time(
      {## ALDEx2 Wilcoxon test
       ALDEx2 <- ALDEx2model(physeq)
       cat("ALDEx2 Wilcoxon test: DONE\n")}
    )
    
    MAST <- system.time(
      {## MAST hurdle models
       MAST <- MASTmodel(physeq)
       cat("MAST lrt tests: DONE\n")}
    )
    
    scde <- system.time(
      {## scde single cell differential expression
        scde <- scdemodel(physeq)
        cat("scde single cell differential expression tests: DONE\n")}
    )
    
    seurat_wilcoxon <- system.time(
      {## Seurat 
        seurat_wilcoxon <- Seuratmodel(physeq)
        cat("Seurat Wilcoxon tests: DONE\n")}
    )
  })
  return(returnList)
}