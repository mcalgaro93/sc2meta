library(RColorBrewer)

# RColorBrewer::display.brewer.all()

cols <- c(
# ALDEx2
"#C0E845",
# 4 DESeq2
brewer.pal(n = 9,"YlOrBr")[c(4,7,7,8,9,5)],
# 5 edgeR
brewer.pal(n = 11, "BrBG")[c(8,9,11)],
brewer.pal(n = 11, "PiYG")[c(9,11,10)],
# 3 limma
brewer.pal(n = 9, "Purples")[c(5,7,9,8)],
# MAST and metagenomeSeq
brewer.pal(n = 11, "PiYG")[c(2,4)],
# scde seurat
brewer.pal(n = 11, "RdGy")[c(8,10)])

methods <- c("ALDEx2",
             "DESeq2_poscounts",
             "DESeq2_poscounts_trueweights",
             "deSeq2_poscounts_trueweights",
             "DESeq2_poscounts_zinbwave",
             "DESeq2_poscounts_zinbwave_epsilon1e14",
             "DESeq2_TMM",
             "edgeR_poscounts_standard",
             "edgeR_TMM_standard",
             "edgeR_TMM_robustDisp",
             "edgeR_TMM_trueweights",
             "edgeR_TMM_zinbwave",
             "edgeR_TMM_zinbwave_epsilon1e14",
             "limma_voom_TMM",
             "limma_voom_TMM_trueweights",
             "limma_voom_TMM_zinbwave",
             "limma_voom_TMM_zinbwave_epsilon1e14",
             "MAST",
             "mgsZig_CSS",
             "scde",
             "seurat_wilcoxon")

names(cols) <- methods
