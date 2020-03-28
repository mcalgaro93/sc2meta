### To aggregate the single output of the job array we use a list

# 16S
Stool_16S_mockDA <- list()
for(i in 1:1000){
  Stool_16S_mockDA[[i]] <- readRDS(paste0("./mockDA/Stool_16S_mockDA_",i,".RDS"))
}
names(Stool_16S_mockDA) <- paste0("Comparison",1:1000)
saveRDS(Stool_16S_mockDA,file = "../data/Stool_16S_mockDA.RDS")

# WMS
Stool_WMS_mockDA <- list()
for(i in 1:1000){
  Stool_WMS_mockDA[[i]] <- readRDS(paste0("./mockDA/Stool_WMS_mockDA_",i,".RDS"))
}
names(Stool_WMS_mockDA) <- paste0("Comparison",1:1000)
saveRDS(Stool_WMS_mockDA,file = "../data/Stool_WMS_mockDA.RDS")

