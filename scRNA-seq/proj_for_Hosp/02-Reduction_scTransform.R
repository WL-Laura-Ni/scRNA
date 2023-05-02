# author information----
### Reduction
### writer: Wanlin Laura NI
### time: 2023-04-11

# working condition setup----
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)
library(scales)
library(RCurl)
setwd("~/Graduation_proj/02-reduction/")

# data loading ----
obj <- readRDS("~/Graduation_proj/01-QC/raw_merged_obj.rds")

# SCTransform for reduction ----
tumor_list <- SplitObject(obj,split.by = "orig.ident")

tumor_list <- lapply(X = tumor_list,FUN = function(x){
  x <- NormalizeData(x)
  x <- CellCycleScoring(x,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes)
  x <- SCTransform(x,
                   vars.to.regress = c("percentMT"))
})

integ_features <- SelectIntegrationFeatures(object.list = tumor_list,
                                            nfeatures = 3000)
split_obj <- PrepSCTIntegration(object.list = tumor_list,
                                anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_obj,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
integ_obj <- IntegrateData(anchorset = integ_anchors,
                           normalization.method = "SCT")
saveRDS(integ_obj,file = paste0("./s-in7.rds"))

#* UMAP reduction ----
integ_obj <- RunPCA(object = integ_obj)

# pca number calculation
pct <- integ_obj[["pca"]]@stdev/sum(integ_obj[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1),decreasing = T)[1]+1
message("------------------------------------------------")
message(paste0("-----------------",co2,"-----------------"))
message("------------------------------------------------")
# calulated PC is 14

# filter BCR clonotype ----
integ_obj$clonotype_id <- paste0(integ_obj$orig.ident,"-",integ_obj$clonotype_id)
integ_obj$clonotype_id[grep("NA$",integ_obj$clonotype_id)] <- NA
tmp <- table(integ_obj$clonotype_id)
tmp <- subset(tmp,tmp > 5)
integ_obj$filter_clonotype <- NA
integ_obj$filter_clonotype[which(integ_obj$clonotype_id %in% names(tmp))] <- integ_obj$clonotype_id[which(integ_obj$clonotype_id %in% names(tmp))]

# plot UMAP ----
integ_obj <- FindNeighbors(object=integ_obj,dims = 1:40)
integ_obj <- FindClusters(object = integ_obj,resolution=c(1.4,1.8,2.2))
integ_obj <- RunUMAP(integ_obj,
                     dims = 1:40,
                     reduction="pca")

Idents(object = integ_obj) <- "integrated_snn_res.1.4"
saveRDS(integ_obj,file = "./integ_7_set_umap_obj.rds")

if (FALSE){
  pdf(file = "./reduction_sct_umap_1.4.pdf",width = 8,height = 6)
  
  DimPlot(integ_obj,label = T)+
    ggtitle("UMAP colored by clusters in sctransform")+
    theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  
  DimPlot(integ_obj,label = T,split.by = "DrugRes")+
    ggtitle("UMAP colored by drug resistance in sctransform")+
    theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  
  DimPlot(integ_obj,label = T,group.by = "filter_clonotype",repel=T)+
    ggtitle("UMAP colored by clonotype in sctransform")+
    theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  
  DimPlot(integ_obj,label = T,group.by = "orig.ident",repel = T)+
    ggtitle("UMAP colored by patient in sctransform")+
    theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  
  FeaturePlot(integ_obj,
              reduction="umap",
              features="nUMI",
              pt.size = 0.1,
              min.cutoff = "q5",
              label = T)+
    theme(legend.position="None")
  
  FeaturePlot(integ_obj,
              reduction="umap",
              features="nGene",
              pt.size = 0.1,
              min.cutoff = "q5",
              label = T)+
    theme(legend.position="None")
  
  FeaturePlot(integ_obj,
              reduction="umap",
              features="percentMT",
              pt.size = 0.1,
              min.cutoff = "q5",
              label = T)+
    theme(legend.position="None")
  
  FeaturePlot(integ_obj,
              reduction="umap",
              features="S.Score",
              pt.size = 0.1,
              min.cutoff = "q5",
              label = T)+
    theme(legend.position="None")
  
  FeaturePlot(integ_obj,
              reduction="umap",
              features="G2M.Score",
              pt.size = 0.1,
              min.cutoff = "q5",
              label = T)+
    theme(legend.position="None")
  
  
  dev.off()
  
  # PC contribution in UMAP ----
  columns <- c(paste0("PC_",1:40),
               "ident",
               "UMAP_1",
               "UMAP_2")
  
  pc_data <- FetchData(integ_obj,vars = columns)
  umap_label <- FetchData(integ_obj,
                          vars = c("ident","UMAP_1","UMAP_2")) %>%
    group_by(ident) %>%
    summarize(x=mean(UMAP_1),y=mean(UMAP_2))
  
  pdf("./PC-in-SCT.pdf",width=8,height=6)
  map(paste0("PC_",1:40),function(pc){
    ggplot(pc_data,aes(UMAP_1,UMAP_2))+
      geom_point(aes_string(color=pc),alpha=0.5,size=0.1)+
      scale_color_gradient(guide = F,low = "grey90",high = "blue")+
      geom_text(data = umap_label,aes(label=ident,x,y))+
      theme_classic()+
      ggtitle(pc)
  }) 
  dev.off()
  
}

rm(list=ls())
gc()
