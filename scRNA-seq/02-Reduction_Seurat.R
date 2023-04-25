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
setwd("~/Graduation_proj/02-reduction_annotation/")

# data loading ----
obj <- readRDS("~/Graduation_proj/01-QC/filtered_merged_obj.rds")

# cell phase contribution to the variance ----
Phase_obj <- NormalizeData(obj)
Phase_obj <- CellCycleScoring(Phase_obj,s.features = cc.genes$s.genes,
                              g2m.features = cc.genes$g2m.genes)
Phase_obj <- FindVariableFeatures(Phase_obj)
Phase_obj <- ScaleData(Phase_obj)
Phase_obj <- RunPCA(Phase_obj)
DimPlot(Phase_obj,
        reduction="pca",
        group.by = "Phase",
        split.by = "Phase")+
  ggtitle("Phase contribution to variance")

# seurat advised workflow for UMAP reduction ----
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures = 2000)
obj <- ScaleData(obj,features = rownames(obj))
obj <- RunPCA(obj,features=VariableFeatures(object=obj))
ElbowPlot(obj,ndim=35)+
  geom_hline(yintercept = 1.5,linetype="dashed")
obj <- FindNeighbors(obj,dims=1:24)
obj <- FindClusters(obj,resolution = 1.8)
obj <- RunUMAP(obj,dims = 1:24)

saveRDS(obj,file = "./obj_seurat_umap.rds")

# filter BCR clonotype ----
obj$clonotype_id <- paste0(obj$orig.ident,"-",obj$clonotype_id)
obj$clonotype_id[grep("NA$",obj$clonotype_id)] <- NA
tmp <- table(obj$clonotype_id)
tmp <- subset(tmp,tmp > 5)
obj$filter_clonotype <- NA
obj$filter_clonotype[which(obj$clonotype_id %in% names(tmp))] <- obj$clonotype_id[which(obj$clonotype_id %in% names(tmp))]

# load cell phase score from the integrated obj ----
integ_obj <- readRDS("./integ_obj_pca.rds")
obj$Phase <- integ_obj$Phase
obj$S.Score <- integ_obj$S.Score
obj$G2M.Score <- integ_obj$G2M.Score

# plot UMAP ----
pdf(file="./reduction_seurat_umap_1.8.pdf",width=8,height=6)

DimPlot(obj,reduction = "umap",group.by = "orig.ident",label = T)+
  ggtitle("UMAP colored by patients in general seurat workflow")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size=20))

DimPlot(obj,reduction = "umap",label = T)+
  ggtitle("UMAP colored by clusters in general seurat workflow")+
  theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=15))

DimPlot(obj,reduction = "umap",group.by = "filter_clonotype",label = T,repel = T)+
  ggtitle("UMAP colored by BCR clonotype in general seurat workflow")+
  theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=15))

dev.off()

# adjust resolution ----
obj <- FindClusters(obj,resolution = 0.6)
Idents(obj) <- "RNA_snn_res.0.6" 

obj.markers <- FindAllMarkers(obj,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
obj.markers %>%
  group_by(cluster) %>%
  slice_max(n=4,order_by = avg_log2FC)

# plot UMAP ----
pdf(file="./reduction_seurat_umap_0.6.pdf",width=8,height=6)

DimPlot(obj,reduction = "umap",label = T)+
  ggtitle("UMAP colored by clusters in general seurat workflow")+
  theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))

DimPlot(obj,reduction = "umap",label = T,split.by = "DrugRes")+
  ggtitle("UMAP colored by clusters in general seurat workflow")+
  theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))

DimPlot(obj,reduction = "umap",label = T,split.by = "Phase")+
  ggtitle("UMAP colored by clusters in general seurat workflow")+
  theme(legend.position="None",plot.title = element_text(hjust = 0.5,face = "bold",size=20))

FeaturePlot(obj,
            reduction="umap",
            features="nUMI",
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

FeaturePlot(obj,
            reduction="umap",
            features="nGene",
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

FeaturePlot(obj,
            reduction="umap",
            features="percentMT",
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

FeaturePlot(obj,
            reduction="umap",
            features="S.Score",
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

FeaturePlot(obj,
            reduction="umap",
            features="G2M.Score",
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

FeaturePlot(obj,
            reduction="umap",
            features=c("nUMI","nGene","percentMT","S.Score","G2M.Score"),
            pt.size = 0.1,
            min.cutoff = "q5",
            label = T)+
  theme(legend.position="None")

dev.off()

# PC contribution in UMAP ----
columns <- c(paste0("PC_",1:24),
             "ident",
             "UMAP_1",
             "UMAP_2")

pc_data <- FetchData(obj,vars = columns)
umap_label <- FetchData(obj,
                        vars = c("ident","UMAP_1","UMAP_2")) %>%
  group_by(ident) %>%
  summarize(x=mean(UMAP_1),y=mean(UMAP_2))

pdf("./PC-in-UMAP.pdf",width=8,height=6)
map(paste0("PC_",1:24),function(pc){
  ggplot(pc_data,aes(UMAP_1,UMAP_2))+
    geom_point(aes_string(color=pc),alpha=0.5,size=0.1)+
    scale_color_gradient(guide = F,low = "grey90",high = "blue")+
    geom_text(data = umap_label,aes(label=ident,x,y))+
    theme_classic()+
    ggtitle(pc)
}) 
dev.off()

# cell marker feature plot ----
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
pdf("./feature_plot_markers.pdf",width = 8,height = 6)

# plasma cells
FeaturePlot(obj,features = "MZB1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "SDC1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "IGHG1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# B cells
FeaturePlot(obj,features = "CD79A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD79B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "MS4A1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# NK cells
FeaturePlot(obj,features = "NKG7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "GNLY",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# CD8 T cells
FeaturePlot(obj,features = "CD8A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD8B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# CD4 T cells
FeaturePlot(obj,features = "CD4",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "IL7R",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CCR7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Macrophages
FeaturePlot(obj,features = "FCGR3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Monocytes
FeaturePlot(obj,features = "CD14",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "LYZ",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "MS4A7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "FCER3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Dendritic cells
FeaturePlot(obj,features = "FCER1A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "CLEC10A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "IL3RA",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "GZMB",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "SERPINF1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(obj,features = "ITM2C",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

dev.off()
gc()