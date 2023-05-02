# author information----
### Annotation with cell marker genes
### writer: Wanlin Laura NI
### time: 2023-05-02

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
setwd("~/Graduation_proj/04-immu/")

# data loading ----
immu <- readRDS("~/Graduation_proj/04-immu/scTrans_immu_umap.rds")

# cell marker feature plot ----
pdf("./feature_plot_markers.pdf",width = 8,height = 6)
DefaultAssay(immu) <- "integrated"
DimPlot(immu,group.by = "integrated_snn_res.1.2",label = T)+NoLegend()

DefaultAssay(immu) <- "RNA"
immu <- NormalizeData(immu)

#* plasma cells ----
FeaturePlot(immu,features = "MZB1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "SDC1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "IGHG1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* B cells ----
FeaturePlot(immu,features = "CD79A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD79B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "MS4A1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* NK cells ----
FeaturePlot(immu,features = "NKG7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "GNLY",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* CD8 T cells ----
FeaturePlot(immu,features = "CD8A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD8B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* CD4 T cells ----
FeaturePlot(immu,features = "CD4",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "IL7R",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CCR7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Macrophages
FeaturePlot(immu,features = "FCGR3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Monocytes
FeaturePlot(immu,features = "CD14",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "LYZ",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "MS4A7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "FCGR3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* Dendritic cells ----
FeaturePlot(immu,features = "FCER1A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "CLEC10A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "IL3RA",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "GZMB",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "SERPINF1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(immu,features = "ITM2C",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

dev.off()

# tumor immune annotation ----
Idents(immu) <- "integrated_snn_res.1.2"
immu$celltype <- NA
immu$celltype[which(immu$integrated_snn_res.1.2 == 1 |
                      immu$integrated_snn_res.1.2 == 3 |
                      immu$integrated_snn_res.1.2 == 4 |
                      immu$integrated_snn_res.1.2 == 13 |
                      immu$integrated_snn_res.1.2 == 15 |
                      immu$integrated_snn_res.1.2 == 21 )] <- "CD4 T cell"
immu$celltype[which(immu$integrated_snn_res.1.2 == 0 |
                      immu$integrated_snn_res.1.2 == 2 |
                      immu$integrated_snn_res.1.2 == 10 |
                      immu$integrated_snn_res.1.2 == 12 |
                      immu$integrated_snn_res.1.2 == 17 |
                      immu$integrated_snn_res.1.2 == 18 |
                      immu$integrated_snn_res.1.2 == 20)] <- "CD8 T cell"
immu$celltype[which(immu$integrated_snn_res.1.2 == 5 |
                      immu$integrated_snn_res.1.2 == 7 |
                      immu$integrated_snn_res.1.2 == 11 |
                      immu$integrated_snn_res.1.2 == 14 |
                      immu$integrated_snn_res.1.2 == 19 )] <- "NK cell"
immu$celltype[which(immu$integrated_snn_res.1.2 == 6 |
                      immu$integrated_snn_res.1.2 == 16 |
                      immu$integrated_snn_res.1.2 == 22)] <- "Monocyte"
immu$celltype[which(immu$integrated_snn_res.1.2 == 8)] <- "B cell"
immu$celltype[which(immu$integrated_snn_res.1.2 == 9)] <- "plasma cell"

# subset cells ----
Idents(immu) <- "celltype"
DimPlot(immu,group.by = "celltype",label = T,repel = T) + NoLegend()
DimPlot(immu,group.by = "celltype",label = T,repel = T,split.by = "orig.ident") + NoLegend()
VlnPlot(immu,features = c("nUMI","nGene","log10GenesPerUMI"),pt.size = 0)+NoLegend()

saveRDS(immu,file = "./immu_anno.rds")

rm(list=ls())
gc()
