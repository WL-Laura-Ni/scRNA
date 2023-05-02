# author information----
### Annotation with cell marker genes
### writer: Wanlin Laura NI
### time: 2023-04-26

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
integ_obj <- readRDS("~/Graduation_proj/02-reduction/integ_7_set_umap_obj.rds")

# cell marker feature plot ----
pdf("./feature_plot_markers.pdf",width = 8,height = 6)
DefaultAssay(integ_obj) <- "integrated"
DimPlot(integ_obj,group.by = "integrated_snn_res.1.8",label = T)+NoLegend()

DefaultAssay(integ_obj) <- "RNA"
integ_obj <- NormalizeData(integ_obj)

#* plasma cells ----
FeaturePlot(integ_obj,features = "MZB1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "SDC1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "IGHG1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* B cells ----
FeaturePlot(integ_obj,features = "CD79A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD79B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "MS4A1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* NK cells ----
FeaturePlot(integ_obj,features = "NKG7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "GNLY",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* CD8 T cells ----
FeaturePlot(integ_obj,features = "CD8A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD8B",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* CD4 T cells ----
FeaturePlot(integ_obj,features = "CD4",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "IL7R",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CCR7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD3E",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CD3D",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Macrophages
FeaturePlot(integ_obj,features = "FCGR3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

# Monocytes
FeaturePlot(integ_obj,features = "CD14",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "LYZ",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "MS4A7",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "FCGR3A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

#* Dendritic cells ----
FeaturePlot(integ_obj,features = "FCER1A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "CLEC10A",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "IL3RA",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "GZMB",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "SERPINF1",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

FeaturePlot(integ_obj,features = "ITM2C",
            pt.size = 0.1,
            min.cutoff = "q10",
            max.cutoff = "q95")+NoLegend()

dev.off()

# tumor immune annotation ----
Idents(integ_obj) <- "integrated_snn_res.1.8"
integ_obj$celltype1 <- "Tumor cells"
integ_obj$celltype1[which(integ_obj$integrated_snn_res.1.8 == 4 |
                            integ_obj$integrated_snn_res.1.8 == 8 |
                            integ_obj$integrated_snn_res.1.8 == 9 |
                            integ_obj$integrated_snn_res.1.8 == 15 |
                            integ_obj$integrated_snn_res.1.8 == 20 |
                            integ_obj$integrated_snn_res.1.8 == 22 |
                            integ_obj$integrated_snn_res.1.8 == 25 |
                            integ_obj$integrated_snn_res.1.8 == 27 |
                            integ_obj$integrated_snn_res.1.8 == 38)] <- "Immune cells"

# subset cells ----
Idents(integ_obj) <- "celltype1"
VlnPlot(integ_obj,features = c("nUMI","nGene","log10GenesPerUMI"),pt.size = 0)+NoLegend()
tumor <- subset(integ_obj,celltype1 == "Tumor cells")
immu <- subset(integ_obj,celltype1 == "Immune cells")

saveRDS(integ_obj,file = "./integ_obj_end.rds")
saveRDS(tumor,file = "./tumor_raw.rds")
saveRDS(immu,file = "./immu_raw.rds")

rm(list=ls())
gc()
