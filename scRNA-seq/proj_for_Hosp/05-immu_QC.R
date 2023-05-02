# author information----
### QC for subset immu data
### writer: Wanlin Laura NI
### time: 2023-04-27

# working condition setup----
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
setwd("~/Graduation_proj/04-immu/")

# data loading ----
immu <- readRDS("~/Graduation_proj/02-reduction/immu_raw.rds")
Idents(immu) <- "orig.ident"

# plot raw data information for QC----
metadata <- as.data.frame(immu@meta.data)

pdf(file = "./immu-QC-plot.pdf",width = 8,height = 6)

# detected cells/patient
metadata %>% 
  ggplot(aes(x=orig.ident,fill=DrugRes))+
  geom_bar() +
  geom_hline(yintercept = 1500,linetype="dashed") +
  geom_hline(yintercept = 3500,linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20)) +
  xlab("Patients") +
  ylab("Cell Count")+
  ggtitle("Numbers of detected cells")

# nGene of the patients
VlnPlot(immu,
        features = "nGene",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 3000,linetype="dashed")+
  geom_hline(yintercept = 300,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# nUMI of the patients
VlnPlot(immu,
        features="nUMI",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 700,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# percent of MT genes
VlnPlot(immu,
        features = "percentMT",
        group.by = "orig.ident",
        pt.size = 0)+
  theme(legend.position = "None",plot.title = element_text(size = 20)) +
  geom_hline(yintercept = 5,linetype="dashed")+
  xlab("Patients")

# detected genes to detect trasncripts, colored by percent MT
metadata %>%
  ggplot(aes(x=nUMI,y=nGene,color=percentMT)) +
  geom_point(size=0.1) +
  scale_color_gradient(low = "grey",high = "red")+
  stat_smooth()+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()+
  geom_vline(xintercept = 700,linetype="dashed")+
  geom_hline(yintercept = 300,linetype="dashed")+
  facet_wrap(~orig.ident)

# genes/UMI to show the complexity
metadata %>%
  ggplot(aes(x=log10GenesPerUMI,color=orig.ident,fill=orig.ident))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(fill="Patients",color="Patients")+
  geom_vline(xintercept = 0.85,linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20)) +
  ggtitle("Transcriptomic complexity of merged objects")

filtered_immu <- subset(immu,
                         subset = (nGene > 300) &
                           (nGene < 3000) &
                           (nUMI > 700) &
                           (percentMT < 5))

saveRDS(filtered_immu,file = "../05-inte/filtered_immu_raw.rds")

# gene level filtering
counts <- GetAssayData(object = filtered_immu, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes,]
filtered_immu <- CreateSeuratObject(counts = filtered_counts,meta.data = filtered_immu@meta.data)

# filter comparison to the raw obj plot
filtered_immu@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI,color=orig.ident,fill=orig.ident))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(fill="Patients",color="Patients")+
  geom_vline(xintercept = 0.85,linetype="dashed")+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20))+
  ggtitle("Transcriptiomic complexity of filtered immu objects")

dev.off()

saveRDS(filtered_immu,file = "./filtered_immu_obj.rds")
rm(list = ls())
gc()

# integration ----
filtered_immu <- readRDS("./filtered_immu_obj.rds")

# remove IGH and IGL related genes
# counts <- GetAssayData(object = filtered_immu,slot = "counts")
# counts <- counts[-grep("^IGH",rownames(counts)),]
# counts <- counts[-grep("^IGL",rownames(counts)),]
# filtered_immu <- CreateSeuratObject(counts = counts,meta.data = filtered_immu@meta.data)

#* scTransform ----

#** phase contributes little ----
phase_obj <- NormalizeData(filtered_immu)
phase_obj <- CellCycleScoring(phase_obj,
                              g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)
phase_obj <- FindVariableFeatures(phase_obj,
                                  selection.method = "vst",
                                  nfeatures = 2000,
                                  verbose = F)
phase_obj <- ScaleData(phase_obj)
phase_obj <- RunPCA(phase_obj)
DimPlot(phase_obj,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

#** scTransform in immu ----
split_obj <- SplitObject(filtered_immu,split.by = "orig.ident")

for (i in paste0("P",c(1:7))){
  split_obj[[i]] <- NormalizeData(split_obj[[i]],verbose = T)
  split_obj[[i]] <- CellCycleScoring(split_obj[[i]],
                                     g2m.features = cc.genes$g2m.genes,
                                     s.features = cc.genes$s.genes)
  split_obj[[i]] <- SCTransform(split_obj[[i]],vars.to.regress = "percentMT")
}

integ_features <- SelectIntegrationFeatures(object.list = split_obj,
                                            nfeatures = 3000)
split_obj <- PrepSCTIntegration(object.list = split_obj,
                                anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_obj,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
integ_immu <-IntegrateData(anchorset = integ_anchors,
                            normalization.method = "SCT")
saveRDS(integ_immu,file = "./scTrans_immu_obj.rds")

integ_immu <- RunPCA(integ_immu)
integ_immu <- FindNeighbors(object = integ_immu,
                             dims=1:40)
integ_immu <- FindClusters(object = integ_immu,
                            resolution=c(0.6,0.9,1.2,1.5))
integ_immu <- RunUMAP(integ_immu,
                       reduction = "pca",
                       dims = 1:40)
metadata <- integ_immu@meta.data

integ_immu@meta.data <- integ_immu@meta.data[,c(-16,-18)]
colnames(integ_immu@meta.data)[16] <- "integ_all_snn_res.1.8"

saveRDS(integ_immu,file = "scTrans_immu_umap.rds")

DimPlot(integ_immu,
        group.by = "integrated_snn_res.0.9",
        label=T)
DimPlot(integ_immu,
        group.by = "orig.ident",
        split.by = "DrugRes",
        label=T,
        repel = T)
DimPlot(integ_immu,
        group.by = "filter_clonotype",
        label=T,
        repel = T)


