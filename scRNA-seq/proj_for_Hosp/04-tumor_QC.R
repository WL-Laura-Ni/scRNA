# author information----
### QC for subset tumor data
### writer: Wanlin Laura NI
### time: 2023-04-27

# working condition setup----
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
setwd("~/Graduation_proj/03-tumor/")

# data loading ----
tumor <- readRDS("~/Graduation_proj/02-reduction/tumor_raw.rds")
Idents(tumor) <- "orig.ident"

if (FALSE) {
# plot raw data information for QC----
metadata <- as.data.frame(tumor@meta.data)

pdf(file = "./Tumor-QC-plot.pdf",width = 8,height = 6)

# detected cells/patient
metadata %>% 
  ggplot(aes(x=orig.ident,fill=DrugRes))+
  geom_bar() +
  geom_hline(yintercept = 10000,linetype="dashed") +
  geom_hline(yintercept = 12500,linetype="dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20)) +
  xlab("Patients") +
  ylab("Cell Count")+
  ggtitle("Numbers of detected cells")

# nGene of the patients
VlnPlot(tumor,
        features = "nGene",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 4000,linetype="dashed")+
  geom_hline(yintercept = 500,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# nUMI of the patients
VlnPlot(tumor,
        features="nUMI",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 2000,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# percent of MT genes
VlnPlot(tumor,
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
  geom_vline(xintercept = 2000,linetype="dashed")+
  geom_hline(yintercept = 500,linetype="dashed")+
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

filtered_tumor <- subset(tumor,
                subset = (nGene > 500) &
                  (nGene < 4000) &
                  (nUMI > 2000) &
                  (percentMT < 5))

saveRDS(filtered_tumor,file = "../05-inte/filtered_tumor_raw.rds")

# gene level filtering
counts <- GetAssayData(object = filtered_tumor, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes,]
filtered_tumor <- CreateSeuratObject(counts = filtered_counts,meta.data = filtered_tumor@meta.data)

# filter comparison to the raw obj plot
filtered_tumor@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI,color=orig.ident,fill=orig.ident))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(fill="Patients",color="Patients")+
  geom_vline(xintercept = 0.85,linetype="dashed")+
  theme(plot.title = element_text(hjust=0.5,face="bold",size=20))+
  ggtitle("Transcriptiomic complexity of filtered tumor objects")

dev.off()

saveRDS(filtered_tumor,file = "./filtered_tumor_obj.rds")
rm(list = ls())
gc()
}

# integration ----
filtered_tumor <- readRDS("./filtered_tumor_obj.rds")

# remove IGH and IGL related genes
counts <- GetAssayData(object = filtered_tumor,slot = "counts")
counts <- counts[-grep("^IGH",rownames(counts)),]
counts <- counts[-grep("^IGL",rownames(counts)),]
filtered_tumor <- CreateSeuratObject(counts = counts,meta.data = filtered_tumor@meta.data)

#* scTransform ----
if (FALSE) {
#** phase contributes little ----
phase_obj <- NormalizeData(filtered_tumor)
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
rm(phase_obj)
}

#** scTransform in tumor ----
tumor_list <- SplitObject(filtered_tumor,split.by = "orig.ident")
tumor_list <- lapply(X = tumor_list,FUN = function(x){
  x <- NormalizeData(x)
  x <- CellCycleScoring(x,
                        g2m.features = cc.genes$g2m.genes,
                        s.features = cc.genes$s.genes)
  x <- SCTransform(x,vars.to.regress = "percentMT")
})

integ_features <- SelectIntegrationFeatures(object.list = tumor_list,
                                            nfeatures = 3000)
tumor_list <- PrepSCTIntegration(object.list = tumor_list,
                                anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = tumor_list,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
integ_tumor <-IntegrateData(anchorset = integ_anchors,
                            normalization.method = "SCT")
saveRDS(integ_tumor,file = "./scTrans_tumor_obj.rds")

integ_tumor <- RunPCA(integ_tumor)
integ_tumor <- FindNeighbors(object = integ_tumor,
                             dims=1:40)
integ_tumor <- FindClusters(object = integ_tumor,
                            resolution=c(0.6,0.9,1.2,1.5))
integ_tumor <- RunUMAP(integ_tumor,
                       reduction = "pca",
                       dims = 1:40)
metadata <- integ_tumor@meta.data

integ_tumor@meta.data <- integ_tumor@meta.data[,c(-16,-18)]
colnames(integ_tumor@meta.data)[16] <- "integ_all_snn_res.1.8"

saveRDS(integ_tumor,file = "scTrans_tumor_umap.rds")

DimPlot(integ_tumor,
        group.by = "integrated_snn_res.0.9",
        label=T)
DimPlot(integ_tumor,
        group.by = "orig.ident",
        split.by = "DrugRes",
        label=T,
        repel = T)
DimPlot(integ_tumor,
        group.by = "filter_clonotype",
        label=T,
        repel = T)

#* seurat ----
if (FALSE) {
#** integrated by experimental set ----

features <- SelectIntegrationFeatures(object.list = tumor_list)
tumor_anchors <- FindIntegrationAnchors(object.list = tumor_list,anchor.features = features)
tumor_comb <- IntegrateData(anchorset = tumor_anchors)

DefaultAssay(tumor_comb) <- "integrated"
tumor_comb <- ScaleData(tumor_comb,verbose = F)
tumor_comb <- RunPCA(tumor_comb,npcs = 40,verbose=F)
ElbowPlot(tumor_comb,ndims = 30)+
  geom_hline(yintercept = 1.5,linetype="dashed")
tumor_comb <- RunUMAP(tumor_comb,reducation="pca",dims=1:24)
tumor_comb <- FindNeighbors(tumor_comb,reduction="pca",dims = 1:24)
tumor_comb <- FindClusters(tumor_comb,resolution = c(0.6,0.9,1.2,1.5))

saveRDS(tumor_comb,file = "./tumor_set_umap.rds")

#** integrated by sample ----
tumor_list <- SplitObject(filtered_tumor,split.by="orig.ident")
tumor_list <- lapply(X = tumor_list,FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method = "vst",nfeatures=2000)
})

features <- SelectIntegrationFeatures(object.list = tumor_list)
tumor_anchors <- FindIntegrationAnchors(object.list = tumor_list,anchor.features = features)
tumor_comb <- IntegrateData(anchorset = tumor_anchors)

DefaultAssay(tumor_comb) <- "integrated"
tumor_comb <- ScaleData(tumor_comb,verbose = F)
tumor_comb <- RunPCA(tumor_comb,npcs = 40,verbose=F)
ElbowPlot(tumor_comb,ndims = 30)+
  geom_hline(yintercept = 1.5,linetype="dashed")
tumor_comb <- RunUMAP(tumor_comb,reducation="pca",dims=1:16)
tumor_comb <- FindNeighbors(tumor_comb,reduction="pca",dims = 1:16)
tumor_comb <- FindClusters(tumor_comb,resolution = c(0.6,0.9,1.2,1.5))

saveRDS(tumor_comb,file = "./tumor_sample_umap.rds")
}
