# author information----
### Annotation and Plot
### writer: Wanlin Laura NI
### time: 2023-04-11

# working condition setup----
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
setwd("~/Graduation_proj/02-reduction_annotation/")

# data loading ----
obj <- readRDS("~/Graduation_proj/02-reduction_annotation/obj_seurat_umap.rds")
Idents(obj) <- "RNA_snn_res.0.6" 

obj.markers <- FindAllMarkers(obj,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
obj.markers %>%
  group_by(cluster) %>%
  slice_max(n=8,order_by = avg_log2FC)

saveRDS(obj.markers,file = "./obj_markers.rds")

# cell annotation ----
obj$cellType <- "Tumor cells"
obj$cellType[which(obj$RNA_snn_res.0.6 == 7|obj$RNA_snn_res.0.6 == 14|
                     obj$RNA_snn_res.0.6 == 15|obj$RNA_snn_res.0.6 == 17)] <- "Immune cells"
DimPlot(obj,group.by = "cellType",label = T)+NoLegend()
DimPlot(obj,group.by = "cellType",label = T,split.by = "DrugRes")+NoLegend()

# subset tumor cells for DEG ----
obj_tumor <- subset(obj,subset= (cellType == "Tumor cells"))
Idents(obj_tumor) <- "DrugRes"
DrugRes_markers <- FindMarkers(obj_tumor,ident.1 = "Drug_Resistant",
                               ident.2 = "Drug_Unresistant")
DEG_result <- DrugRes_markers
# remove the IG related genes due to the specificity of MM
DEG_result <- DEG_result[-grep("^IG",rownames(DEG_result)),]

# diagonal volcano plot ----
pic_cells <- as.data.frame(log1p(AverageExpression(obj_tumor,verbose = F)$RNA))
pic_cells$genes <- rownames(pic_cells)
pic_up <- rownames(DEG_result)[which(DEG_result$avg_log2FC > 0)]
pic_down <- rownames(DEG_result)[which(DEG_result$avg_log2FC < 0)]
pic_label <- rownames(subset(DEG_result,abs(DEG_result$avg_log2FC)>1.5))
pic_cells$label <- "NoSig"
pic_cells$label[which(pic_cells$gene %in% pic_up)] <- "Up"
pic_cells$label[which(pic_cells$gene %in% pic_down)] <- "Down"
p1 <- ggplot(pic_cells,aes(Drug_Unresistant,Drug_Resistant)) +
  geom_point(size=0.2,aes(color=label)) +
  scale_color_manual(values = c("blue","grey","red")) +
  xlim(0,5) +
  ylim(0,5) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  ggtitle("DEG")

LabelPoints(plot = p1,points = pic_label,repel=T,size=3,segment.size=0.1,segment.alpha=0.6)

# bar plot to show the proportion of immune cells ----
n_cells <- FetchData(obj,vars = c("orig.ident","cellType")) %>%
  dplyr::count(orig.ident,cellType) %>%
  tidyr::spread(orig.ident,n)
n_cells <- reshape2::melt(n_cells,id.vars="cellType",variable.name="numbers")

n_t_cells <- FetchData(obj,vars = c("orig.ident","cellType")) %>%
  dplyr::count(orig.ident,cellType) %>%
  tidyr::spread(cellType,n)
n_t_cells$cells <- n_t_cells$`Immune cells`+n_t_cells$`Tumor cells`                                                                                           
n_t_cells$immu_prop <- signif(n_t_cells$`Immune cells`/n_t_cells$cells*100,4)

n_cells$cells <- 0
n_cells$immu_prop <- 0
n_cells$cells[c(1,3,5,7,9,11,13)] <- n_t_cells$cells
n_cells$cells[c(2,4,6,8,10,12,14)] <- n_t_cells$cells
n_cells$immu_prop[c(2,4,6,8,10,12,14)] <- n_t_cells$immu_prop
n_cells$immu_prop[c(1,3,5,7,9,11,13)] <- n_t_cells$immu_prop


ggplot(n_cells)+
  geom_bar(aes(x=numbers,y=value,fill=cellType),
           stat = "identity")+
  geom_text(aes(x=numbers,y=cells,label=paste0(immu_prop,"%")),vjust=-0.5)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=15))+
  xlab("Patient")+
  ylab("Cell Counts")+
  ggtitle("Immune cell proportion in patients")

gc()

