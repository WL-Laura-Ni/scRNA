# author information----
### QC for raw data
### writer: Wanlin Laura NI
### time: 2023-03-26

# working condition setup----
library(dplyr)
library(cowplot)
library(ggplot2)
library(Seurat)
library(patchwork)
setwd("~/Graduation_proj/01-QC/")

# pre quality control by sample ----
for (patient in c(paste0("P",c(1:7)))){
  # data loading----
  bmmc.data <- Read10X(data.dir = paste0("~/MultipleMyeloma/Sample/",patient,"/filtered_feature_bc_matrix/"))
  bmmc <- CreateSeuratObject(counts = bmmc.data,
                             project = patient,
                             min.cells = 3, 
                             min.features = 100)
  
  # loading BCR seq information
  BCR_data <- read.csv(file = paste0("~/MultipleMyeloma/Sample/",patient,"/filtered_contig_annotations.csv"))
  BCR_data <- BCR_data[!duplicated(BCR_data$barcode),]
  names(BCR_data)[names(BCR_data) == "raw_clonotype_id"] <- "clonotype_id"
  BCR_clone <- read.csv(file = paste0("~/MultipleMyeloma/Sample/",patient,"/clonotypes.csv"))
  BCR_data <- merge(BCR_data,BCR_clone)
  rownames(BCR_data) <- BCR_data[,2]
  bmmc <- AddMetaData(object = bmmc,metadata = BCR_data)
  
  saveRDS(bmmc,file = paste0("./BCR_add_raw_obj_",patient,".rds"))
  
  # BCR information is reduced here to decrease the memory size
  bmmc@meta.data <- bmmc@meta.data[,1:4]
  
  assign(patient,bmmc)
}

# merge objects and remove unnecessary variables ----
merged_obj <- merge(x=P1,
                    y=list(P2,P3,P4,P5,P6,P7),
                    add.cell.ids = c(paste0("P",c(1:7))))

rm(list = setdiff(ls(),ls(pattern="merged_obj")))
gc()

merged_obj$log10GenesPerUMI <- log10(merged_obj$nFeature_RNA)/log10(merged_obj$nCount_RNA)
merged_obj$percentMT <- PercentageFeatureSet(object = merged_obj,patter="^MT-")

merged_obj$DrugRes <- "NA"
merged_obj$DrugRes[which(merged_obj$orig.ident == "P1" |
                         merged_obj$orig.ident == "P2" |
                         merged_obj$orig.ident == "P3")] <- "Drug_Unresistant"
merged_obj$DrugRes[which(merged_obj$orig.ident == "P4" |
                         merged_obj$orig.ident == "P5" |
                         merged_obj$orig.ident == "P6" |
                         merged_obj$orig.ident == "P7" )] <- "Drug_Resistant"
colnames(merged_obj@meta.data)[2:3] <- c("nUMI","nGene") 
saveRDS(object = merged_obj,file = "./raw_merged_obj.rds")

# plot raw data information for QC----
merged_obj <- readRDS(file = "raw_merged_obj.rds")
metadata <- as.data.frame(merged_obj@meta.data)

pdf(file = "./QC-plot.pdf",width = 8,height = 6)

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
VlnPlot(merged_obj,
        features = "nGene",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 4000,linetype="dashed")+
  geom_hline(yintercept = 300,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# nUMI of the patients
VlnPlot(merged_obj,
        features="nUMI",
        group.by = "orig.ident",
        pt.size = 0)+
  geom_hline(yintercept = 2000,linetype="dashed")+
  theme(legend.position = "None",plot.title = element_text(size = 20))+
  xlab("Patients")

# percent of MT genes
VlnPlot(merged_obj,
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


# cell level filtering
filtered_obj <- subset(merged_obj,
                     subset= (nGene > 300 ) &
                       (nGene < 4000) &
                       (nUMI > 2000) &
                       (percentMT < 5))

# gene level filtering
counts <- GetAssayData(object = filtered_obj,slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes,]
filtered_obj <- CreateSeuratObject(counts = filtered_counts,meta.data = filtered_obj@meta.data)

# filter comparison to the raw obj plot
filtered_obj@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI,color=orig.ident,fill=orig.ident))+
  geom_density(alpha=0.2)+
  theme_classic()+
  labs(fill="Patients",color="Patients")+
  geom_vline(xintercept = 0.85,linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 20)) +
  ggtitle("Transcriptomic complexity of filtered merged objects")

# data frame construction
compare_mt <- data.frame("1","2","3","4","5")
colnames(compare_mt) <- c("object","cell_numbers","cell_prop","gene_numbers","gene_prop")
compare_mt[1,] <- c("merged_obj",ncol(merged_obj),NA,nrow(merged_obj),NA)
compare_mt[2,] <- c("filtered_obj",ncol(filtered_obj),NA,nrow(filtered_obj),NA)
compare_mt[,1] <- factor(compare_mt[,1],levels=c("merged_obj","filtered_obj"))
compare_mt[,2] <- as.numeric(compare_mt[,2])
compare_mt[,3] <- signif(as.numeric(compare_mt[,2]/compare_mt[1,2])*100,4)
compare_mt[,4] <- as.numeric(compare_mt[,4])
compare_mt[,5] <- signif(as.numeric(compare_mt[,4]/compare_mt[1,4])*100,4)

compare_mt %>%
  ggplot()+
  geom_bar(aes(x=object,y=cell_numbers,fill=object),
           stat="identity")+
  geom_text(aes(x=object,y=cell_numbers,label=paste0(cell_prop,"%")),vjust=-0.2,size=5)+
  xlab("object")+
  ylab("Cell numbers")+
  theme_classic()+
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5,face="bold",size=13))+
  ggtitle("Cell number change after the filter") -> p1

compare_mt %>%
  ggplot()+
  geom_bar(aes(x=object,y=gene_numbers,fill=object),
           stat="identity")+
  geom_text(aes(x=object,y=gene_numbers,label=paste0(gene_prop,"%")),vjust=-0.2,size=5)+
  xlab("object")+
  ylab("Gene numbers")+
  theme_classic()+
  theme(legend.position = "None",
        plot.title = element_text(hjust = 0.5,face="bold",size=13))+
  ggtitle("Gene number change after the filter") ->p2

p1+p2

dev.off()

saveRDS(filtered_obj,file = "./filtered_merged_obj.rds")

rm(list = ls())
gc()
