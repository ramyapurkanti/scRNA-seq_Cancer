library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tibble)
library(magrittr)

setwd("Cancer_scRNAseq_pancreas/Data")

#Load Oh et al., Nat Commun 2023 GSE231535
data <- Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE231535_RAW/PDAC1")
Oh_PDAC1 <- CreateSeuratObject(counts = data, project = "Oh_pdac1", min.cells = 3, min.features = 200)
data <- Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE231535_RAW/PDAC2")
Oh_PDAC2 <- CreateSeuratObject(counts = data, project = "Oh_pdac2", min.cells = 3, min.features = 200)
Oh<-merge(Oh_PDAC1,y=c(Oh_PDAC2))
rm(Oh_PDAC1,Oh_PDAC2)
# Load the data from Qadir et al., 2020
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE131886_RAW/HPD1/")
Qadir_HPD1<- CreateSeuratObject(counts = data, project = "Qadir_hpd1", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE131886_RAW/HPD2/")
Qadir_HPD2<- CreateSeuratObject(counts = data, project = "Qadir_hpd2", min.cells = 3, min.features = 200)
#HPD3 raw matrix file is corrupted
#data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE131886_RAW/HPD3/")
#Qadir_HPD3<- CreateSeuratObject(counts = data, project = "Qadir_hpd3", min.cells = 3, min.features = 200)
Qadir<-merge(Qadir_HPD1,y=c(Qadir_HPD2))
rm(Qadir_HPD1,Qadir_HPD2)
# Load the data from Lin et al., 2020
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T2/")
Lin_T2<- CreateSeuratObject(counts = data, project = "Lin_t2", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T3/")
Lin_T3<- CreateSeuratObject(counts = data, project = "Lin_t3", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T4/")
Lin_T4<- CreateSeuratObject(counts = data, project = "Lin_t4", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T5/")
Lin_T5<- CreateSeuratObject(counts = data, project = "Lin_t5", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T6/")
Lin_T6<- CreateSeuratObject(counts = data, project = "Lin_t6", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T8/")
Lin_T8<- CreateSeuratObject(counts = data, project = "Lin_t8", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T9/")
Lin_T9<- CreateSeuratObject(counts = data, project = "Lin_t9", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/T10/")
Lin_T10<- CreateSeuratObject(counts = data, project = "Lin_t10", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/K16733/")
Lin_K16733<- CreateSeuratObject(counts = data, project = "Lin_k16733", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00006/")
Lin_Y00006<- CreateSeuratObject(counts = data, project = "Lin_y00006", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00008/")
Lin_Y00008<- CreateSeuratObject(counts = data, project = "Lin_y00008", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00013/")
Lin_Y00013<- CreateSeuratObject(counts = data, project = "Lin_y00013", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00014/")
Lin_Y00014<- CreateSeuratObject(counts = data, project = "Lin_y00014", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00016/")
Lin_Y00016<- CreateSeuratObject(counts = data, project = "Lin_y00016", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00019/")
Lin_Y00019<- CreateSeuratObject(counts = data, project = "Lin_y00019", min.cells = 3, min.features = 200)
data <-Read10X(data.dir = "Cancer_scRNAseq_pancreas/Data/GSE154778_RAW/Y00027/")
Lin_Y00027<- CreateSeuratObject(counts = data, project = "Lin_y00027", min.cells = 3, min.features = 200)
Lin<-merge(Lin_T2,y=c(Lin_T3,Lin_T4,Lin_T5,Lin_T6,Lin_T8,Lin_T9,Lin_T10,Lin_K16733,Lin_Y00006,Lin_Y00008,Lin_Y00013,Lin_Y00014,Lin_Y00016,Lin_Y00019,Lin_Y00027))
rm(Lin_T2,Lin_T3,Lin_T4,Lin_T5,Lin_T6,Lin_T8,Lin_T9,Lin_T10,Lin_K16733,Lin_Y00006,Lin_Y00008,Lin_Y00013,Lin_Y00014,Lin_Y00016,Lin_Y00019,Lin_Y00027)
#Lin<-merge(Lin_T2,y=c(Lin_T3,Lin_T4,Lin_T5,Lin_T6,Lin_T8,Lin_T9,Lin_T10,Lin_Y00008,Lin_Y00013,Lin_Y00014,Lin_Y00016,Lin_Y00019,Lin_Y00027))
#rm(Lin_T2,Lin_T3,Lin_T4,Lin_T5,Lin_T6,Lin_T8,Lin_T9,Lin_T10,Lin_Y00008,Lin_Y00013,Lin_Y00014,Lin_Y00016,Lin_Y00019,Lin_Y00027)
# Load the data from Moncada et al., 2020
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405527_PDAC-A-indrop3.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACA_S3<- CreateSeuratObject(counts = data, project = "Moncada_PDACA_3", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405528_PDAC-A-indrop4.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACA_S4<- CreateSeuratObject(counts = data, project = "Moncada_PDACA_4", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405529_PDAC-A-indrop5.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACA_S5<- CreateSeuratObject(counts = data, project = "Moncada_PDACA_5", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405530_PDAC-A-indrop6.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACA_S6<- CreateSeuratObject(counts = data, project = "Moncada_PDACA_6", min.cells = 3, min.features = 200)

data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405531_PDAC-B-indrop1.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACB_S1<- CreateSeuratObject(counts = data, project = "Moncada_PDACB_1", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405532_PDAC-B-indrop2.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACB_S2<- CreateSeuratObject(counts = data, project = "Moncada_PDACB_2", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM3405533_PDAC-B-indrop3.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACB_S3<- CreateSeuratObject(counts = data, project = "Moncada_PDACB_3", min.cells = 3, min.features = 200)

data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM4100717_PDAC-C-indrop1.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACC_S1<- CreateSeuratObject(counts = data, project = "Moncada_PDACC_1", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM4100718_PDAC-C-indrop2.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACC_S2<- CreateSeuratObject(counts = data, project = "Moncada_PDACC_2", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM4100719_PDAC-C-indrop3.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACC_S3<- CreateSeuratObject(counts = data, project = "Moncada_PDACC_3", min.cells = 3, min.features = 200)
data <-data.table::fread("Cancer_scRNAseq_pancreas/Data/GSE111672_RAW/GSM4100720_PDAC-C-indrop4.tsv.gz")
data <- data %>% column_to_rownames("Genes")
colnames(data) = gsub("_", "-", colnames(data))
Moncada_PDACC_S4<- CreateSeuratObject(counts = data, project = "Moncada_PDACC_4", min.cells = 3, min.features = 200)
Moncada<-merge(Moncada_PDACA_S3,y=c(Moncada_PDACA_S4,Moncada_PDACA_S5,Moncada_PDACA_S6,Moncada_PDACB_S1,Moncada_PDACB_S2,Moncada_PDACB_S3,Moncada_PDACC_S1,Moncada_PDACC_S2,Moncada_PDACC_S3,Moncada_PDACC_S4))
rm(Moncada_PDACA_S3,Moncada_PDACA_S4,Moncada_PDACA_S5,Moncada_PDACA_S6,Moncada_PDACB_S1,Moncada_PDACB_S2,Moncada_PDACB_S3,Moncada_PDACC_S1,Moncada_PDACC_S2,Moncada_PDACC_S3,Moncada_PDACC_S4)

# Loading data for Muraro et al., 2016 (Transcriptome Atlas)
#4 donors, 8 cultures per donor so 32 samples
#There are fractional values suggesting that the data has already been normalized
#data <- data.table::fread("GSE85241_RAW/GSE85241_cellsystems_dataset_4donors_updated.csv")
# Gene names are in the first column so we need to move them to rownames
#data <- data %>% column_to_rownames("V1")
#GSE85241 <- CreateSeuratObject(counts = data, project = "Muraro", min.cells = 3, min.features = 200)

# Load the data from Segerstolpe et al., 2016
#The first two columns are gene names and then 3514 columns are RPKM values (ignore) and the last 3514 columns are the counts
#data <- data.table::fread("EMTAB-5061/pancreas_refseq_rpkms_counts_3514sc.txt")
# Gene names are in the first column so we need to move them to rownames
#data <- data %>% column_to_rownames("#samples")
#EMTAB_5061 <- CreateSeuratObject(counts = data, project = "Segerstolpe", min.cells = 3, min.features = 200)

pancreas<-merge(Oh,y=c(Qadir,Lin,Moncada))

pancreas <- PercentageFeatureSet(pancreas, pattern = "^MT-", col.name = "percent.mt")
view(pancreas@meta.data)
VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
pancreas <- subset(pancreas, subset = nFeature_RNA > 300 & nCount_RNA < 50000 & percent.mt < 10)

# Seurat Standard work flow (SWF)
pancreas_SWF <- NormalizeData(pancreas) 
pancreas_SWF <- FindVariableFeatures(pancreas_SWF) 
pancreas_SWF <- ScaleData(pancreas_SWF) 
pancreas_SWF <- RunPCA(pancreas_SWF, verbose = FALSE)
ElbowPlot(pancreas_SWF, ndims = 50, reduction = "pca")
pancreas_SWF <- RunUMAP(pancreas_SWF, dims = 1:50)
pancreas_SWF <- FindNeighbors(pancreas_SWF, dims = 1:50) 
pancreas_SWF <- FindClusters(pancreas_SWF)
FeaturePlot(pancreas_SWF, features = c("GCG", "INS", "SST", "PPY","PRSS1","KRT19","COL1A1","GHRL","ESAM"))
DimPlot(pancreas_SWF, group.by = c("orig.ident", "seurat_clusters"), ncol = 2)


pancreas_SWF.list <- SplitObject(pancreas_SWF, split.by = 'orig.ident')

# Run_Harmony
library(harmony)

pancreas_HAR <- NormalizeData(pancreas)
pancreas_HAR <- FindVariableFeatures(pancreas_HAR) 
# Identify the 8000 most highly variable genes
top8000 <- head(VariableFeatures(pancreas_HAR), 8000)
pancreas_HAR <- ScaleData(pancreas_HAR) 
pancreas_HAR <- RunPCA(pancreas_HAR, verbose = FALSE, features=top8000)
pancreas_HAR <- RunHarmony(pancreas_HAR, group.by.vars = "orig.ident",plot_convergence=TRUE,max_iter=50)
ElbowPlot(pancreas_HAR, ndims = 50, reduction = "harmony")
pancreas_HAR <- RunUMAP(pancreas_HAR, reduction = "harmony", dims = 1:50)
pancreas_HAR <- FindNeighbors(pancreas_HAR, reduction = "harmony", dims = 1:50) 
pancreas_HAR <- FindClusters(pancreas_HAR,resolution=0.05)
DimPlot(pancreas_HAR, group.by = c("orig.ident", "seurat_clusters"), ncol = 2)
DimPlot(pancreas_HAR, group.by = "seurat_clusters",label=T) + NoLegend()
DimPlot(pancreas_HAR, group.by = "orig.ident") + NoLegend()
pancreas_HAR.markers <- FindAllMarkers(pancreas_HAR,only.pos = TRUE,features=top8000)
pancreas_HAR.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pancreas_HAR, features = top10$gene) + NoLegend()
#library(clustermole)
## Make subset of normal cells as from Qadir dataset
normalcells <- subset(pancreas_HAR, subset = orig.ident %in% c("Qadir_hpd1","Qadir_hpd2"))
## Get cell names
normalcellNames <- rownames(normalcells@meta.data)
## Mutate a column in original metadata
pancreas_HAR$condition<- ifelse(colnames(pancreas_HAR) %in% normalcellNames, "Normal", "Tumor")
DimPlot(pancreas_HAR, label=F, group.by="orig.ident", cells.highlight= list(normalcellNames), cols.highlight = c("grey"), cols= "darkred")

FeaturePlot(pancreas_HAR, features = c("EPCAM", "AMBP", "MUC1", "COL1A1","PECAM1","RGS5","AIF1","MS4A1","CD3D"))
VlnPlot(pancreas_HAR, features = c("EPCAM", "AMBP", "MUC1", "COL1A1","PECAM1","RGS5","AIF1","MS4A1","CD3D"))
#Canonical markers I found in the markers list
FeaturePlot(pancreas_HAR, features = c("KRT19", "SOX9","AMBP","COL1A1", "ESAM", "CD3D","PRSS1","RGS5","GCG"))
new.cluster.ids <- c("Duct cells", "Mesenchyme/Fibroblasts", "2", "Endothelial cells","4", "T-lymphocytes","Acinar", "Pancreatic Stellate cells", "Alpha cells","9","10","11")
names(new.cluster.ids) <- levels(pancreas_HAR)
pancreas_HAR <- RenameIdents(pancreas_HAR, new.cluster.ids)
DimPlot(pancreas_HAR, label = TRUE, pt.size = 0.5) + NoLegend()

#Now I will do differential gene expression analysis between tumor and normal within duct cell cluster (cluster 0)
pancreas_HAR$cluster.cond <- paste(pancreas_HAR$seurat_clusters, pancreas_HAR$condition, sep = "_")
Idents(pancreas_HAR) <- "cluster.cond"
top10_tum_norm_pos_clust0<- head(FindMarkers(pancreas_HAR, ident.1 = "0_Tumor", ident.2 = "0_Normal",only.pos = TRUE, verbose = FALSE,features=top8000),10)

VlnPlot(pancreas_HAR, features <- c('ATP5J2','TFF3'), idents = c("0_Tumor", "0_Normal"), group.by = "condition") 
#Doing pseudo bulking
# pseudobulk the counts based on donor-condition-celltype
pseudo_pancreas <- AggregateExpression(pancreas_HAR, assays = "RNA", return.seurat = T, group.by = c("seurat_clusters","condition","orig.ident"))
# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo_pancreas_HAR$seurat_clusters <- sapply(strsplit(Cells(pseudo_pancreas_HAR), split = "_"), "[", 1)
pseudo_pancreas_HAR$condition <- sapply(strsplit(Cells(pseudo_pancreas_HAR), split = "_"), "[", 2)
pseudo_pancreas_HAR$sample <- paste(sapply(strsplit(Cells(pseudo_pancreas_HAR), split = "_"), "[", 3),sapply(strsplit(Cells(pseudo_pancreas_HAR), split = "_"), "[", 4),sep="_")
pseudo_pancreas_HAR$seurat_clusters.condition <- paste(pseudo_pancreas_HAR$seurat_clusters, pseudo_pancreas_HAR$condition, sep = "_")

Idents(pseudo_pancreas_HAR) <- "seurat_clusters.condition"
bulk.duct.de <- FindMarkers(object = pseudo_pancreas_HAR, 
                            ident.1 = "0_Tumor", 
                            ident.2 = "0_Normal",
                            test.use = "DESeq2",
                            min.cells.group = 2)
head(bulk.duct.de, n = 15)

# Examine and visualize PCA results a few different ways
print(pancreas_HAR[["harmony"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pancreas_HAR, dims = 1:2, reduction = "harmony")

DimPlot(pancreas_HAR, reduction = "harmony") + NoLegend()
DimHeatmap(pancreas_HAR, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pancreas_HAR, dims = 1:20, cells = 500, balanced = TRUE)
# Look at cluster IDs of the first 5 cells
head(Idents(pancreas_HAR), 5)
pbmc <- RunUMAP(pancreas_HAR, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pancreas_HAR, reduction = "umap")
# find all markers of cluster 2
cluster2016.markers <- FindMarkers(pancreas_HAR, ident.1 = c(20,16))
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pancreas_HAR, ident.1=5, ident.2 = c(0, 3))
tumor.markers <- FindMarkers(pancreas_HAR, condition = "tumor")
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
cluster0.markers <- FindMarkers(pancreas_HAR, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
