# scRNA-seq_Cancer
Meta analysis of pancreatic cancer scRNA-seq datasets using Seurat and Harmony

# Data Sources
I downloaded count data from the following sources and loaded as dataframes in R
```
Oh et al., Nat Commun 2023
Qadir et al., 2020
Lin et al., 2020
Moncada et al., 2020
```
# Analysis steps:
(1) Merged the datasets and removed mitochondrial and doublet cells
```
    subset = nFeature_RNA > 300 & nCount_RNA < 50000 & percent.mt < 10
```
(2) Normalized and scaled the data, identified clusters using Seurat
```
pancreas_SWF <- NormalizeData(pancreas) 
pancreas_SWF <- FindVariableFeatures(pancreas_SWF) 
pancreas_SWF <- ScaleData(pancreas_SWF) 
pancreas_SWF <- RunPCA(pancreas_SWF, verbose = FALSE)
ElbowPlot(pancreas_SWF, ndims = 50, reduction = "pca")
pancreas_SWF <- RunUMAP(pancreas_SWF, dims = 1:50)
pancreas_SWF <- FindNeighbors(pancreas_SWF, dims = 1:50) 
pancreas_SWF <- FindClusters(pancreas_SWF)
```
(3) Used Harmony to merge the datasets as much as possible
```
pancreas_HAR <- RunUMAP(pancreas_HAR, reduction = "harmony", dims = 1:50)
```
Note: I tried different resolutions to see which one made sense using the visual PCA plot

(4) Nominate normal and tumor cells
I did not have a matched normal for every tumor sample (which is ideal). However, Qadir et al., 2020 dataset is of normal pancreatic cells hence labeled them as normal with the rest as tumor. 
```
normalcells <- subset(pancreas_HAR, subset = orig.ident %in% c("Qadir_hpd1","Qadir_hpd2"))
## Get cell names
normalcellNames <- rownames(normalcells@meta.data)
## Mutate a column in original metadata
pancreas_HAR$condition<- ifelse(colnames(pancreas_HAR) %in% normalcellNames, "Normal", "Tumor")
```

(5) Identify gene markers for tumor versus normal cells within duct cell cluster 
To identify the cluster corresponding to duct cells, I visualized the expression of KRT19 and SOX9 genes.
```
FeaturePlot(pancreas_HAR, features = c("KRT19", "SOX9"))
```
I found it to be cluster 0 so now finding markers within cluster0 which ditinguish between normal and tumor cells.
```
pancreas_HAR$cluster.cond <- paste(pancreas_HAR$seurat_clusters, pancreas_HAR$condition, sep = "_")
Idents(pancreas_HAR) <- "cluster.cond"
top10_tum_norm_pos_clust0<- head(FindMarkers(pancreas_HAR, ident.1 = "0_Tumor", ident.2 = "0_Normal",only.pos = TRUE, verbose = FALSE,features=top8000),10)
```





