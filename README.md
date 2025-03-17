# Aging Mouse Colon ScRNA-seq Analysis
This is a single-cell RNA seq analysis workflow using a dataset of epithelial cells from the mouse colon (can be found in the GEO repository with accession numbers GSE168448). 
The objective of this analysis is to find genes which is differentially expressed in aging colon compared to a young colon. 

## Contents of the repository

- The R script `Colon scRNAseq merged.R` that was used for this analysis.

## Table of contents

- [Step 0: Set up](#step-0-set-up)
- [Step 1: Read Input files](#step-1-read-input-files)
- [Step 2: Filtering](#step-2-filtering)
- [Step 3: Normalization](#step-3-normalization)
- [Step 4: Add aneuploidy information](#step-4-Find-Highly-Variable-Genes)
- [Step 5: Higly variable genes](#step-5-Scale-Data)
- [Step 6: Scale data](#step-6-Dimensional-Reduction)
- [Step 7: Dimensionality reduction](#step-7-Find-Cluster-Markers)
- [Step 8: Cluster markers](#step-8-Merge-Clusters-10-and-12-and-Focus-on-this-Clsuter)
- [Step 9: Scores](#step-9-Enrichment-Analysis)


## Pipeline
### Step 0: Set up

```
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggrepel)
```

### Step 1: Read Input Files
```
files <- list.files(path = "/Users/chiarasantoso/Desktop/GSE168448_RAW", pattern = "*.csv", full.names = TRUE)
```
Separate files by age group (young vs old)
```
young_files <- files[grepl("Young", files)]
old_files <- files[grepl("Old", files)]

young_data <- do.call(cbind, lapply(young_files, read.csv, row.names = 1))
colnames(young_data) <- paste0("Young_", colnames(young_data)) 

old_data <- do.call(cbind, lapply(old_files, read.csv, row.names = 1))
colnames(old_data) <- paste0("Old_", colnames(old_data))
```
Merge young and old
```
merged_data <- merge(young_data, old_data, by = "row.names", all = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[, -1]
```
Create seurat object
```
merged_seurat <- CreateSeuratObject(counts = merged_data)
```


### Step 2: Filtering
```
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")
```
Visualize QC metrics
```
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
Filtering based on QC metrics
```
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 10)
```

### Step 3: Normalization
```
merged_seurat <- NormalizeData(merged_seurat)
```


### Step 4: Find Highly Variable Genes
```
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)
```

### Step 5: Scale Data
```
all.genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes)
```

### Step 6: Dimensional Reduction

Principal Component Analysis (PCA)
```
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))
```
Cluster cells
```
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:12)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)
```
UMAP
```
merged_seurat <- RunUMAP(merged_seurat, dims = 1:12)

merged_seurat$group <- ifelse(grepl("Young", merged_seurat$orig.ident), "Young", "Old")
Idents(merged_seurat) <- "group"
```
Plot the distribution of young cells vs old cells
```
DimPlot(merged_seurat, reduction = "umap", group.by = "group", pt.size = 0.5, alpha = 0.5 ) +
  ggtitle("Young vs Old")
```

Plot the clusters
```
DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)
 ```
Plot the QC metrics on the UMAP
```
FeaturePlot(merged_seurat, features = "nFeature_RNA", reduction = "umap") +
  ggtitle("nFeature_RNA")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
FeaturePlot(merged_seurat, features = "nCount_RNA", reduction = "umap") +
  ggtitle("nCount_RNA")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
FeaturePlot(merged_seurat, features = "percent.mt", reduction = "umap") + 
  ggtitle("percent.mt")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
```

Calculate percentage of old/young in each cluster 
```
metadata <- merged_seurat@meta.data
cluster_group_counts <- metadata %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(count = n()) %>%
  ungroup()
cluster_group_percent <- cluster_group_counts %>%
  group_by(seurat_clusters) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()
```

Bar plot of Percentage of Young vs. Old Cells per Cluster 
```
ggplot(cluster_group_percent, aes(x = seurat_clusters, y = percent, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Percentage", fill = "Group", 
       title = "Percentage of Young vs. Old Cells per Cluster") +
  theme_minimal()
```


calculate mitochondrial percentage for each cluster (to make sure no clusters represents dying cells)
```
mt_percent_by_cluster <- metadata %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_mt_percent = mean(percent.mt, na.rm = TRUE),
    median_mt_percent = median(percent.mt, na.rm = TRUE)
  )

print(mt_percent_by_cluster)
ggplot(mt_percent_by_cluster, aes(x = seurat_clusters, y = mean_mt_percent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Mean Mitochondrial Percentages by Cluster", x = "Cluster", y = "Mean MT Percentage")
```


### Step 7: Find Cluster Markers
```
merged.markers <- FindAllMarkers(merged_seurat)
```
Rank based on log2fc
```
merged.markers <- merged.markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(merged.markers, "merged_cluster_markers.csv")
merged.markers <- read.csv("merged_cluster_markers.csv")
```

Overlays
```
VlnPlot(merged_seurat, features = c("Lpar1"), split.by = "orig.ident")
FeaturePlot(merged_seurat, features = c("Sox2", "Lgr5", "Aldh1", "Epcam", "Mki67"))
```

After visualizing overlays of certain genes, we hypothesized that clusters 10 and 12 were stem cells and we decided to focus on these clusters.



### Step 8: Merge Clusters 10 and 12, and Focus on this Clsuter

Combine clusters 10 and 12 to create one stem cell cluster
```
merged_seurat$combined_cluster <- as.character(merged_seurat$seurat_clusters)
merged_seurat$combined_cluster[merged_seurat$seurat_clusters %in% c("10", "12")] <- "Stem_cell_cluster"
merged_seurat$combined_cluster <- factor(merged_seurat$combined_cluster) 
merged_seurat@meta.data$original_cluster <- merged_seurat$seurat_clusters
Stem_cell_cluster <- subset(merged_seurat, subset = combined_cluster == "Stem_cell_cluster")
```
Find markers for the stem cell cluster 
```
StemCells_YvsO_markers <- FindMarkers(cluster_10_12, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")
```
7871 markers were found. Only keep markers with p_val <= 0.05, pct.1 > 0.1 and pct.2 > 0.1, and arrange it by abs_avg_log2FC
```
StemCells_YvsO_markers <- StemCells_YvsO_markers %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

dim(StemCells_YvsO_markers) # 2230 markers remain after filtering

StemCells_YvsO_markers$gene <- rownames(StemCells_YvsO_markers)
```

Visualize the top 20 markers of the stem cell cluster 
```
top20 <- StemCells_YvsO_markers %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(StemCells_YvsO_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 10 + 12")
```
Visualize overlays of certain markers within the stem cell cluster 
```
FeaturePlot(Stem_cell_cluster, features = c("Hspa1a"))
VlnPlot(Stem_cell_cluster, features = c("Ano6"), group.by = "original_cluster", split.by = "orig.ident")
```

### Step 9: Enrichment Analysis

Required libraries
```
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
```

Select Top 200 Genes
```
top_genes <- StemCells_YvsO_markers %>%
  head(200)

gene_list <- top_genes$gene
```

Map Gene Symbols to Entrez IDs
```
gene_entrez <- bitr(gene_list, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
gene_entrez <- gene_entrez[!duplicated(gene_entrez$ENTREZID), ]
```

Perform GO Enrichment Analysis
```
ego <- enrichGO(gene         = gene_entrez$ENTREZID,
                OrgDb        = org.Mm.eg.db,
                keyType      = "ENTREZID",
                ont          = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable     = TRUE)


ego_df <- as.data.frame(ego)
```

Visualize GO Results
```
barplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment")

```

Perform KEGG Enrichment Analysis
```
ekegg <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                    organism     = 'mmu',
                    keyType      = 'kegg',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 0.5)

ekegg_df <- as.data.frame(ekegg)
```

Visualize KEGG Results
```
barplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")
```

