library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggrepel)
files <- list.files(path = "/Users/chiarasantoso/Desktop/GSE168448_RAW", pattern = "*.csv", full.names = TRUE)

# separate by age group
young_files <- files[grepl("Young", files)]
old_files <- files[grepl("Old", files)]

young_data <- do.call(cbind, lapply(young_files, read.csv, row.names = 1))
colnames(young_data) <- paste0("Young_", colnames(young_data)) 

old_data <- do.call(cbind, lapply(old_files, read.csv, row.names = 1))
colnames(old_data) <- paste0("Old_", colnames(old_data))

#merge young and old
merged_data <- merge(young_data, old_data, by = "row.names", all = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[, -1]
write.csv(merged_data, "/Users/chiarasantoso/Desktop/merged_data.csv")

#create seurat object
merged_seurat <- CreateSeuratObject(counts = merged_data)
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")

#QC
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filter
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 7500 & percent.mt < 10)

VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)

#find variable features
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered, features = all.genes)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:12)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:12)
merged_seurat_filtered <- RunTSNE(merged_seurat_filtered, dims = 1:12)
saveRDS(merged_seurat_filtered, file = "/Users/chiarasantoso/Desktop/merged_seurat_filtered.rds")
saveRDS(merged_seurat, file = "/Users/chiarasantoso/Desktop/merged_seurat.rds")

merged_seurat_filtered$group <- ifelse(grepl("Young", merged_seurat_filtered$orig.ident), "Young", "Old")
Idents(merged_seurat_filtered) <- "group"

#plots
DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "group", pt.size = 0.5, alpha = 0.5 ) +
  ggtitle("Young vs Old")
DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "seurat_clusters", split.by = "group", pt.size = 0.5, label = TRUE) +
  ggtitle("t-SNE Plot by Group (Young vs Old)")

DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)
        
        
FeaturePlot(merged_seurat_filtered, features = "nFeature_RNA", reduction = "umap") +
  ggtitle("nFeature_RNA")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
FeaturePlot(merged_seurat_filtered, features = "nCount_RNA", reduction = "umap") +
  ggtitle("nCount_RNA")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
FeaturePlot(merged_seurat_filtered, features = "percent.mt", reduction = "umap") + 
  ggtitle("percent.mt")+
  scale_colour_gradientn(colours = c("blue", "green", "yellow", "red"))
                         
#calculate percentage of old/young in each cluster
metadata <- merged_seurat_filtered@meta.data
cluster_group_counts <- metadata %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(count = n()) %>%
  ungroup()
cluster_group_percent <- cluster_group_counts %>%
  group_by(seurat_clusters) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()

#bar plot Percentage of Young vs. Old Cells per Cluster
ggplot(cluster_group_percent, aes(x = seurat_clusters, y = percent, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Percentage", fill = "Group", 
       title = "Percentage of Young vs. Old Cells per Cluster") +
  theme_minimal()

# make a plot for each cluster
clusters <- levels(merged_seurat_filtered$seurat_clusters)
for (cluster in clusters) {
  cluster_subset <- subset(merged_seurat_filtered, subset = seurat_clusters == cluster)
  p <- DimPlot(cluster_subset, reduction = "umap", group.by = "orig.ident", alpha = 0.5) +
    ggtitle(paste("Cluster", cluster, "Young vs Old"))
  ggsave(filename = paste0("Cluster_", cluster, "_Young_vs_Old_umap.png"), plot = p, width = 6, height = 5)
}

#calculate mitochondrial percentage for each cluster 
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

#Find markers
merged.markers <- FindAllMarkers(merged_seurat_filtered)

#rank based on log2fc
merged.markers <- merged.markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(merged.markers, "merged_cluster_markers.csv")
merged.markers <- read.csv("merged_cluster_markers.csv")

#overlays
VlnPlot(merged_seurat_filtered, features = c("Lpar1"), split.by = "orig.ident")
FeaturePlot(merged_seurat_filtered, features = c("Sox2", "Lgr5", "Aldh1", "Epcam", "Mki67"))

#annotation
new.cluster.ids <- c("", " ", "", "", "", "","", "", "")
names(new.cluster.ids) <- levels(merged_seurat_filtered)
merged_seurat_filtered <- RenameIdents(merged_seurat_filtered, new.cluster.ids)
DimPlot(merged_seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#Combine clusters 10 and 12
merged_seurat_filtered$combined_cluster <- as.character(merged_seurat_filtered$seurat_clusters)
merged_seurat_filtered$combined_cluster[merged_seurat_filtered$seurat_clusters %in% c("10", "12")] <- "Cluster_10_12"
merged_seurat_filtered$combined_cluster <- factor(merged_seurat_filtered$combined_cluster) 
merged_seurat_filtered@meta.data$original_cluster <- merged_seurat_filtered$seurat_clusters
cluster_10_12 <- subset(merged_seurat_filtered, subset = combined_cluster == "Cluster_10_12")

StemCells_YvsO_markers <- FindMarkers(cluster_10_12, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")

StemCells_YvsO_markers <- StemCells_YvsO_markers %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

write.csv(StemCells_YvsO_markers, "filtered_StemCells_YvsO_markers.csv")
dim(StemCells_YvsO_markers) #7871 markers, 2230 after filtered

StemCells_YvsO_markers$gene <- rownames(StemCells_YvsO_markers)
top20 <- StemCells_YvsO_markers %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(StemCells_YvsO_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 10 + 12")

FeaturePlot(cluster_10_12, features = c("Hspa1a"))
VlnPlot(cluster_10_12, features = c("Ano6"), group.by = "original_cluster", split.by = "orig.ident")
table(Idents(merged_seurat_filtered))
table(merged_seurat_filtered$seurat_clusters)


#Combine clusters 0, 10, 12
merged_seurat_filtered$combined_cluster <- as.character(merged_seurat_filtered$seurat_clusters)
merged_seurat_filtered$combined_cluster[merged_seurat_filtered$seurat_clusters %in% c("0", "10", "12")] <- "Cluster_0_10_12"
merged_seurat_filtered$combined_cluster <- factor(merged_seurat_filtered$combined_cluster)
merged_seurat_filtered@meta.data$original_cluster <- merged_seurat_filtered$seurat_clusters
cluster_0_10_12 <- subset(merged_seurat_filtered, subset = combined_cluster == "Cluster_0_10_12")

StemCells_YvsO_markers_2 <- FindMarkers(cluster_0_10_12, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")

StemCells_YvsO_markers_2 <- StemCells_YvsO_markers_2 %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

StemCells_YvsO_markers_2$gene <- rownames(StemCells_YvsO_markers_2)
write.csv(StemCells_YvsO_markers_2, "filtered _StemCells_YvsO_markers_+0.csv")
dim(StemCells_YvsO_markers_2)#8913 markers, 4877 after filtered

StemCells_YvsO_markers_2$gene <- rownames(StemCells_YvsO_markers_2)
top20_2 <- StemCells_YvsO_markers_2 %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(StemCells_YvsO_markers_2, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20_2, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 0 + 10 + 12")

FeaturePlot(cluster_0_10_12, features = c())
VlnPlot(cluster_0_10_12, features = c("Tanc1", "Satb2", "Tet2", "Tgbfr1","Lpar1","Yap1",
                                      "Ptger4"), group.by = "original_cluster", split.by = "orig.ident")+
  theme(legend.position = "right")

###################### individual clusters

# Cluster 0
Cluster_0 <- subset(merged_seurat_filtered, subset = original_cluster == 0)

YvsO_markers_0 <- FindMarkers(Cluster_0, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")

YvsO_markers_0 <- YvsO_markers_0 %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

YvsO_markers_0$gene <- rownames(YvsO_markers_0)
write.csv(YvsO_markers_0, "filtered_cluster0_YvsO_markers.csv")
dim(YvsO_markers_0)#3939 after filtered

top20_0 <- YvsO_markers_0 %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(YvsO_markers_0, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20_0, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 0")

FeaturePlot(Cluster_0, features = c())
VlnPlot(Cluster_0, features = c("Tanc1", "Satb2", "Galnt12", "Tgfbr1","Lpar1","Yap1",
                                      "Ptger4"), split.by = "orig.ident")

#Cluster 10
Cluster_10 <- subset(merged_seurat_filtered, subset = original_cluster == 10)

YvsO_markers_10 <- FindMarkers(Cluster_10, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")

YvsO_markers_10 <- YvsO_markers_10 %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

YvsO_markers_10$gene <- rownames(YvsO_markers_10)
write.csv(YvsO_markers_10, "filtered_cluster10_YvsO_markers.csv")
dim(YvsO_markers_10)#2521 after filtered

top20_10 <- YvsO_markers_10 %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(YvsO_markers_10, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20_10, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 10")

FeaturePlot(Cluster_10, features = c())
VlnPlot(Cluster_10, features = c("Satb2", "Cdkn2b", "Galnt12", "Ptger4"), split.by = "orig.ident")

#Cluster 12
Cluster_12 <- subset(merged_seurat_filtered, subset = original_cluster == 12)

YvsO_markers_12 <- FindMarkers(Cluster_12, ident.1 = "Young", ident.2 = "Old", group.by = "orig.ident")

YvsO_markers_12 <- YvsO_markers_12 %>%
  mutate(abs_avg_log2FC = abs(avg_log2FC)) %>% 
  filter(!(pct.1 < 0.1 & pct.2 < 0.1),  
         p_val <= 0.05) %>%  
  arrange(desc(abs_avg_log2FC))

YvsO_markers_12$gene <- rownames(YvsO_markers_12)
write.csv(YvsO_markers_12, "filtered_cluster12_YvsO_markers.csv")
dim(YvsO_markers_12)#1158 after filtered

top20_12 <- YvsO_markers_12 %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(YvsO_markers_12, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  geom_text_repel(data = top20_12, aes(label = gene), size = 3) +
  theme_minimal() +
  ggtitle("Young vs Old DEGs in Cluster 12")

FeaturePlot(Cluster_12, features = c())
VlnPlot(Cluster_12, features = c("Zfp992", "Armcx3", "Muc16", "Cox6b2","Slc37a2","Parp3",
                                 "AY761184", "Ptpn6", "Clic6", "Defa30", "Bex1", "Trim35"), split.by = "orig.ident")




#####################

######enrichment analysis 
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr) 

###cluster 10,12
# Step 1: Select Top 200 Genes
top_genes_1012 <- StemCells_YvsO_markers %>%
  head(200)

gene_list_1012 <- top_genes_1012$gene

###cluster 0,10,12
# Step 1: Select Top 200 Genes
top_genes_01012 <- StemCells_YvsO_markers_2 %>%
  head(200)

gene_list_01012 <- top_genes_01012$gene

###cluster 0
# Step 1: Select Top 200 Genes
top_genes_0 <- YvsO_markers_0 %>%
  head(200)

gene_list_0 <- top_genes_0$gene

###cluster 10
# Step 1: Select Top 200 Genes
top_genes_10 <- YvsO_markers_10 %>%
  head(200)

gene_list_10 <- top_genes_10$gene

###cluster 12
# Step 1: Select Top 200 Genes
top_genes_12 <- YvsO_markers_12 %>%
  head(200)

gene_list_12 <- top_genes_12$gene

### continue
# Step 2: Map Gene Symbols to Entrez IDs
gene_entrez <- bitr(gene_list_01012, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
gene_entrez <- gene_entrez[!duplicated(gene_entrez$ENTREZID), ]

# Step 3: Perform GO Enrichment Analysis
ego <- enrichGO(gene         = gene_entrez$ENTREZID,
                OrgDb        = org.Mm.eg.db,
                keyType      = "ENTREZID",
                ont          = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable     = TRUE)

ego_df <- as.data.frame(ego)
write.csv(ego_df, "Cluster01012_eGO")

# Barplot for GO BP
barplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment")

# Step 4: Perform KEGG Enrichment Analysis
ekegg <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                    organism     = 'mmu',
                    keyType      = 'kegg',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 0.5)

# Visualize KEGG Results
barplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")

ekegg_df <- as.data.frame(ekegg)





###################### individual samples
# Load each file individually into a list
young_data_list <- lapply(young_files, function(file) {
  read.csv(file, row.names = 1)
})
names(young_data_list) <- young_files

old_data_list <- lapply(old_files, function(file) {
  read.csv(file, row.names = 1)
})
names(old_data_list) <- old_files

#Create Seurat Objects for Each File
young_seurat_list <- lapply(young_data_list, function(data) {
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  return(seurat_obj)
})
names(young_seurat_list) <- gsub("/Users/chiarasantoso/Desktop/GSE168448_RAW/", "", names(young_data_list))

old_seurat_list <- lapply(old_data_list, function(data) {
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  return(seurat_obj)
})
names(old_seurat_list) <- gsub("/Users/chiarasantoso/Desktop/GSE168448_RAW/", "", names(old_seurat_list))

# Combine all Young and Old Seurat objects into one list
all_seurat_list <- c(young_seurat_list, old_seurat_list)

# Create a metadata column for sample identity
for (name in names(all_seurat_list)) {
  all_seurat_list[[name]]$sample <- name
}
for (name in names(young_seurat_list)) {
  young_seurat_list[[name]]$sample <- name
}
for (name in names(old_seurat_list)) {
  old_seurat_list[[name]]$sample <- name
}

# Merge all objects into one for visualization
merged_seurat_individual_samples <- merge(all_seurat_list[[1]], y = all_seurat_list[-1], add.cell.ids = names(all_seurat_list))
young_seurat_individual_samples <- merge(young_seurat_list[[1]], y = young_seurat_list[-1], add.cell.ids = names(young_seurat_list))
old_seurat_individual_samples <- merge(old_seurat_list[[1]], y = old_seurat_list[-1], add.cell.ids = names(old_seurat_list))

# QC
VlnPlot(young_seurat_individual_samples, features = "percent.mt", group.by = "sample") +
  ggtitle("Mitochondrial Gene Percentage by Sample (Young)")
VlnPlot(old_seurat_individual_samples, features = "percent.mt", group.by = "sample") +
  ggtitle("Mitochondrial Gene Percentage by Sample (Old)")

VlnPlot(young_seurat_individual_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3)
VlnPlot(old_seurat_individual_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", ncol = 3)


##########################

