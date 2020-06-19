library(dplyr)
library(Seurat)
library(patchwork)


#load pbmc dataset
pbmc.data <- Read10X(data.dir ="filtered_gene_bc_matrices/hg19")
dense.size <- object.size(as.matrix(pbmc.data))


#initialize the seurat object with the raw (non-normalized data).

pbmc <-CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells =3,min.features=200)
print(pbmc)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#to print in a r script you have to make an object of the plot and than call print on it
p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)

#filter cells that have unique feature counts over 2,500 or less than 200 and cells that have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#after the unwanted cells have been filtered, the data must be normalized
#by using LogNormalize method, it normalizes the feature expression measures for each cell by the total expression
#then multiply by a scale factor of 10,000
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize" , scale.factor = 10000)


#calculate subset of features that exhibit high cell to cell variation y are highly expressed in some cells, and lowly expressed in others
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#show the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# #add the plots together so we can see both, call it p2 and print it so it displays (Weird hack not sure why R needs this)
# p2  <- plot1 + plot2
# print(p2)


#downstream analysis continue with scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features= all.genes)


#perform pca

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#print plots
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
print(VizDimLoadings(pbmc, dims = 1:2, reduction = "pca"))

print(DimPlot(pbmc, reduction = "pca"))

#this prints the pbmc pca but highlights the selected feature
print(FeaturePlot(object = pbmc, features = "MS4A1"))

#elbow plot
print(ElbowPlot(pbmc))


#cluster cells based on nearest neighbours algo
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#run non-linear dimensional reduction uMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

print(DimPlot(pbmc, reduction = "umap"))

#save object so we dont have to reload / rerurn everything we did above
# saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
print(DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())

