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

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# add the plots together so we can see both, call it p2 and print it so it displays (Weird hack not sure why R needs this)
p2  <- plot1 + plot2
print(p2)


#downstream analysis continue with scaling the data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features= all.genes)


#perform pca

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#print plots
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
print(VizDimLoadings(pbmc, dims = 1:2, reduction = "pca"))

print(DimPlot(pbmc, reduction = "pca"))