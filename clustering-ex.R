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

##to print in a r script you have to make an object of the plot and than call print on it
p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)








