library(SeuratData)
library(SeuratDisk)
library(png)
library(Seurat)
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
Convert("/Users/guojuanru/Desktop/Petti/data/spot_allen_cortex.h5ad", dest = "h5seurat", 
        assay = "Spatial", overwrite = TRUE)
spotdata <- LoadH5Seurat("/Users/guojuanru/Desktop/Petti/data/spot_allen_cortex.h5seurat")
image = readPNG("/Users/guojuanru/Desktop/Petti/data/tissue_lowres_image.png")
spotdata@images = brain@images
names(spotdata@images) = "Anndata"
spotdata@images[["Anndata"]]@image = image
spotdata@images[["Anndata"]]@key = 'Anndata_'
spotdata@images[["Anndata"]]@scale.factors[["lowres"]] = 1/23
spotdata@images[["Anndata"]]@scale.factors[["hires"]] = spotdata@images[["Anndata"]]@scale.factors[["lowres"]]*(10/3)

position = read.csv("/Users/guojuanru/Desktop/Petti/data/position.csv",header = FALSE)
colnames(position) = c("row","col","imagerow","imagecol")
spotdata@images[["Anndata"]]@coordinates = position


# here we create the seurate object called spotdata
# Let's take a look of this data
spotdata <- SCTransform(spotdata, assay = "Spatial", verbose = FALSE)
spotdata <- RunPCA(spotdata, assay = "SCT", verbose = FALSE, npcs = 30 )
spotdata <- FindNeighbors(spotdata, reduction = "pca", dims = 1:20)
spotdata <- FindClusters(spotdata, verbose = FALSE)
spotdata <- RunUMAP(spotdata, reduction = "pca", dims = 1:20)



DimPlot(spotdata , reduction = "umap")

SpatialDimPlot(spotdata,group.by = "seurat_clusters")

SpatialPlot(spotdata, features = c("Gm9925"))

truth= read.csv("/Users/guojuanru/Desktop/Petti/data/truth.csv",row.names=1)
spotlight = spotlight_ls[[2]][,colnames(truth)]


# here are some methods to calculate the errors
calculate_error<- function(answer = answer, key = key){
  library(fBasics)
  mse = mean(rowMeans((answer-key)^2))
  correlation = diag(cor(answer,key))
  return( list("mse" = mse, "cor" = correlation))
}

answer = spotlight
key = truth/rowSums(truth)
calculate_error(spotlight,truth/rowSums(truth))

