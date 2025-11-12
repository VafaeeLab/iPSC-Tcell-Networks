#in this script, we are loading the ishi 2d dataset, annotating any
#populations of interest for downstream analysis and then setting it up
#for hdWGCNA optimisation

library(Seurat)

#this function takes the name of the dataset, and the file paths for its count and metadata
#then it returns the seurat obj for that data as output
create_seurat <- function(dataset_name, countdata_file, metadata_file) {
  countdata <- read.table(file = countdata_file)
  metadata <- read.table(file = metadata_file,
                         sep = ",",
                         header = TRUE,
                         row.names = 1)
  
  countdata <- t(countdata) #cols must be cells, rows must be genes
  
  ncells <- length(colnames(countdata))
  obj <- CreateSeuratObject(counts = countdata, 
                            meta.data = metadata, 
                            project = dataset_name,
                            min.cells = 0.05*ncells) #filtering for genes expressed in >5% of cells
  return(obj)
}

## Firstly load up the file paths
ishi2d_countdata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/countdata.tsv"
ishi2d_metadata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/metadata.csv"

seurat_obj <- create_seurat(dataset_name = "ishi2d", countdata_file = ishi2d_countdata_f, metadata_file = ishi2d_metadata_f)

##Visualise dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, group.by = "celltypist_high_predictions")

#no groups of cells which we want to annotate separately
#no batch covariates to add
#therefore we will save as an RDS

saveRDS(seurat_obj, file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")
