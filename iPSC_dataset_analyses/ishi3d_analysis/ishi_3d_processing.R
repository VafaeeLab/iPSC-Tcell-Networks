#in this script, we are loading the ishi 3d dataset, annotating any
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
ishi3d_countdata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/countdata.tsv"
ishi3d_metadata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/metadata.csv"

seurat_obj <- create_seurat(dataset_name = "ishi3d", countdata_file = ishi3d_countdata_f, metadata_file = ishi3d_metadata_f)

##Visualise dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

DimPlot(seurat_obj, group.by = "celltypist_high_predictions", label = TRUE)
DimPlot(seurat_obj, group.by = "broad_celltype", label = TRUE)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.01)
DimPlot(seurat_obj)

#Making a new "broad_celltype" metadata label titled "broad_celltype_revised"
#with the new clusters (0.01 resolution)
broad_map <- c(
  "0" = "lymphoid",
  "1" = "erythroid",
  "2" = "myeloid"
)

# Add the new metadata column
new_meta <- broad_map[ as.character(seurat_obj$seurat_clusters) ]
names(new_meta) <- colnames(seurat_obj)   # ensure names match cell barcodes
seurat_obj <- AddMetaData(seurat_obj, metadata = new_meta, col.name = "broad_celltype_revised")

DimPlot(seurat_obj, group.by = "broad_celltype_revised")

saveRDS(seurat_obj, file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")
