#in this script, we are loading the ATO and ATO_tnf datasets, merging them
#annotating any populations of interest for downstream analysis and adding batch
#metdata as necessary

library(Seurat)
library(harmony)

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
EZ_T_countdata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/countdata.tsv"
EZ_T_metadata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/metadata.csv"

seurat_obj <- create_seurat(dataset_name = "EZ_T", EZ_T_countdata_f, EZ_T_metadata_f)

##Visualise dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

DimPlot(seurat_obj, group.by = "celltypist_high_predictions")
DimPlot(seurat_obj, group.by = c("replicate","treatment"))
DimPlot(seurat_obj, group.by = c("orig.ident", "batch"))

seurat_obj <- RunHarmony(seurat_obj, group.by.vars = c("orig.ident"))
ElbowPlot(seurat_obj, reduction = "harmony")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15, reduction = "harmony", reduction.name = "harmony.umap")

DimPlot(seurat_obj, group.by = "celltypist_high_predictions", reduction = "harmony.umap")
DimPlot(seurat_obj, group.by = c("replicate","treatment"), reduction = "harmony.umap")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.025)
DimPlot(seurat_obj, reduction = "harmony.umap")
DimPlot(seurat_obj)

maturity_map <- c(
  "0" = "mature_T", #mature T cells
  "1" = "immature_T", #immature wild type cells
  "2" = "mature_activated_T" #mature T
)

# Add the new metadata column
new_meta <- maturity_map[ as.character(seurat_obj$seurat_clusters) ]
names(new_meta) <- colnames(seurat_obj)   # ensure names match cell barcodes
seurat_obj <- AddMetaData(seurat_obj, metadata = new_meta, col.name = "phenotype")

DimPlot(seurat_obj, reduction = "harmony.umap", group.by = "phenotype")

saveRDS(seurat_obj, file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")
