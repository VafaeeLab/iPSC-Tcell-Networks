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
ATO_T6_countdata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/g2g_ipsc/ATO_T06/ATO_T06_countdata.tsv"
ATO_T6_metadata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/g2g_ipsc/ATO_T06/ATO_T06_metadata.csv"

ATO_tnf_countdata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ATO_tnf/countdata.tsv"
ATO_tnf_metadata_f <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ATO_tnf/metadata.csv"

ato_t6_obj <- create_seurat(dataset_name = "ATO_T6", countdata_file = ATO_T6_countdata_f, metadata_file = ATO_T6_metadata_f)



ato_tnf_obj <- create_seurat(dataset_name = "ATO_tnf", countdata_file = ATO_tnf_countdata_f, metadata_file = ATO_tnf_metadata_f)


#adding cell_line covariate
# TNF = KOLF
# ATO_T6 = mix of KOLF and FIAJ

#assigning cell lines based on orig.ident
ato_t6_obj$cell_line <- ifelse(
  ato_t6_obj$batch == "1",
  "FIAJ",            # value if condition is TRUE
  "KOLF"    # value if condition is FALSE
)

ato_tnf_obj$cell_line <- "KOLF"
#also adding a batch value for integration
ato_tnf_obj$batch <- 3


#also fixing the celltypist metadata labels
ato_tnf_obj$celltypist_high_predictions <- ato_tnf_obj$celltypist_high_level_predictions
ato_tnf_obj$celltypist_low_predictions <- ato_tnf_obj$celltypist_low_level_predictions

ato_tnf_obj$celltypist_high_level_predictions <- NULL
ato_tnf_obj$celltypist_low_level_predictions <- NULL



#finally adding a treatment field, in case orig.ident is lost or altered
#in subsequent merging
ato_tnf_obj$treatment <- "TNF"
ato_t6_obj$treatment <- "WT" #WT for wildtype

seurat_obj <- merge(ato_t6_obj, y = ato_tnf_obj)
seurat_obj <- JoinLayers(seurat_obj)

##Visualise dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

DimPlot(seurat_obj, group.by = c("orig.ident", "batch"))

DimPlot(seurat_obj, group.by = c("cell_line"))

seurat_obj <- RunHarmony(seurat_obj, group.by.vars = c("cell_line", "orig.ident"))
#seurat_obj <- RunHarmony(seurat_obj, group.by.vars = c("cell_line"))
ElbowPlot(seurat_obj, reduction = "harmony")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:10)

DimPlot(seurat_obj, group.by = c("orig.ident", "cell_line"), reduction = "harmony.umap")
DimPlot(seurat_obj, group.by = c("batch"), reduction = "harmony.umap")
DimPlot(seurat_obj, group.by = c("celltypist_high_predictions"), reduction = "harmony.umap")


seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.025)
DimPlot(seurat_obj, reduction = "harmony.umap")

maturity_map <-  c(
  "0" = "mature_T", #mature wild type cells
  "1" = "immature_T", #immature wild type cells
  "2" = "mature_T" #all TNF cells appear mature, we will subset to them later with our orig.ident field
)

# Add the new metadata column
new_meta <- maturity_map[ as.character(seurat_obj$seurat_clusters) ]
names(new_meta) <- colnames(seurat_obj)   # ensure names match cell barcodes
seurat_obj <- AddMetaData(seurat_obj, metadata = new_meta, col.name = "maturity")

DimPlot(seurat_obj, group.by = c("treatment", "maturity"), reduction = "harmony.umap")
saveRDS(seurat_obj, file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")
