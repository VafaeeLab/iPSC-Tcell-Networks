#in this script, we are going to load up our type 1 innate T and 
#conventional T cell datasets, pseudobulk them and produce DEGs

# Load necessary libraries
library(Seurat)
library(DESeq2)

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

#### Accessing the count and metadata files ####
# 1. Define the base directory
base_dir <- "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/reference_data"

# 2. Find all countdata files recursively and sort them
count_files <- sort(
  list.files(
    path = base_dir,
    pattern = "_countdata\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
)

# 3. Find all metadata files recursively and sort them
meta_files <- sort(
  list.files(
    path = base_dir,
    pattern = "_metadata\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
)

# 4. Combine the two sorted vectors into a list of pairs
#    mapply applies a function (here, `c` to combine) to corresponding elements
#    of the input vectors. SIMPLIFY=FALSE ensures the output is a list.
file_paths <- mapply(c, count = count_files, meta = meta_files, SIMPLIFY = FALSE)

# 5. (Optional) Set the names of the list elements based on the sample name
#    We use gsub to remove the path and file suffixes, leaving just the sample ID.
names(file_paths) <- gsub("_countdata\\.tsv$", "", basename(count_files))

conv_T_obj <- create_seurat("conv_T", count_files[2], meta_files[2])
innate_T_obj <- create_seurat("innate_T", count_files[5], meta_files[5])

#### merging, pseudobulking and producing DEGs ####
merged_seurat <- merge(conv_T_obj, y = innate_T_obj, add.cell.ids = c("conv_T", "innate_T"), project = "merged_project")

# pseudobulk the counts based on donor-condition-celltype
pseudobulk_seurat <- AggregateExpression(merged_seurat, assays = "RNA", return.seurat = T, group.by = c("orig.ident", #either conv_T or innate T
                                                                                                        "celltype_annotation",
                                                                                                        "donor", #batch covariates as well
                                                                                                        "method"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudobulk_seurat))

#set identity class to be conv vs innate
Idents(pseudobulk_seurat) <- pseudobulk_seurat$celltype_annotation

bulk.mono.de <- FindMarkers(object = pseudobulk_seurat, 
                            ident.1 = c("CD4+T", "CD8+T"), 
                            ident.2 = "TYPE-1-INNATE-T",
                            test.use = "DESeq2")

write.csv(bulk.mono.de, file = "iPSC_analysis/reference_data_analysis/conv_vs_innate_DEGs.csv")
