#in this script we are loading up our processed reference and ishi2d rds files, 
#using the previously determined optimal params to produce metacells and then
#performing module preservation analysis using our modules from the reference

library(Seurat)
library(hdWGCNA)

setwd("~/ipsc_project/R_scripts/iPSC_analysis/iPSC_dataset_processing/ishi_2d")

#our in vivo ref data, containing conv. T derived modules
seurat_ref <- readRDS(file = "~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/annotated_reference_obj.rds")
#our iPSC query data
seurat_query <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

gene_list <- readLines("~/ipsc_project/R_scripts/iPSC_analysis/conv_gene_list.txt")

#filtering gene list to those in our dataset
filtered_gene_list <- intersect(gene_list, rownames(seurat_query))

#### Prepare query obj for module projection ####

## Firstly set up for WGCNA
seurat_query <- SetupForWGCNA(
  seurat_query,
  wgcna_name = "conv_T_mod_projection",
  features = filtered_gene_list #otherwise it will use var features
)

k_val <- 20 #selected using optimisation
max_shared_val <- 10
num_cells <- length(Cells(seurat_query))

##Construct metacells in the query object
seurat_query <- MetacellsByGroups(
  seurat_obj = seurat_query,
  group.by = c("orig.ident"), #no covariates or groups of interest to group by
  k = k_val,
  max_shared = max_shared_val,
  ident.group = 'orig.ident',
  target_metacells = num_cells, #making it so it is not capped at making 1000 metacells
  reduction = "pca" #using pca is suitable when we are also selecting per sample
)

# normalize metacell expression matrix:
seurat_query <- NormalizeMetacells(seurat_query)

#set metacell expression matrix as our hdWGCNA expression mtx
#this means this exp data will be used in module preservation analysis
seurat_query <- SetDatExpr(
  seurat_query,
  group_name = "ishi2d", # the name of the group of interest in the group.by column
  group.by='orig.ident', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

#verifying that the metacell exp mtx has already been set in the ref obj
exp <- GetDatExpr(seurat_ref)

#### Perform module projection and module preservation ####
#project reference modules on to query object
seurat_query <- ProjectModules(
  seurat_obj = seurat_query,
  seurat_ref = seurat_ref,
  wgcna_name = "reference_analysis",
  wgcna_name_proj="conv_T_mod_projection",
  assay="RNA" # assay for query dataset
)

#perform module preservation analysis
seurat_query <- ModulePreservation(
  seurat_query,
  seurat_ref = seurat_ref,
  name="conv_T_preservation",
  verbose=3,
  seed = 12345,
  parallel = FALSE,
  n_permutations=250 # n_permutations=20 used for the tutorial
)


saveRDS(seurat_query, file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

seurat_query <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")
plot_list <- PlotModulePreservation(
  seurat_query,
  name="conv_T_preservation",
  statistics = "summary"
)

patchwork::wrap_plots(plot_list, ncol=2)

plot_list <- PlotModulePreservation(
  seurat_query,
  name="conv_T_preservation",
  statistics = "all",
  plot_labels=TRUE
)
patchwork::wrap_plots(plot_list, ncol=6)

mod_pres <- GetModulePreservation(seurat_query, wgcna_name = "conv_T_mod_projection", mod_name = "conv_T_preservation")
