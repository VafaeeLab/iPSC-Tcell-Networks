#in this script, we are going to analyse modules present in our conventional
#T cells, comparing their expression across other in vivo cell types to
#determine which modules are specific to T cell phenotype, immune cell phenotype
#etc. and then project them onto our iPSC query datasets

#We will run with HVGs from the conv. T cells, as well as 
#DEGs between conventional and innate like T cells, and we are also
#removing cell cycle genes. Unlike previous runs were are only going to remove
#cell cycle genes using cc.genes provided by Seurat

# Load necessary libraries
library(Seurat)
library(hdWGCNA)
library(dplyr)
library(ggplot2)
library(harmony)
library(enrichR)

setwd("~/ipsc_project/R_scripts/iPSC_analysis")

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


#makes enrichment plots
make_enrichment_dot_plots <- function(seurat_obj, dbs) {
  
  plots <- list()
  for (db in dbs) {
    # enrichr dotplot
    p <- EnrichrDotPlot(
      seurat_obj,
      mods = "all", # use all modules (default)
      database = db, # this must match one of the dbs used previously
      n_terms=2, # number of terms per module
      term_size=8, # font size for the terms
      p_adj = FALSE, # show the p-val or adjusted p-val?
      p_cutoff = 0.1
    )  + scale_color_stepsn(colors=rev(viridis::magma(256)))
    plots[[db]] <- p
  }
  return(plots)
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

#### loop through the file path list, creating Seurat objs for each cell type ####

for (dataset in names(file_paths)) {
  print(paste("Processing:", dataset))
  # Correctly access the file paths using the 'dataset' variable
  current_paths <- file_paths[[dataset]]
  
  obj <- create_seurat(dataset_name = dataset, 
                       countdata_file = current_paths["count"], 
                       metadata_file = current_paths["meta"])
  
  obj_name <- paste0(dataset, "_obj")
  assign(obj_name, obj)
}

#Create a named list of our individual objs
obj_list <- list(conv_T = conv_T_obj, 
                 t1_innate_T = t1_innate_T_obj, 
                 monocyte = monocyte_obj,
                 B = B_obj,
                 fibroblast = fibroblast_obj)

#merge the individual datasets into one seurat obj
seurat_obj <- merge(obj_list[['conv_T']], 
                    c(obj_list[['t1_innate_T']], 
                      obj_list[['monocyte']], 
                      obj_list[['B']], 
                      obj_list[['fibroblast']]))

seurat_obj <- JoinLayers(seurat_obj)

#new feature: covariates
#covariates will be the vector of metadata fields we will integrate against
covariates <- c("donor", "method")

#process and visualise the dataset
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

ElbowPlot(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

DimPlot(seurat_obj, group.by = c("method", "donor", "celltype_annotation"))



seurat_obj <- RunHarmony(object = seurat_obj, 
                         group.by.vars = covariates, 
                         reduction.use = "pca")

ElbowPlot(seurat_obj, reduction = "harmony")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:10, reduction.name = "harmony.umap")

DimPlot(seurat_obj, group.by = c("method", "donor", "celltype_annotation"), reduction = "harmony.umap")

#saveRDS(seurat_obj, file = "conv_T.RDS")
#seurat_obj <- readRDS(file = "conv_T.RDS")


#### -- producing gene list -- ####


gene_list <- readLines("conv_gene_list.txt")

#create the gene list to use on the cd8_obj (if not using merged)
filtered_gene_list <- intersect(gene_list, rownames(seurat_obj))

### now that we have our gene list, we can perform network analysis with hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  wgcna_name = "conv_T_analysis",
  features = filtered_gene_list #otherwise it will use var features
)

## Extract ids of reference cells
#we will pass this to cells.use
ref_ids <- WhichCells(seurat_obj, expression = orig.ident == "conv_T") 

#will now run with k = 30, ms = 10
k_val <- 25
max_shared_val <- 10
covariates <- c("donor", "method")
# a) Create metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("orig.ident","celltype_annotation", covariates),
  k = k_val,
  max_shared = max_shared_val,
  ident.group = 'celltype_annotation',
  cells.use = ref_ids, #this subsets to our reference cells, as we are only inferring networks from them
  reduction = "pca" #using pca is suitable when we are also selecting per sample
)

# b) Normalize metacells
metacell_obj <- GetMetacellObject(seurat_obj)
metacell_obj <- NormalizeData(metacell_obj)

### Optional section: processing and examining the metacells
#optional: process and visualise the metacells, investigate their properties

metacell_obj <- FindVariableFeatures(metacell_obj)
metacell_obj <- ScaleData(metacell_obj)
metacell_obj <- RunPCA(metacell_obj)
metacell_obj <- RunUMAP(metacell_obj, dims = 1:15)
metacell_obj <- RunHarmony(metacell_obj, group.by.vars = covariates)
metacell_obj <- RunUMAP(metacell_obj, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:10)
DimPlot(metacell_obj, group.by = c("celltype_annotation", covariates)) + NoLegend()
DimPlot(metacell_obj, reduction = "harmony.umap", group.by = c("celltype_annotation", covariates))


### End of optional section

#this is how we backdoor in our normalised metacell data (we get errors,
#presumedly because we used the cells.use param)
metacell_obj_subset <- subset(metacell_obj, features = filtered_gene_list) #subset to genes of interest
exp <- SeuratObject::LayerData(metacell_obj_subset, assay="RNA", layer="data") #data layer gets normalised metacell data
exp <- t(exp) #transpose to set cols as genes and rows as cells (metacells)

# c) Set up the expression matrix for our group of interest
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "conv_reference",
  mat = exp
)

# d) Select soft-power threshold
seurat_obj <- TestSoftPowers(seurat_obj)
# run soft power test


# e) Construct the co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = paste0("TOM_k", k_val, "_ms", max_shared_val),
  overwrite_tom = TRUE,
)
PlotDendrogram(seurat_obj)

# f) Compute module eigengenes, connectivity, and enrichment
seurat_obj <- ModuleEigengenes(seurat_obj, 
                               group.by.vars = covariates) #this will give us harmonised and non-harmonised MEs


seurat_obj <- ModuleConnectivity(seurat_obj, 
                                 group.by = "celltype_annotation", 
                                 group_name = c("CD8+T", "CD4+T"),
                                 harmonized = TRUE)

saveRDS(seurat_obj, file = "conv_T.RDS")
seurat_obj <- readRDS("conv_T.RDS")

# Optional: Module Cohesion Score
modules <- GetModules(seurat_obj)
kME_table <- modules %>%
  dplyr::filter(module != 'grey') %>%
  rowwise() %>%
  mutate(kME_own_module = get(paste0("kME_", module))) %>%
  ungroup() %>%
  group_by(module) %>%
  summarise(avg_kME = mean(kME_own_module), .groups = 'drop')

avg_module_cohesion_val <- mean(kME_table$avg_kME)
median_module_cohesion_val <- median(kME_table$avg_kME)

#plotting correlations of MEs to examine any likely shared functions
eig <- GetMEs(seurat_obj)
eig <- eig[,colnames(eig) != "grey"] #remove grey module
cor_mat <- cor(eig, method = "pearson") 
pheatmap::pheatmap(cor_mat, clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   main = "Module–module eigengene correlation")


#### ---- EXAMINING EXPRESSION ACROSS CELL TYPES ---- ####

cell_types <- seurat_obj$celltype_annotation
hMEs <- GetMEs(seurat_obj, harmonized = TRUE)
hMEs <- hMEs[,colnames(hMEs) != "grey"] #remove grey module


hME_by_type <- cbind(cell_type = celltypes, as.data.frame(hMEs)) %>%
  group_by(cell_type) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

library(pheatmap)

palette_length <- 100
my_color <- colorRampPalette(c("darkblue", "white", "red"))(palette_length)
my_breaks <- seq(min(hME_by_type), max(hME_by_type), length.out = palette_length)

pheatmap(hME_by_type,
         color = my_color,
         breaks = my_breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         main = "Mean module eigengene per cell type")


DimPlot(seurat_obj, reduction = "harmony.umap", group.by = "orig.ident")
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  reduction = "harmony.umap",
  order=TRUE # order so the points with highest hMEs are on top
)
# stitch together with patchwork
patchwork::wrap_plots(plot_list, ncol=3)

### Some dME analysis as well ###


## immune cells vs non-immune cells
group1 <- seurat_obj@meta.data %>% subset(orig.ident == c('conv_T', "t1_innate_T", 'B', "monocyte")) %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 'fibroblast') %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)

non_immune_mods <- c("purple", "black", "brown", "red")


##next identify modules unique to T vs other immune cells
group1 <- seurat_obj@meta.data %>% subset(orig.ident == c('conv_T', "t1_innate_T")) %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == c('B', "monocyte")) %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)

non_T_mods <- c(non_immune_mods, "yellow", "green", "red")
T_mods <- c("pink", "turquoise", "magenta", "blue")

##firstly identify modules unique to conv T cells vs innate T
group1 <- seurat_obj@meta.data %>% subset(orig.ident == 'conv_T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 't1_innate_T') %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)

#subset to modules unique to T-cells
T_DMEs <- DMEs[T_mods,]

PlotDMEsLollipop(
  seurat_obj, 
  T_DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  T_DMEs,
  wgcna_name = 'conv_T',
  
)

##firstly identify modules unique to CD4 vs CD8 T cells
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD8+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD4+T') %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)
#subset to modules unique to T-cells
T_DMEs <- DMEs[T_mods,]

PlotDMEsLollipop(
  seurat_obj, 
  T_DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  T_DMEs,
  wgcna_name = 'conv_T',
  
)


##innate vs CD8 T
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD8+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 't1_innate_T') %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)

#subset to modules unique to T-cells
T_DMEs <- DMEs[T_mods,]

PlotDMEsLollipop(
  seurat_obj, 
  T_DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  T_DMEs,
  wgcna_name = 'conv_T',
  
)


##innate vs CD4
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD4+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 't1_innate_T') %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)
PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-8,5),
  DMEs,
  wgcna_name = 'conv_T',
  
)

#subset to modules unique to T-cells
T_DMEs <- DMEs[T_mods,]

PlotDMEsLollipop(
  seurat_obj, 
  T_DMEs, 
  wgcna_name = "conv_T",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  T_DMEs,
  wgcna_name = 'conv_T',
  
)

#### ---- ENRICHMENT ANALYSIS ---- ####

# g) Run functional enrichment analysis
dbs <- c('GO_Biological_Process_2023',
         'GO_Cellular_Component_2023',
         'GO_Molecular_Function_2023', 
         'CellMarker_2024',
         'KEGG_2021_Human',
         'MSigDB_Hallmark_2020')

seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs = dbs,
  max_genes = 100
)

dot_plots <- make_enrichment_dot_plots(seurat_obj, dbs)

dot_plots[['GO_Biological_Process_2023']] + 
  dot_plots[['GO_Cellular_Component_2023']] + 
  dot_plots[['GO_Molecular_Function_2023']]

dot_plots[['CellMarker_2024']] + 
  dot_plots[['KEGG_2021_Human']] + 
  dot_plots[['MSigDB_Hallmark_2020']]

modules <- GetModules(seurat_obj) #examining modules manually
magenta <- subset(modules, module == "magenta")
magenta <- magenta[,c("gene_name","kME_magenta")]

#we will now perform GSEA with the turquoise module to investigate further
# load the GO Biological Pathways file (downloaded from EnrichR website)
pathways <- fgsea::gmtPathways('../GO_Biological_Process_2023.txt')
# optionally, remove the GO term ID from the pathway names to make the downstream plots look cleaner
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")
head(pathways)

# get the modules table and remove grey genes
modules <- GetModules(seurat_obj)

# rank all genes in this kME
cur_mod <- "pink" 
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]
# run fgsea to compute enrichments
gsea_df <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 10,
  maxSize = 500
)
head(gsea_df)

top_pathways <- gsea_df %>% 
  subset(pval < 0.05) %>% 
  slice_max(order_by=NES, n=25) %>% 
  .$pathway

fgsea::plotGseaTable(
  pathways[top_pathways], 
  ranks, 
  gsea_df, 
  gseaParam=0.5,
  colwidths = c(10, 4, 1, 1, 1)
)

# get the modules table and remove grey genes
modules <- GetModules(seurat_obj) %>% subset(module == cur_mod)

# rank all genes in this kME
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]
# run fgsea to compute enrichments
gsea_df <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 10,
  maxSize = 500
)
head(gsea_df)

top_pathways <- gsea_df %>% 
  subset(pval < 0.05) %>% 
  slice_max(order_by=NES, n=25) %>% 
  .$pathway

fgsea::plotGseaTable(
  pathways[top_pathways], 
  ranks, 
  gsea_df, 
  gseaParam=0.5,
  colwidths = c(10, 4, 1, 1, 1)
)

#we will now perform GSEA with CELLULAR COMPONENTS
pathways <- fgsea::gmtPathways('GO_Cellular_Component_2025.txt')
# optionally, remove the GO term ID from the pathway names to make the downstream plots look cleaner
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")
head(pathways)

# get the modules table and remove grey genes
modules <- GetModules(seurat_obj)  %>% subset(module != 'grey')

# rank all genes in this kME
cur_mod <- "black" 
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]
# run fgsea to compute enrichments
gsea_df <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 5, #need to have a smaller minSize to capture TCR
  maxSize = 500
)
head(gsea_df)

top_pathways <- gsea_df %>% 
  subset(pval < 0.05) %>% 
  slice_max(order_by=NES, n=25) %>% 
  .$pathway

fgsea::plotGseaTable(
  pathways[top_pathways], 
  ranks, 
  gsea_df, 
  gseaParam=0.5,
  colwidths = c(10, 4, 1, 1, 1)
)

##based on all the above analysis we will rename the modules
rename_list <- list(
  "magenta" = "TCR",
  "turquoise" = "T_cell_activation",
  "pink" = "all_T_cytokine_signalling",
  "blue" = "conv_T_cytokine_signalling",
  "yellow" = "immune_specific_housekeeper1",
  "green" = "immune_specific_housekeeper2",
  "black" = "ubiquitous_housekeeper",
  "red" = "non_immune_housekeeper1",
  "brown" = "non_immune_housekeeper2",
  "purple" = "non_immune_housekeeper3"
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = rename_list
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  reduction = "harmony.umap",
  order=TRUE # order so the points with highest hMEs are on top
)
# stitch together with patchwork
patchwork::wrap_plots(plot_list, ncol=5)

#generating a similar plot using ucell scores
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

ggplot2::theme_set(ggplot2::theme_grey())   # restore the default theme
ggplot2::scale_colour_continuous()          # ensures no global continuous colour scale
ggplot2::scale_fill_continuous()

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the ucell scores
  reduction = "harmony.umap",
  ucell = TRUE,
  order=TRUE # order so the points with highest hMEs are on top
)

DimPlot(seurat_obj, reduction = "harmony.umap", group.by = "celltype_annotation")
patchwork::wrap_plots(plot_list, ncol = 5)


saveRDS(seurat_obj, file = "reference_data_analysis/annotated_reference_obj.rds")
seurat_obj <- readRDS(file = "reference_data_analysis/annotated_reference_obj.rds")
