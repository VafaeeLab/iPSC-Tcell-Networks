#in this script we are exploring the space of k-vals from 10 to 70, keeping
#ms constant at 10
#unlike our previous optimisation, we are going to store mean AND median
#cohesion. We are also going to store the connectivity of the TOM and the EGADs
#scores

# Load necessary libraries
library(Seurat)
library(hdWGCNA)
library(dplyr)
library(ggplot2)
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



setwd("~/ipsc_project/R_scripts/iPSC_analysis")

seurat_obj <- readRDS("reference_analysis.RDS")


#### -- producing gene list -- ####

conv_T_obj <- create_seurat("conv_T", count_files[[2]], meta_files[[2]])

#firstly get HVGs from conventional T cells
conv_T_obj <- NormalizeData(conv_T_obj)
conv_T_obj <- FindVariableFeatures(conv_T_obj)
gene_list <- VariableFeatures(conv_T_obj)

#now load up DEGs between conv. and innate T cells
DEG_list <- read.csv(file = "reference_data_analysis/conv_vs_innate_DEGs.csv", row.names = 1)

# IMPORTANT NOTE: upregulated genes have a positive FC in the DEG list
# Filter for significant and upregulated genes
significant_upregulated_genes <- DEG_list[DEG_list$p_val_adj < 0.05 & DEG_list$avg_log2FC > 0, ]
significant_upregulated_genes <- rownames(significant_upregulated_genes)

gene_list <- union(gene_list, significant_upregulated_genes)

covariates <- c("donor", "method")

#remove biologically irrevelant genes (rRNA and mt)
# Identify mitochondrial genes
# For human data (e.g., "MT-ND1")
is_mt <- grepl("^MT-", gene_list, ignore.case = FALSE)
# Identify ribosomal genes (genes starting with RPS or RPL)
is_rp <- grepl("^RP[SL]", gene_list, ignore.case = FALSE)

cc_genes <- union(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

is_cc <- gene_list %in% cc_genes

# Combine the logical vectors to identify any gene that is mitochondrial OR ribosomal
genes_to_exclude <- is_mt | is_rp | is_cc

# Create the final list of genes to keep for WGCNA
gene_list <- gene_list[!genes_to_exclude]

#create the gene list to use on the cd8_obj (if not using merged)
filtered_gene_list <- intersect(gene_list, rownames(seurat_obj))


#### EGAD setup ####
library(EGAD)
library(biomaRt)
data(biogrid)
data(GO.human)

GO <- GO.human

# Make your gene list and the network 
genelist <- make_genelist(biogrid)

# Store your annotation matrix
goterms <- unique(GO.human[,3])
annotations <- make_annotations(GO.human[,c(2,3)],genelist,goterms)

# Select the Ensembl genes mart
ensembl_mart <- useEnsembl(biomart = "genes")
# Select the human dataset
# To see all available datasets: datasets <- listDatasets(ensembl_mart)
human_dataset <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl_mart)


# Get GO terms for the genes in your network
go_annotations <- getBM(
  attributes = c('hgnc_symbol', 'go_id', 'namespace_1003'),
  filters = 'hgnc_symbol',
  values = filtered_gene_list,
  mart = human_dataset
)

# It's a good practice to remove blank GO IDs
go_annotations <- go_annotations[go_annotations$go_id != '', ]

# Prepare the data for make_annotations
# The first column should be gene names, the second should be GO IDs
annotations_data <- go_annotations[, c("hgnc_symbol", "go_id")]

# Get the unique gene names and GO terms
gogenes <- unique(annotations_data$hgnc_symbol)
goterms <- unique(annotations_data$go_id)

# Create the binary annotation matrix
annotations <- make_annotations(annotations_data, gogenes, goterms)


### now that we have our gene list, we can perform network analysis with hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  wgcna_name = "conventional_optimisation",
  features = filtered_gene_list #otherwise it will use var features
)

## Extract ids of reference cells
#we will pass this to cells.use
ref_ids <- WhichCells(seurat_obj, expression = orig.ident == "conv_T") 

## Params and results setup
# Create a parameter grid data frame to define the combinations to test
# We will tune 'k' (number of nearest neighbors) and 'max_shared'
param_grid <- expand.grid(
  k = c(10, 20, 25, 30, 35, 40, 50, 60, 70), #tutorial advises between 20 and 75, but higher for large datasets
  max_shared = c(10) #just running with ms = 10
)

# Initialize an empty data frame to store the results
results_df <- data.frame(
  k = integer(),
  max_shared = integer(),
  num_metacells = integer(),
  soft_power = integer(),
  EGADs_auroc = numeric(),
  mean_connectivity = numeric(),
  median_connectivity = numeric(),
  num_modules = integer(),
  avg_module_cohesion = numeric(),
  median_module_cohesion = numeric()
)

cat("Starting hdWGCNA parameter optimization loop...\n")

# Loop through each row of the parameter grid
for (i in 1:nrow(param_grid)) {
  
  # Get the current set of parameters
  k_val <- param_grid$k[i]
  max_shared_val <- param_grid$max_shared[i]
  
  cat(paste0("\nRunning iteration ", i, "/", nrow(param_grid), ": k = ", k_val, ", max_shared = ", max_shared_val, "\n"))
  
  # Use tryCatch to handle potential errors
  new_row <- tryCatch({
    
    # --- 2.1: Reset Object ---
    seurat_obj_loop <- seurat_obj
    
    # --- 2.2: Run hdWGCNA Pipeline ---
    
    # a) Create metacells
    seurat_obj_loop <- MetacellsByGroups(
      seurat_obj = seurat_obj_loop,
      group.by = c("celltype_annotation", covariates),
      k = k_val,
      max_shared = max_shared_val,
      ident.group = 'celltype_annotation',
      cells.use = ref_ids, #this subsets to our reference cells, as we are only inferring networks from them
      reduction = "pca" #using pca is suitable when we are also selecting per sample
    )
    
    # b) Normalize metacells
    metacell_obj <- GetMetacellObject(seurat_obj_loop)
    metacell_obj <- NormalizeData(metacell_obj)
    
    #Getting number of metacells
    num_metacells_val <- length(Cells(metacell_obj))
    
  
    
    
    # c) Set up the expression matrix for our group of interest
    #this is how we backdoor in our normalised metacell data (we get errors,
    #presumedly because we used the cells.use param)
    metacell_obj_subset <- subset(metacell_obj, features = filtered_gene_list) #subset to genes of interest
    exp <- SeuratObject::LayerData(metacell_obj_subset, assay="RNA", layer="data") #data layer gets normalised metacell data
    exp <- t(exp) #transpose to set cols as genes and rows as cells (metacells)
    seurat_obj_loop <- SetDatExpr(
      seurat_obj_loop,
      group_name = "conv_reference",
      mat = exp
    )
    
    # d) Select soft-power threshold
    seurat_obj_loop <- TestSoftPowers(seurat_obj_loop)
    
    # e) Construct the co-expression network
    tom_name <- paste0("TOM_k", k_val, "_ms", max_shared_val)
    seurat_obj_loop <- ConstructNetwork(
      seurat_obj_loop,
      tom_name = tom_name,
      overwrite_tom = TRUE
    )
    
    #extract the min soft power required to achieve 0.8 scale free topology
    soft_power_val = seurat_obj_loop@misc$conventional_optimisation$wgcna_params$power
    
    #examining the connectivity
    connectivity <- WGCNA::softConnectivity(as.matrix(exp), power = soft_power_val)
    
    mean_k <- mean(connectivity)
    median_k <- median(connectivity)
    
    
    # 3.1: EGADs AUROC
    tom <- GetTOM(seurat_obj_loop)
    genes_in_tom <- rownames(tom)
    common_genes <- intersect(genes_in_tom, rownames(annotations))
    annotations_subset <- annotations[common_genes, ]
    filtered_tom <- tom[common_genes, common_genes]
    
    egads_results <- run_GBA(tom, annotations_subset, min = 10, max = Inf)
    egads_auroc_val <- egads_results[[3]]
    
    
    # f) Compute module eigengenes, connectivity, and enrichment
    seurat_obj_loop <- ModuleEigengenes(seurat_obj_loop,  
                                        group.by.vars = covariates, 
                                        verbose = FALSE) #this will give us harmonised and non-harmonised MEs)
    seurat_obj_loop <- ModuleConnectivity(seurat_obj_loop,
                                          group.by = "orig.ident", 
                                          group_name = "conv_T",
                                          harmonized = TRUE)
    
    modules <- GetModules(seurat_obj_loop)
    
    num_modules_val <- (length(unique(modules$color)) - 1) #subtract 1 to remove grey module
    
    # 3.2: Module Cohesion Score
    kME_table <- modules %>%
      dplyr::filter(module != 'grey') %>%
      rowwise() %>%
      mutate(kME_own_module = get(paste0("kME_", module))) %>%
      ungroup() %>%
      group_by(module) %>%
      summarise(avg_kME = mean(kME_own_module), .groups = 'drop')
    
    avg_module_cohesion_val <- mean(kME_table$avg_kME)
    median_module_cohesion_val <- median(kME_table$avg_kME)
    
    
    # deleting TOM
    TOM_path <- paste0("TOM/", tom_name, "_TOM.rda")
    file.remove(TOM_path)
    
    # --- 4: Storing Results ---
    data.frame(
      k = k_val,
      max_shared = max_shared_val,
      num_metacells = num_metacells_val,
      soft_power = as.integer(soft_power_val),
      EGADs_auroc = egads_auroc_val,
      mean_connectivity = mean_k,
      median_connectivity = median_k,
      num_modules = num_modules_val,
      avg_module_cohesion = avg_module_cohesion_val,
      median_module_cohesion = median_module_cohesion_val
    )
    
  }, error = function(e) {
    # This function is executed if an error occurs in the code above
    cat(paste0("Error in iteration ", i, ": ", e$message, "\n"))
    
    # Return a data frame with NA values for the metrics
    data.frame(
      k = k_val,
      max_shared = max_shared_val,
      num_metacells = NA,
      soft_power = NA,
      EGADs_auroc = NA,
      mean_connectivity = NA,
      median_connectivity = NA,
      num_modules = NA,
      avg_module_cohesion = NA,
      median_module_cohesion = NA
    )
  })
  
  # Add the new row (either with results or NA) to the results data frame
  results_df <- rbind(results_df, new_row)
  
  # Optionally, print a summary of the iteration's outcome
  if (!is.na(new_row$avg_module_cohesion)) {
    cat(paste0("Iteration ", i, " complete. Cohesion: ", round(new_row$avg_module_cohesion, 3)))
  } else {
    cat(paste0("Iteration ", i, " failed. Storing NA values.\n"))
  }
}
print(results_df)



write.csv(results_df, file = "reference_data_analysis/conv_T_optimisation/conv_T_opt_2_cc_results.csv")

results_df <- read.csv("reference_data_analysis/conv_T_optimisation/conv_T_opt_2_cc_results.csv")

print(results_df)
results_df <- results_df[, -1]
