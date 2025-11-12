#in this script we are exploring the space of k-vals from 10 to 70, keeping
#ms constant at 10
#we are omitting egad and focusing on network and module features

# Load necessary libraries
library(Seurat)
library(hdWGCNA)
library(dplyr)
library(harmony)

seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

num_cells <- length(Cells(seurat_obj))

setwd("~/ipsc_project/R_scripts/iPSC_analysis")

#create the gene list to use
gene_list <- readLines("conv_gene_list.txt")

#filtering gene list to those in our dataset
filtered_gene_list <- intersect(gene_list, rownames(seurat_obj))

#change dir to ishi3d processing dir
cat("Changing dir to ishi3d processing dir...\n")
setwd("iPSC_dataset_processing/ishi_3d")

### now that we have our gene list, we can perform network analysis with hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  wgcna_name = "ishi3d_optimisation",
  features = filtered_gene_list #otherwise it will use var features
)

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
      group.by = c("orig.ident", "broad_celltype_revised"), #no covariates, one field of interest to group by (broad_celltype)
      k = k_val,
      max_shared = max_shared_val,
      target_metacells = num_cells, #setting it so maximum metacells is not capped at 1000
      ident.group = 'broad_celltype_revised',
      reduction = "pca" #using pca is suitable when we are also selecting per sample
    )
    
    metacell_obj <- GetMetacellObject(seurat_obj_loop)
    num_metacells_val <- length(Cells(metacell_obj))
    
    # normalize metacell expression matrix:
    seurat_obj_loop <- NormalizeMetacells(seurat_obj_loop)
    
    seurat_obj_loop <- SetDatExpr(
      seurat_obj_loop,
      group_name = "ishi3d", # the name of the group of interest in the group.by column
      group.by='orig.ident', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
      assay = 'RNA', # using RNA assay
      layer = 'data' # using normalized data
    )
    
    # Test different soft powers:
    seurat_obj_loop <- TestSoftPowers(
      seurat_obj_loop
    )
    
    # e) Construct the co-expression network
    tom_name <- paste0("TOM_k", k_val, "_ms", max_shared_val)
    seurat_obj_loop <- ConstructNetwork(
      seurat_obj_loop,
      tom_name = tom_name,
      overwrite_tom = TRUE
    )
    
    #extract the min soft power required to achieve 0.8 scale free topology
    soft_power_val = seurat_obj_loop@misc$ishi3d_optimisation$wgcna_params$power
    
    exp <- GetDatExpr(seurat_obj_loop, wgcna_name = "ishi3d_optimisation")
    
    #examining the connectivity
    connectivity <- WGCNA::softConnectivity(as.matrix(exp), power = soft_power_val)
    
    mean_k <- mean(connectivity)
    median_k <- median(connectivity)
    
    # f) Compute module eigengenes, cohesion
    seurat_obj_loop <- ModuleEigengenes(seurat_obj_loop) #there are no batch covariates, therefore no need for hMEs
    seurat_obj_loop <- ModuleConnectivity(seurat_obj_loop,
                                          group.by = "orig.ident", 
                                          group_name = "ishi3d",
                                          harmonized = FALSE)
    
    modules <- GetModules(seurat_obj_loop)
    
    num_modules_val <- (length(unique(modules$color)) - 1) #subtract 1 to remove grey module
    
    # 3.2: Module Cohesion Scores
    kME_table <- modules %>%
      dplyr::filter(module != 'grey') %>%
      rowwise() %>%
      mutate(kME_own_module = get(paste0("kME_", module))) %>%
      ungroup() %>%
      group_by(module) %>%
      summarise(avg_kME = mean(kME_own_module), .groups = 'drop')
    
    avg_module_cohesion_val <- mean(kME_table$avg_kME)
    median_module_cohesion_val <- median(kME_table$avg_kME)
    
    # deleting TOM to preserve local storage space
    #TOM_path <- paste0("TOM/", tom_name, "_TOM.rda")
    #file.remove(TOM_path)
    
    # --- 4: Storing Results ---
    data.frame(
      k = k_val,
      max_shared = max_shared_val,
      num_metacells = num_metacells_val,
      soft_power = as.integer(soft_power_val),
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

write.csv(results_df, file = "ishi3d_opt_results.csv")

results_df <- read.csv(file = "ishi3d_opt_results.csv")
results_df <- results_df[, -1]
