#this script represents our current method of producing a gene list
#for co-expression analysis of conventional T cells

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



#saving gene_list
write(gene_list, file = "conv_gene_list.txt")