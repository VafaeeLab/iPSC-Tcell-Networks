#in this script we are visualising the GSEA results for our modules of interest

library(Seurat)
library(hdWGCNA)
library(fgsea)

seurat_obj <- readRDS(file = "~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/annotated_reference_obj.rds")

#we will now perform GSEA with the turquoise module to investigate further
# load the GO Biological Pathways file (downloaded from EnrichR website)

pathways <- fgsea::gmtPathways('~/ipsc_project/R_scripts/GO_Biological_Process_2025.txt')

pathways <- fgsea::gmtPathways('~/ipsc_project/R_scripts/GO_Cellular_Component_2025.txt')

pathways <- fgsea::gmtPathways('~/ipsc_project/R_scripts/GO_Molecular_Function_2025.txt')

# optionally, remove the GO term ID from the pathway names to make the downstream plots look cleaner
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")
head(pathways)

# get the modules table and remove grey genes
modules <- GetModules(seurat_obj)
mod_names <- unique(modules$module)

# rank all genes in this kME
cur_mod <- "T_cell_activation" 
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]
# run fgsea to compute enrichments
gsea_df <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 5,
  maxSize = 500
)
head(gsea_df)

top_pathways <- gsea_df %>% 
  subset(pval < 0.05) %>% 
  slice_max(order_by=NES, n=25) %>% 
  .$pathway

ranked_gsea_df <- gsea_df %>% 
  subset(padj < 0.05) %>% 
  slice_max(order_by=NES, n=25)

# Select relevant columns
table_out <- ranked_gsea_df[, .(pathway, NES, pval, padj)]

# Print nicely
print(table_out)

fgsea::plotGseaTable(
  pathways[top_pathways], 
  ranks, 
  gsea_df, 
  gseaParam=0.5,
  colwidths = c(10, 4, 1, 1, 1)
)

# name of the pathway to plot 
selected_pathway <- 'Positive Regulation of Cytokine Production'

plotEnrichment(
  pathways[[selected_pathway]],
  ranks
) + labs(title=selected_pathway)

#### Subsetting to JUST the module of interest  ####
# get the modules table and remove grey genes
cur_mod <- "T_cell_activation" 
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

# name of the pathway to plot 
selected_pathway <- 'Regulation of Cytokine Production'

plotEnrichment(
  pathways[[selected_pathway]],
  ranks
) + labs(title=selected_pathway)
