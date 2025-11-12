#in this script, we are examining the expression of hMEs across
#our reference cell populations
library(Seurat)
library(hdWGCNA)
library(patchwork)
setwd("~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis")

seurat_obj <- readRDS(file = "annotated_reference_obj.rds")

seurat_obj <- ConstructNetwork(seurat_obj)

PlotDendrogram(seurat_obj)

seurat_obj <- readRDS(file = "annotated_reference_obj.rds")
#getting conv T cells
cd8_cells <- WhichCells(seurat_obj, expression = celltype_annotation == "CD8+T")
cd4_cells <- WhichCells(seurat_obj, expression = celltype_annotation == "CD4+T")

DimPlot(seurat_obj,
        reduction = "harmony.umap",
        group.by = "celltype_annotation",
        label = TRUE,
        label.size = 6,
        label.box = FALSE,
        repel = TRUE, 
        pt.size = 0.5, 
        cells.highlight = list("CD8+T" = cd8_cells, "CD4+T" = cd4_cells),
        cols.highlight = c("yellow","turquoise"),
        sizes.highlight = 0.5) + ggtitle("In vivo reference cells") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 18)
  )

DimPlot(seurat_obj,
        reduction = "harmony.umap",
        group.by = "celltype_annotation",
        label = TRUE,
        label.size = 6,
        label.box = FALSE,
        repel = TRUE,
        #pt.size = 2,
        cells = conv_T_cells) + ggtitle("In vivo reference conventional T cells")

#plotting ref population
DimPlot(seurat_obj,
        reduction = "harmony.umap",
        group.by = "celltype_annotation",
        label = TRUE,
        label.size = 6,
        label.box = FALSE,
        repel = TRUE)



plot_list <- ModuleFeaturePlot(seurat_obj, reduction = "harmony.umap", module_names = c("T_cell_activation", "conv_T_cytokine_signalling", "all_T_cytokine_signalling", "TCR", "immune_specific_housekeeper1", "immune_specific_housekeeper2", "ubiquitous_housekeeper", "non_immune_housekeeper1", "non_immune_housekeeper2", "non_immune_housekeeper3"))

wrap_plots(plot_list, ncol = 2)
wrap_plots(plot_list, nrow = 2)

small_plot_list <- ModuleFeaturePlot(seurat_obj, reduction = "harmony.umap", module_names = c("T_cell_activation", "conv_T_cytokine_signalling", "all_T_cytokine_signalling", "TCR"))

wrap_plots(small_plot_list, nrow=1)

#### ---- EXAMINING EXPRESSION ACROSS CELL TYPES ---- ####

cell_types <- seurat_obj$celltype_annotation
hMEs <- GetMEs(seurat_obj, harmonized = TRUE)
hMEs <- hMEs[,colnames(hMEs) != "grey"] #remove grey module


hME_by_type <- cbind(cell_type = cell_types, as.data.frame(hMEs)) %>%
  group_by(cell_type) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

library(pheatmap)

palette_length <- 100
my_color <- colorRampPalette(c("darkblue", "white", "red"))(palette_length)
my_breaks <- seq(min(hME_by_type), max(hME_by_type), length.out = palette_length)



custom_col_order <- c("T_cell_activation", "conv_T_cytokine_signalling", "all_T_cytokine_signalling", "TCR", "immune_specific_housekeeper1", "immune_specific_housekeeper2", "ubiquitous_housekeeper", "non_immune_housekeeper1", "non_immune_housekeeper2", "non_immune_housekeeper3")
custom_row_order <- c("CD4+T", "CD8+T", "TYPE_1_INNATE_T", "MATURE_B", "MONOCYTE_II_CCR2", "FIBROBLAST_I")

#default heatmap
pheatmap(hME_by_type,
         color = my_color,
         breaks = my_breaks,
         cluster_rows = TRUE, 
         angle_col = 45,
         cluster_cols = TRUE,
         border_color = NA,
         main = "Mean module eigengene per cell type")

pheatmap(hME_by_type[custom_row_order, custom_col_order],
         color = my_color,
         breaks = my_breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45,
         border_color = NA,
         main = "Mean module eigengene per cell type")

modules <- GetModules(seurat_obj)         

col_map <- c(T_cell_activation = "turquoise", 
             conv_T_cytokine_signalling = "blue", 
             TCR = "magenta",
             all_T_cytokine_signalling = "#E75480",
             immune_specific_housekeeper1 = "#C2B000",
             immune_specific_housekeeper2 = "#2E8B57",
             ubiquitous_housekeeper = "black",
             non_immune_housekeeper1 = "red",
             non_immune_housekeeper2 = "brown",
             non_immune_housekeeper3 = "purple")

library(ComplexHeatmap)

draw(
  Heatmap(
    hME_by_type[custom_row_order, custom_col_order],
    name = "Mean hME",
    column_title = "Mean hME across in vivo cell types",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_labels = names(col_map),
    column_names_gp = gpar(
      col = unname(col_map),  # <-- color each label according to col_map
      fontsize = 20
    ),
    row_names_gp = gpar(fontsize = 16),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 20, fontface = "bold")
  ),padding = unit(c(20, 25, 5, 5), "mm")  # top, right, bottom, left
  
)



T_mods <- c("TCR", "T_cell_activation", "conv_T_cytokine_signalling", "all_T_cytokine_signalling")
##firstly identify modules unique to CD4 vs CD8 T cells
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD8+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD4+T') %>% rownames
DME_cd4_vs_cd8 <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

rownames(DME_cd4_vs_cd8) <- gsub("-", "_", rownames(DME_cd4_vs_cd8))
DME_cd4_vs_cd8$module <- gsub("-", "_", DME_cd4_vs_cd8$module)

PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name = "conv_T_analysis",
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj, 
  #xlim_range = c(-4,4),
  DMEs,
  wgcna_name = 'conv_T',
  
)
#subset to modules unique to T-cells
DME_cd4_vs_cd8 <- DME_cd4_vs_cd8[T_mods,]



##innate vs CD8 T
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD8+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 't1_innate_T') %>% rownames
DME_t1_vs_CD8 <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)

rownames(DME_t1_vs_CD8) <- gsub("-", "_", rownames(DME_t1_vs_CD8))
DME_t1_vs_CD8$module <- gsub("-", "_", DME_t1_vs_CD8$module)

#subset to modules unique to T-cells
DME_t1_vs_CD8 <- DME_t1_vs_CD8[T_mods,]



##innate vs CD4
group1 <- seurat_obj@meta.data %>% subset(celltype_annotation == 'CD4+T') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(orig.ident == 't1_innate_T') %>% rownames
DME_t1_vs_CD4 <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox'
)
rownames(DME_t1_vs_CD4) <- gsub("-", "_", rownames(DME_t1_vs_CD4))
DME_t1_vs_CD4$module <- gsub("-", "_", DME_t1_vs_CD4$module)


#subset to modules unique to T-cells
DME_t1_vs_CD4 <- DME_t1_vs_CD4[T_mods,]

PlotDMEsVolcano(
  seurat_obj, 
  xlim_range = c(-15,15),
  DME_t1_vs_CD4,
  wgcna_name = 'conv_T_analysis',
  
)

p1 <- PlotDMEsVolcano(
  seurat_obj, 
  xlim_range = c(-15,15),
  DME_cd4_vs_cd8,
  wgcna_name = 'conv_T_analysis',
  label_size = 6,
  mod_point_size = 6
) + ggtitle(label = "CD4+ T (left) vs CD8+ T (right)") 
p2 <- PlotDMEsVolcano(
  seurat_obj, 
  xlim_range = c(-15,15),
  DME_t1_vs_CD8,
  wgcna_name = 'conv_T_analysis',
  label_size = 6,
  mod_point_size = 6
) + ggtitle(label = "Type 1 innate T (left) vs CD8+ T (right)") 

p3 <- PlotDMEsVolcano(
  seurat_obj, 
  xlim_range = c(-15,15),
  DME_t1_vs_CD4,
  wgcna_name = 'conv_T_analysis',
  label_size = 6,
  mod_point_size = 6
) + ggtitle(label = "Type 1 innate T (left) vs CD4+ T (right)") 
(p1 / p2 / p3)