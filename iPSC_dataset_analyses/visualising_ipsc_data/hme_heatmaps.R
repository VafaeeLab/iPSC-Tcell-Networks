#in this script we are measuring the correlations of the modules across our
#ipsc datasets
library(Seurat)
library(hdWGCNA)
library(ComplexHeatmap)
library(circlize)

#sets the gradient for the heatmap vals from -1 to 1
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

#col map
col_map <- c(T_cell_activation = "turquoise",
             TCR = "magenta",
             all_T_cytokine_signalling = "#E75480",
             conv_T_cytokine_signalling = "blue")


#getting the cor matrix from the EZ T dataset
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

T_mods <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling") #omitting non-preserved mods

#extracting hME vals
eig <- GetMEs(seurat_obj)
#subsetting to preserved mods
eig <- eig[,colnames(eig) %in% T_mods]

ez_t_cor_mat <- cor(eig, method = "pearson") 

#getting the cor matrix from the EZ T dataset
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

T_mods <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling", "conv_T_cytokine_signalling") #omitting non-preserved mods

#extracting hME vals
eig <- GetMEs(seurat_obj)
#subsetting to preserved mods
eig <- eig[,colnames(eig) %in% T_mods]

g2gato_cor_mat <- cor(eig, method = "pearson") 


#getting the cor matrix from the ishi2d dataset
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

T_mods <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling") #omitting non-preserved mods

#extracting hME vals
eig <- GetMEs(seurat_obj)
#subsetting to preserved mods
eig <- eig[,colnames(eig) %in% T_mods]

ishi2d_cor_mat <- cor(eig, method = "pearson") 


#getting the cor matrix from the ishi3d dataset
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

T_mods <- c( "conv_T_cytokine_signalling", "all_T_cytokine_signalling") #omitting non-preserved mods

#extracting hME vals
eig <- GetMEs(seurat_obj)
#subsetting to preserved mods
eig <- eig[,colnames(eig) %in% T_mods]

ishi3d_cor_mat <- cor(eig, method = "pearson") 


#use this to copy/paste into custom orders for each heatmap
custom_col_order <- c("T_cell_activation", "TCR", "conv_T_cytokine_signalling", "all_T_cytokine_signalling")
custom_row_order <- c("T_cell_activation", "TCR", "conv_T_cytokine_signalling", "all_T_cytokine_signalling")

##Creating heatmaps

# Create the EZ_T/Jing2D heatmap

custom_col_order <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling")
custom_row_order <- c("T_cell_activation", "TCR","all_T_cytokine_signalling")

draw(
  Heatmap(
    ez_t_cor_mat[custom_row_order, custom_col_order],
    name = "Jing_2D hME Correlation",
    col = col_fun,
    rect_gp = gpar(col = "grey55", lwd = 0.4),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    column_title = "Jing_2D hME correlation",
    
    # Font and layout styling (matching your example)
    column_names_gp = gpar(fontsize = 20,
                           col = unname(col_map)),
    row_names_gp = gpar(fontsize = 20,
                        col = unname(col_map)),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    show_heatmap_legend = FALSE
  ),
  padding = unit(c(5, 25, 5, 20), "mm")  # top, right, bottom, left
)


# Create the g2gato/sumanaweera3d heatmap

custom_col_order <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling", "conv_T_cytokine_signalling")
custom_row_order <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling", "conv_T_cytokine_signalling")

draw(
  Heatmap(
    g2gato_cor_mat[custom_row_order, custom_col_order],
    name = "Sumanaweera_3D hME Correlation",
    col = col_fun,
    rect_gp = gpar(col = "grey55", lwd = 0.4),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    column_title = "Sumanaweera_3D hME correlation",
    
    # Font and layout styling (matching your example)
    column_names_gp = gpar(fontsize = 18,
                           col = unname(col_map)),
    row_names_gp = gpar(fontsize = 18,
                        col = unname(col_map)),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    show_heatmap_legend = FALSE
  ),
  padding = unit(c(5, 25, 5, 25), "mm")  # top, right, bottom, left
)

# Create the ishi2d heatmap

custom_col_order <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling")
custom_row_order <- c("T_cell_activation", "TCR", "all_T_cytokine_signalling")

draw(
  Heatmap(
    ishi2d_cor_mat[custom_row_order, custom_col_order],
    name = "Ishiguro_2D hME Correlation",
    col = col_fun,
    rect_gp = gpar(col = "grey55", lwd = 0.4),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    column_title = "Ishiguro_2D hME correlation",
    
    # Font and layout styling (matching your example)
    column_names_gp = gpar(fontsize = 20,
                           col = unname(col_map)),
    row_names_gp = gpar(fontsize = 20,
                        col = unname(col_map)),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    show_heatmap_legend = FALSE
  ),
  padding = unit(c(5, 25, 5, 20), "mm")  # top, right, bottom, left
)

# Create the ishi3d heatmap

custom_col_order <- c("all_T_cytokine_signalling", "conv_T_cytokine_signalling")
custom_row_order <- c("all_T_cytokine_signalling", "conv_T_cytokine_signalling")

draw(
  Heatmap(
    ishi3d_cor_mat[custom_row_order, custom_col_order],
    name = "Ishiguro_3D hME Correlation",
    col = col_fun,
    rect_gp = gpar(col = "grey55", lwd = 0.4),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    column_title = "Ishiguro_3D hME correlation",
    
    # Font and layout styling (matching your example)
    column_names_gp = gpar(fontsize = 18,
                           col = unname(col_map)),
    row_names_gp = gpar(fontsize = 20,
                        col = unname(col_map)),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 18, fontface = "bold"),
    show_heatmap_legend = FALSE
  ),
  padding = unit(c(5, 30, 5, 25), "mm")  # top, right, bottom, left
)



# Create a legend object only once
# Create a horizontal legend with larger dimensions
legend_obj <- Legend(
  col_fun = col_fun,
  title = "hME Correlation",
  at = seq(-1, 1, by = 0.5),
  labels = seq(-1, 1, by = 0.5),
  direction = "horizontal",       # ← makes it horizontal
  legend_width = unit(8, "cm"),   # ← makes the colour bar longer
  grid_height = unit(0.8, "cm"),  # ← makes the colour bar thicker
  title_gp = gpar(fontsize = 20, fontface = "bold"),
  labels_gp = gpar(fontsize = 16)
)

# Draw it by itself as a separate plot
grid.newpage()
draw(
  legend_obj,
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  just = "center"
)
