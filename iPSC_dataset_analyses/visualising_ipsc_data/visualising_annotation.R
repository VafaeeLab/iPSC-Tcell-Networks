library(Seurat)
library(hdWGCNA)
library(patchwork)
library(ggplot2)
library(dplyr)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

#add a phenotype field that just copies maturity
g2gato_obj$phenotype <- g2gato_obj$maturity

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_obj$phenotype <- recode(
  EZ_obj$phenotype,
  "mature_activated_T" = "mature_T_activated_subset"
)

EZ_obj$treatment <- recode(
  EZ_obj$treatment,
  "EZ" = "EZH1-KD (pre-activation)",
  "EZ_activated" = "EZH1-KD (post-activation)"
)


#next ishi3d
ishi3d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_obj$phenotype <- ishi3d_obj$broad_celltype_revised

#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_obj$phenotype <- "mature_T"


##Creating colour maps for our various fields

# Define colors manually or use a palette
library(scales)
library(RColorBrewer)

# Function to make a color map for a metadata field
make_color_map <- function(values, palette = "Set2") {
  values <- unique(values)
  n <- length(values)
  cols <- brewer.pal(min(n, brewer.pal.info[palette, "maxcolors"]), palette)
  
  # If more groups than colors, recycle palette
  if (n > length(cols)) cols <- rep(cols, length.out = n)
  
  setNames(cols, values)
}

make_custom_color_map <- function(values, palette) {
  values <- unique(values)
  n <- length(values)
  
  # Use the custom palette
  cols <- palette
  
  # If more groups than colors, recycle the palette
  if (n > length(cols)) {
    cols <- rep(cols, length.out = n)
  } else {
    cols <- cols[1:n]
  }
  
  setNames(cols, values)
}

# Get all unique levels across all objects you’ll plot
celltypist_levels <- unique(c(
  g2gato_obj$celltypist_high_predictions,
  ishi3d_obj$celltypist_high_predictions,
  ishi2d_obj$celltypist_high_predictions,
  EZ_obj$celltypist_high_predictions
))

phenotype_levels <- unique(c(
  g2gato_obj$phenotype,
  ishi3d_obj$phenotype,
  ishi2d_obj$phenotype,
  EZ_obj$phenotype
))

treatment_levels <- unique(c(
  g2gato_obj$treatment,
  EZ_obj$treatment
))

# 1. Define your custom color palette
celltypist_palette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
                        "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
                        "#cab2d6", "#6a3d9a", "#FFD700", "#b15928",
                        "#8dd3c7", "#bebada") 

# Create a named vector
celltypist_colors <- make_custom_color_map(celltypist_levels, palette = celltypist_palette)
phenotype_colors <- make_color_map(phenotype_levels, palette = "Dark2")
treatment_colors <- make_color_map(treatment_levels, palette = "Set1")

label_size <- 9

#MAKING G2G ATO PLOTs
G2G_ato_plots <- list()

G2G_ato_plots[["treatment"]] <- DimPlot(g2gato_obj, 
                                        group.by = "treatment",
                                        label = FALSE,
                                        reduction = "harmony.umap",
                                        cols = treatment_colors) + 
  ggtitle(label = "Sumanaweera_3D") +
  theme(legend.position = "top", legend.text = element_text(size = 18))

G2G_ato_plots[["celltypist"]] <- DimPlot(g2gato_obj, 
        group.by = "celltypist_high_predictions",
        label = TRUE,
        repel = TRUE,
        label.size = label_size,
        label.box = FALSE,
        reduction = "harmony.umap",
        cols = celltypist_colors) + 
  ggtitle(label = "Sumanaweera_3D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))




G2G_ato_plots[["clusters"]] <- DimPlot(g2gato_obj, 
                                       label = FALSE,
                                       group.by = "seurat_clusters",
                                       reduction = "harmony.umap")+
  theme(legend.position = "top", legend.text = element_text(size = 18))


G2G_ato_plots[["phenotype"]] <- DimPlot(g2gato_obj, 
                                        group.by = "phenotype",
                                        label = TRUE,
                                        repel = TRUE,
                                        label.size = label_size,
                                        label.box = FALSE,
                                        reduction = "harmony.umap",
                                        cols = phenotype_colors) + 
  ggtitle(label = "Sumanaweera_3D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))




##MAKING EZ T CELLS PLOTS
EZ_plots <- list()

EZ_plots[["treatment"]] <- DimPlot(EZ_obj, 
                                   group.by = "treatment",
                                   label = FALSE,
                                   reduction = "harmony.umap",
                                   cols = treatment_colors) + 
  ggtitle(label = "Jing_2D") +
  theme(legend.position = "top", legend.text = element_text(size = 18))

EZ_plots[["celltypist"]] <- DimPlot(EZ_obj, 
                                    group.by = "celltypist_high_predictions",
                                    label = TRUE,
                                    repel = TRUE,
                                    label.size = label_size,
                                    label.box = FALSE,
                                    reduction = "harmony.umap",
                                    cols = celltypist_colors)  + 
  ggtitle("Jing_2D") +
  theme(legend.position = "top", legend.text = element_text(size = 15))




EZ_plots[["clusters"]] <- DimPlot(EZ_obj, 
                                       label = FALSE,
                                       group.by = "seurat_clusters",
                                       reduction = "harmony.umap")+
  theme(legend.position = "top", legend.text = element_text(size = 18))


EZ_plots[["phenotype"]] <- DimPlot(EZ_obj, 
                                        group.by = "phenotype",
                                        label = TRUE,
                                        repel = TRUE,
                                        label.size = label_size,
                                        label.box = FALSE,
                                        reduction = "harmony.umap",
                                        cols = phenotype_colors) + 
  ggtitle(label = "Jing_2D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))






##MAKING Ishi2d CELLS PLOTS
#ishi2d_obj <- FindNeighbors(ishi2d_obj)
#ishi2d_obj <- FindClusters(ishi2d_obj, resolution = 0.1)

ishi2d_plots <- list()

ishi2d_plots[["treatment"]] <- DimPlot(ishi2d_obj, 
                                       group.by = "orig.ident",
                                       label = FALSE,
                                       repel = TRUE,
                                       label.size = label_size,
                                       label.box = FALSE,
                                       reduction = "umap") + 
  ggtitle(label = "Ishiguro_2D") +
  NoLegend()

ishi2d_plots[["celltypist"]] <- DimPlot(ishi2d_obj,
                                        group.by = "celltypist_high_predictions",
                                        label = TRUE,
                                        repel = TRUE,
                                        label.size = label_size,
                                        label.box = FALSE,
                                        reduction = "umap",
                                        cols = celltypist_colors) + 
  ggtitle("Ishiguro_2D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))




ishi2d_plots[["clusters"]] <- DimPlot(ishi2d_obj, 
                                  label = FALSE,
                                  group.by = "seurat_clusters",
                                  reduction = "umap")+
  theme(legend.position = "top", legend.text = element_text(size = 18))


ishi2d_plots[["phenotype"]] <- DimPlot(ishi2d_obj, 
                                   group.by = "phenotype",
                                   label = TRUE,
                                   repel = TRUE,
                                   label.size = label_size,
                                   label.box = FALSE,
                                   reduction = "umap",
                                   cols = phenotype_colors) + 
  ggtitle("Ishiguro_2D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))





##ISHI 3D plots
ishi3d_plots <- list()

ishi3d_plots[["treatment"]] <- DimPlot(ishi3d_obj, 
                                       group.by = "orig.ident",
                                       label = FALSE,
                                       repel = TRUE,
                                       label.size = label_size,
                                       label.box = FALSE,
                                       reduction = "umap") + 
  ggtitle(label = "Ishiguro_3D") +
  NoLegend()

ishi3d_plots[["celltypist"]] <- DimPlot(ishi3d_obj, 
                                        group.by = "celltypist_high_predictions",
                                        label = TRUE,
                                        repel = TRUE,
                                        label.size = label_size,
                                        label.box = FALSE,
                                        reduction = "umap",
                                        cols = celltypist_colors) + 
  ggtitle("Ishiguro_3D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))




ishi3d_plots[["clusters"]] <- DimPlot(ishi3d_obj, 
                                      label = FALSE,
                                      group.by = "seurat_clusters",
                                      reduction = "umap")+
  theme(legend.position = "top", legend.text = element_text(size = 18))


ishi3d_plots[["phenotype"]] <- DimPlot(ishi3d_obj,
                                       group.by = "phenotype",
                                       label = TRUE,
                                       repel = TRUE,
                                       label.size = label_size,
                                       label.box = FALSE,
                                       reduction = "umap",
                                       cols = phenotype_colors) + 
  ggtitle("Ishiguro_3D")+
  theme(legend.position = "top", legend.text = element_text(size = 18))

## Making treatment plots
G2G_ato_plots[["treatment"]]
EZ_plots[["treatment"]]
ishi2d_plots[["treatment"]]
ishi3d_plots[["treatment"]]

wrap_plots(list(G2G_ato_plots[["treatment"]], EZ_plots[["treatment"]],ishi2d_plots[["treatment"]], ishi3d_plots[["treatment"]]), ncol = 2)


## Making celltypist plots
G2G_ato_plots[["celltypist"]]
EZ_plots[["celltypist"]]
ishi2d_plots[["celltypist"]]
ishi3d_plots[["celltypist"]]

wrap_plots(list(EZ_plots[["celltypist"]],G2G_ato_plots[["celltypist"]],ishi2d_plots[["celltypist"]], ishi3d_plots[["celltypist"]]), ncol = 2) & NoLegend()

#plus a combined legend
# Create a dummy dataframe that includes all categories
dummy_df <- data.frame(
  celltype = factor(names(celltypist_colors), levels = names(celltypist_colors)),
  x = 1, y = 1
)


# Make a plot that uses these colors
legend_plot <- ggplot(dummy_df, aes(x, y, color = celltype)) +
  geom_point(size = 8) +  # or geom_bar() if you prefer
  scale_color_manual(values = celltypist_colors) +
  theme_void() +  # remove axes and background
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.direction = "horizontal"
  )

legend_only <- cowplot::get_legend(legend_plot)
# Display the legend
cowplot::plot_grid(legend_only)


##Making phenotype plots
G2G_ato_plots[["phenotype"]] + NoLegend()
EZ_plots[["phenotype"]] + NoLegend()
ishi2d_plots[["phenotype"]] + NoLegend()
ishi3d_plots[["phenotype"]] + NoLegend()
#plus a combined legend
# Create a dummy dataframe that includes all categories
dummy_df <- data.frame(
  phenotype = factor(names(phenotype_colors), levels = names(phenotype_colors)),
  x = 1, y = 1
)


####making supplementary material plots
wrap_plots(G2G_ato_plots, nrow = 2)
wrap_plots(EZ_plots, nrow = 2)
wrap_plots(ishi2d_plots, nrow = 2)
wrap_plots(ishi3d_plots, nrow = 2)


# Make a plot that uses these colors
legend_plot <- ggplot(dummy_df, aes(x, y, color = phenotype)) +
  geom_point(size = 5) +  # or geom_bar() if you prefer
  scale_color_manual(values = phenotype_colors) +
  theme_void() +  # remove axes and background
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

legend_only <- cowplot::get_legend(legend_plot)
# Display the legend
cowplot::plot_grid(legend_only)

wrap_plots(ishi3d_plots[1:2], nrow = 1) &
  theme(
    legend.position = "top",   # top-right inside plot (x, y from 0 to 1)
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",     
    legend.key.size = unit(0.4, "cm"), # smaller boxes
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11)
  ) 

G2G_ato_plots[[4]] +
  ggtitle(label = "Sumanaweera-3D: manual annotations") +
  theme(
    legend.position = "top",
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",     
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

wrap_plots(c(EZ_plots[[1]], G2G_ato_plots[[1]], ishi2d_plots[[1]], ishi3d_plots[[1]]), nrow = 2)

# Create the master scale object
shared_scale <- scale_color_manual(
  name = "CellTypist predictions",  # Or whatever you want the legend title to be
  values = celltypist_palette,
  limits = celltypist_levels 
)


plot1 <- EZ_plots[["celltypist"]] + shared_scale
plot2 <- G2G_ato_plots[["celltypist"]] + shared_scale
plot3 <- ishi2d_plots[["celltypist"]] + shared_scale
plot4 <- ishi3d_plots[["celltypist"]] + shared_scale

# Now, combining them will produce a single, de-duplicated legend
plot1 + plot2 + plot3 + plot4 + 
  plot_layout(guides = "collect")

EZ_plots[["celltypist"]] + G2G_ato_plots[["celltypist"]] + ishi2d_plots[["celltypist"]] + ishi3d_plots[["celltypist"]] & 
  NoLegend() 




EZ_plots[["phenotype"]] + G2G_ato_plots[["phenotype"]] + ishi2d_plots[["phenotype"]] + ishi3d_plots[["phenotype"]] & 
  NoLegend() 


# Create a dummy dataframe that includes all categories
dummy_df <- data.frame(
  phenotype = factor(names(phenotype_colors), levels = names(phenotype_colors)),
  x = 1, y = 1
)


# Make a plot that uses these colors
legend_plot <- ggplot(dummy_df, aes(x, y, color = phenotype)) +
  geom_point(size = 4) +  # or geom_bar() if you prefer
  scale_color_manual(values = phenotype_colors) +
  theme_void() +  # remove axes and background
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

legend_only <- cowplot::get_legend(legend_plot)

# Display the legend
cowplot::plot_grid(legend_only)

EZ_plots[[4]] +
  ggtitle(label = "Jing-2D: manual annotations") +
  theme(
    legend.position = "top",
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",     
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11)
  )

ishi2d_plots[[3]] +
  ggtitle(label = "Ishiguro-2D: manual annotations") +
  theme(
    legend.position = "top",
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",     
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

ishi3d_plots[[3]] +
  ggtitle(label = "Ishiguro-3D: manual annotations") +
  theme(
    legend.position = "top",
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",     
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
