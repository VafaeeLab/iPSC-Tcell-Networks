#in this script, we are visualising the TCR module eigengene expression across
#all of our iPSC datasets, once across the whole dataset and once subsetted
#to the most mature population

library(Seurat)
library(hdWGCNA)
library(patchwork)
library(ggplot2)
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


#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

g2gato_obj$phenotype <- g2gato_obj$maturity


#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_obj <- SetActiveWGCNA(EZ_obj, wgcna_name = "conv_T_mod_projection")

EZ_obj$phenotype <- recode(
  EZ_obj$phenotype,
  "mature_activated_T" = "mature_T_activated_subset"
)

#next ishi3d
ishi3d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_obj <- SetActiveWGCNA(ishi3d_obj, wgcna_name = "conv_T_mod_projection")

ishi3d_obj$phenotype <- ishi3d_obj$broad_celltype_revised


#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_obj <- SetActiveWGCNA(ishi2d_obj, wgcna_name = "conv_T_mod_projection")
ishi2d_obj$phenotype <- "mature_T"

#setting up color map for phenotypes
phenotype_levels <- unique(c(
  g2gato_obj$phenotype,
  ishi3d_obj$phenotype,
  ishi2d_obj$phenotype,
  EZ_obj$phenotype
))

label_size <- 8

#now visualising each module in each dataset
phenotype_colors <- make_color_map(phenotype_levels, palette = "Dark2")

#g2g_ato
g2g_ato_treatment <- DimPlot(g2gato_obj, group.by = "treatment", reduction = "harmony.umap", label = FALSE)
g2g_ato_phenotype <- DimPlot(g2gato_obj, 
                             group.by = "phenotype", 
                             reduction = "harmony.umap", 
                             label = TRUE,
                             repel = TRUE,
                             label.size = label_size,
                             cols = phenotype_colors)

g2g_ato_TCR <- ModuleFeaturePlot(g2gato_obj, module_names = "TCR", reduction = "harmony.umap")
g2g_ato_conv_t_cyt <- ModuleFeaturePlot(g2gato_obj, module_names = "conv_T_cytokine_signalling", reduction = "harmony.umap")
g2g_ato_all_t_cyt <- ModuleFeaturePlot(g2gato_obj, module_names = "all_T_cytokine_signalling", reduction = "harmony.umap")
g2g_ato_t_activation <- ModuleFeaturePlot(g2gato_obj, module_names = "T_cell_activation", reduction = "harmony.umap")

#EZ_T
EZ_T_treatment <- DimPlot(EZ_obj, group.by = "treatment", reduction = "harmony.umap", label = FALSE)
EZ_T_phenotype <- DimPlot(EZ_obj,  
                          group.by = "phenotype", 
                          reduction = "harmony.umap", 
                          label = TRUE,
                          repel = TRUE,
                          label.size = 6,
                          cols = phenotype_colors)

EZ_T_TCR <- ModuleFeaturePlot(EZ_obj, module_names = "TCR", reduction = "harmony.umap")
EZ_T_conv_t_cyt <- ModuleFeaturePlot(EZ_obj, module_names = "conv_T_cytokine_signalling", reduction = "harmony.umap")
EZ_T_all_t_cyt <- ModuleFeaturePlot(EZ_obj, module_names = "all_T_cytokine_signalling", reduction = "harmony.umap")
EZ_T_t_activation <- ModuleFeaturePlot(EZ_obj, module_names = "T_cell_activation", reduction = "harmony.umap")

#ishi3d
ishi_3d_phenotype <- DimPlot(ishi3d_obj,  
                             group.by = "phenotype",  
                             label = TRUE,
                             repel = TRUE,
                             label.size = 6,
                             cols = phenotype_colors)
ishi_3d_TCR <- ModuleFeaturePlot(ishi3d_obj, module_names = "TCR", reduction = "umap")
ishi_3d_conv_t_cyt <- ModuleFeaturePlot(ishi3d_obj, module_names = "conv_T_cytokine_signalling", reduction = "umap")
ishi_3d_all_t_cyt <- ModuleFeaturePlot(ishi3d_obj, module_names = "all_T_cytokine_signalling", reduction = "umap")
ishi_3d_t_activation <- ModuleFeaturePlot(ishi3d_obj, module_names = "T_cell_activation", reduction = "umap")

#ishi2d
ishi_2d_phenotype <- DimPlot(ishi2d_obj, 
                             group.by = "phenotype", 
                             label = TRUE,
                             repel = TRUE,
                             label.size = label_size,
                             cols = phenotype_colors)
ishi_2d_TCR <- ModuleFeaturePlot(ishi2d_obj, module_names = "TCR", reduction = "umap")
ishi_2d_conv_t_cyt <- ModuleFeaturePlot(ishi2d_obj, module_names = "conv_T_cytokine_signalling", reduction = "umap")
ishi_2d_all_t_cyt <- ModuleFeaturePlot(ishi2d_obj, module_names = "all_T_cytokine_signalling", reduction = "umap")
ishi_2d_t_activation <- ModuleFeaturePlot(ishi2d_obj, module_names = "T_cell_activation", reduction = "umap")


treatment_list <- list(
  "EZ_T" <- EZ_T_treatment,
  "g2g_ato" <- g2g_ato_treatment
)

phenotype_list <- list(
  "EZ_T" = EZ_T_phenotype,
  "g2g_ato" = g2g_ato_phenotype,
  "ishi_3d" = ishi_3d_phenotype,
  "ishi_2d" = ishi_2d_phenotype
)

TCR_list <- list(
  "EZ_T" = EZ_T_TCR,
  "g2g_ato" = g2g_ato_TCR,
  "ishi_3d" = ishi_3d_TCR,
  "ishi_2d" = ishi_2d_TCR
)

conv_t_cyt_list <- list(
  "EZ_T" = EZ_T_conv_t_cyt,
  "g2g_ato" = g2g_ato_conv_t_cyt,
  "ishi_3d" = ishi_3d_conv_t_cyt,
  "ishi_2d" = ishi_2d_conv_t_cyt
)

all_t_cyt_list <- list(
  "EZ_T" = EZ_T_all_t_cyt,
  "g2g_ato" = g2g_ato_all_t_cyt,
  "ishi_3d" = ishi_3d_all_t_cyt,
  "ishi_2d" = ishi_2d_all_t_cyt
)

t_activation_list <- list(
  "EZ_T" = EZ_T_t_activation,
  "g2g_ato" = g2g_ato_t_activation,
  "ishi_3d" = ishi_3d_t_activation,
  "ishi_2d" = ishi_2d_t_activation
)


themed_phenotype_list <- list()
for (dataset in names(phenotype_list)) {
  plot <- phenotype_list[[dataset]]
  themed_plot <- plot + 
    theme_void() +
    coord_fixed() +
    theme(#legend.position = "left",
          #legend.direction = "vertical",
          plot.title = element_text(hjust = 0.5,
                                    margin = margin(b = 20))) +
    NoLegend()
  themed_phenotype_list[[dataset]] <- themed_plot
}


#now visualising each module adjacent to a phenotype and treatment plot

#first making our top row, which plots EZ T

top <- wrap_plots(list(themed_phenotype_list[["EZ_T"]], 
                all_t_cyt_list[["EZ_T"]] + labs(title = "all T cytokine signalling"), 
                conv_t_cyt_list[["EZ_T"]]+ labs(title = "conventional T cytokine signalling"),
                TCR_list[["EZ_T"]]+ labs(title = "TCR") , 
                t_activation_list[["EZ_T"]]+ labs(title = "T cell activation")),
           ncol = 5) & theme(plot.title = element_text(size = 20))

top

top_mid <- wrap_plots(list(themed_phenotype_list[["g2g_ato"]]+ labs(title = NULL), 
                           all_t_cyt_list[["g2g_ato"]] + labs(title = NULL), 
                           conv_t_cyt_list[["g2g_ato"]]+ labs(title = NULL),
                           TCR_list[["g2g_ato"]]+ labs(title = NULL) , 
                           t_activation_list[["g2g_ato"]]+ labs(title = NULL)),
                      ncol = 5)

top_mid

bot_mid <- wrap_plots(list(themed_phenotype_list[["ishi_2d"]]+ labs(title = NULL), 
                           all_t_cyt_list[["ishi_2d"]] + labs(title = NULL), 
                           conv_t_cyt_list[["ishi_2d"]]+ labs(title = NULL),
                           TCR_list[["ishi_2d"]]+ labs(title = NULL) , 
                           t_activation_list[["ishi_2d"]]+ labs(title = NULL)),
                      ncol = 5)

bot <- wrap_plots(list(themed_phenotype_list[["ishi_3d"]]+ labs(title = NULL), 
                           all_t_cyt_list[["ishi_3d"]] + labs(title = NULL), 
                           conv_t_cyt_list[["ishi_3d"]]+ labs(title = NULL),
                           TCR_list[["ishi_3d"]]+ labs(title = NULL) , 
                           t_activation_list[["ishi_3d"]]+ labs(title = NULL)),
                      ncol = 5)



wrap_plots(list(top, top_mid, bot_mid, bot), nrow = 4)

top
top_mid
bot_mid
bot

all_t_cyt_list[["EZ_T"]] + 
  labs(title = NULL) +
  theme(legend.position = "none")

for (dataset in names(phenotype_list)) {
  plot <- wrap_plots(list(themed_phenotype_list[[dataset]] + labs(title = NULL), 
                  all_t_cyt_list[[dataset]] + labs(title = NULL), 
                  conv_t_cyt_list[[dataset]]+ labs(title = NULL),
                  TCR_list[[dataset]]+ labs(title = NULL), 
                  t_activation_list[[dataset]]+ labs(title = NULL)),
             ncol = 1)
  print(plot)
}


create_twin_plots <- function(list1, list2) {
  
  # Ensure the lists have the same names for proper matching
  if (!identical(sort(names(list1)), sort(names(list2)))) {
    stop("Input lists must have the same set of names.")
  }
  
  # Get the common names to iterate over
  plot_names <- names(list1)
  
  # Use lapply to iterate through each plot name and create the combined plot
  combined_plots_list <- lapply(plot_names, function(name) {
    
    # Access the corresponding plots from each list
    plot1 <- list1[[name]]
    plot2 <- list2[[name]]
    
    # Apply the specific theme to the first plot
    themed_plot1 <- plot1 + 
      theme_void() +
      coord_fixed() +
      theme(legend.position = "left",
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5,
                                      margin = margin(b = 20)))
    
    # Combine the plots
    combined_plot <- themed_plot1 + plot2
    
    # Use the '&' operator to apply coord_fixed() to BOTH subplots.
    # This is the key step to guarantee identical ratios and sizes.
    final_plot <- combined_plot
    
    return(combined_plot)
  })
  
  # Assign the original names to the new list of combined plots
  names(combined_plots_list) <- plot_names
  
  return(combined_plots_list)
}

create_triplet_plots <- function(list1, list2, list3) {
  
  # Ensure the lists have the same names for proper matching
  #if (!identical(sort(names(list1)), sort(names(list2))) || !identical(sort(names(list1)), sort(names(list3)))) {
  #  stop("Input lists must have the same set of names.")
  #}
  
  # Get the common names to iterate over
  plot_names <- names(list1)
  
  # Use lapply to iterate through each plot name and create the combined plot
  combined_plots_list <- lapply(plot_names, function(name) {
    
    # Access the corresponding plots from each list
    plot1 <- list1[[name]]
    plot2 <- list2[[name]]
    plot3 <- list3[[name]]
    
    # Apply the specific theme to the first two plots
    themed_plot1 <- plot1 + 
      theme_void() +
      coord_fixed() +
      theme(legend.position = "left",
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5,
                                      margin = margin(b = 20)))
    
    themed_plot2 <- plot2 + 
      theme_void() +
      coord_fixed() +
      theme(legend.position = "left",
            legend.direction = "vertical",
            plot.title = element_text(hjust = 0.5,
                                      margin = margin(b = 20)))
    
    # Combine the three plots side-by-side
    combined_plot <- themed_plot1 + themed_plot2 + plot3
    
    # Use the '&' operator to apply coord_fixed() to ALL subplots.
    # This is a key step to guarantee identical ratios and sizes.
    final_plot <- combined_plot & coord_fixed()
    
    return(final_plot)
  })
  
  # Assign the original names to the new list of combined plots
  names(combined_plots_list) <- plot_names
  
  return(combined_plots_list)
}

phenotype_TCR_list <- create_twin_plots(phenotype_list, TCR_list)
phenotype_TCR_list[1] 
phenotype_TCR_list[2]
phenotype_TCR_list[3] 
phenotype_TCR_list[4]

wrap_plots(phenotype_TCR_list, ncol = 2)

phenotype_all_t_cyt <- create_twin_plots(phenotype_list, all_t_cyt_list)
phenotype_all_t_cyt[1]
phenotype_all_t_cyt[2]
phenotype_all_t_cyt[3]
phenotype_all_t_cyt[4]

wrap_plots(phenotype_all_t_cyt, ncol = 2)

phenotype_conv_t_cyt <- create_twin_plots(phenotype_list, conv_t_cyt_list)
phenotype_conv_t_cyt[1]
phenotype_conv_t_cyt[2]
phenotype_conv_t_cyt[3]
phenotype_conv_t_cyt[4]

wrap_plots(c(phenotype_conv_t_cyt[1], phenotype_conv_t_cyt[2]), ncol = 2)

phenotype_activation <- create_twin_plots(phenotype_list, t_activation_list)
phenotype_activation[1]
phenotype_activation[2]
phenotype_activation[3]
phenotype_activation[4]

wrap_plots(phenotype_activation[c(1,2,4)], ncol = 2) #omitting ishi 3d
