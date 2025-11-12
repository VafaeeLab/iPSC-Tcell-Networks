#in this script we are going to load up each iPSC obj sequentially
#and extract the module preservation statistics (using most mature cells each time)

#NOTE: We used gemini here to create a custom method for generating plots
#but in retrospect, it would actually just be easier to generate plot_lists
#from each obj one at a time, the second index in each list is our pres stats
#dumb dumb

library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(tibble)

#starting with ISHI 3D
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_pres <- GetModulePreservation(seurat_obj, mod_name = "conv_T_preservation_lymphoid_only")

#now ISHI 2d
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_pres <- GetModulePreservation(seurat_obj, mod_name = "conv_T_preservation")

#now EZ T
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_T_pres <- GetModulePreservation(seurat_obj, mod_name = "conv_T_preservation_mature", wgcna_name = "conv_T_mod_projection_mature")

#lasly g2g_ato
seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated_NO_TNF.rds")

g2g_ato_pres <- GetModulePreservation(seurat_obj, mod_name = "conv_T_preservation_mature")


## NOW we extract data and plots ##

PlotModulePreservationCustom <- function(
    zsummary_qual_df,
    zsummary_pres_list,
    module_info_df,
    plot_labels = TRUE,
    label_size = 4,
    mod_point_size = 4,
    modules_to_plot = NULL  # ← new argument
){
  
  # zsummary_qual_df: A dataframe with 'module' and 'Zsummary.qual' columns.
  # zsummary_pres_list: A named list of dataframes. Each dataframe should have 'module' and 'Zsummary.pres' columns.
  # module_info_df: A dataframe with 'module', 'color', and 'moduleSize' columns.
  # modules_to_plot: Optional vector of module names to include (e.g., c("TCR", "all_T_cytokine_signalling")).
  
  plot_list <- list()
  
  # Filter for specified modules (if provided)
  if(!is.null(modules_to_plot)){
    zsummary_qual_df <- zsummary_qual_df %>% filter(module %in% modules_to_plot)
    module_info_df  <- module_info_df %>% filter(module %in% modules_to_plot)
    zsummary_pres_list <- lapply(zsummary_pres_list, function(df) {
      df %>% filter(module %in% modules_to_plot)
    })
  }
  
  # 1. Generate the plot for Zsummary.qual
  #-----------------------------------------
  
  plot_df_qual <- inner_join(zsummary_qual_df, module_info_df, by = "module") %>%
    subset(!(module %in% c('grey', 'gold')))
  
  qual_p <- plot_df_qual %>% 
    ggplot(aes(x = moduleSize, y = Zsummary.qual, fill = module, color = module)) +
    geom_rect(
      data = plot_df_qual[1,],
      aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 2), 
      fill = 'grey75', alpha = 0.8, color = NA
    ) +
    geom_rect(
      data = plot_df_qual[1,],
      aes(xmin = 0, xmax = Inf, ymin = 2, ymax = 10), 
      fill = 'grey92', alpha = 0.8, color = NA
    ) +
    geom_point(size = mod_point_size, pch = 21, color = 'black') +
    scale_fill_manual(values = plot_df_qual$color) +
    scale_color_manual(values = plot_df_qual$color) +
    scale_x_continuous(trans = 'log10') +
    ylab("Zsummary.qual") +
    xlab("Module Size") +
    ggtitle("Zsummary.qual") +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(plot_labels){
    qual_p <- qual_p + ggrepel::geom_text_repel(
      label = plot_df_qual$module, 
      size = label_size, 
      max.overlaps = Inf, 
      color = 'black'
    )
  }
  
  plot_list[['Zsummary.qual']] <- qual_p
  
  # 2. Loop through and plot each Zsummary.pres
  #----------------------------------------------
  
  for(name in names(zsummary_pres_list)){
    cur_zsummary_pres <- zsummary_pres_list[[name]]
    
    plot_df_pres <- inner_join(cur_zsummary_pres, module_info_df, by = "module") %>%
      subset(!(module %in% c('grey', 'gold')))
    
    pres_p <- plot_df_pres %>% 
      ggplot(aes(x = moduleSize, y = Zsummary.pres, fill = module, color = module)) +
      geom_rect(
        data = plot_df_pres[1,],
        aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 2), 
        fill = 'grey75', alpha = 0.8, color = NA
      ) +
      geom_rect(
        data = plot_df_pres[1,],
        aes(xmin = 0, xmax = Inf, ymin = 2, ymax = 10), 
        fill = 'grey92', alpha = 0.8, color = NA
      ) +
      geom_point(size = mod_point_size, pch = 21, color = 'black') +
      scale_fill_manual(values = plot_df_pres$color) +
      scale_color_manual(values = plot_df_pres$color) +
      scale_x_continuous(trans = 'log10') +
      ylab("Zsummary.pres") +
      xlab("Module Size") +
      ggtitle(name) +
      NoLegend() +
      theme(plot.title = element_text(hjust = 0.5))
    
    if(plot_labels){
      pres_p <- pres_p + ggrepel::geom_text_repel(
        label = plot_df_pres$module, 
        size = label_size, 
        max.overlaps = Inf, 
        color = 'black'
      )
    }
    
    plot_list[[name]] <- pres_p
  }
  
  # Return the list of ggplot objects
  if(length(plot_list) == 1){
    return(plot_list[[1]])
  }
  
  return(plot_list)
}

z_df_ref <- EZ_T_pres$Z
zsummary_qual_df <- data.frame(
  module = rownames(z_df_ref),
  Zsummary.qual = z_df_ref$Zsummary.qual
)

# Extract Z dataframes from each query
G2G_ATO_df1 <- g2g_ato_pres$Z
EZ_T_df2 <- EZ_T_pres$Z
Ishi_2D_df3 <- ishi2d_pres$Z
Ishi_3D_df4 <- ishi3d_pres$Z

# Create the list
zsummary_pres_list <- list(
  Jing_2D = data.frame(module = rownames(EZ_T_df2), Zsummary.pres = EZ_T_df2$Zsummary.pres),
  Sumanaweera_3D = data.frame(module = rownames(G2G_ATO_df1), Zsummary.pres = G2G_ATO_df1$Zsummary.pres),
  Ishiguro_2D = data.frame(module = rownames(Ishi_2D_df3), Zsummary.pres = Ishi_2D_df3$Zsummary.pres),
  Ishiguro_3D = data.frame(module = rownames(Ishi_3D_df4), Zsummary.pres = Ishi_3D_df4$Zsummary.pres)
)

# Get module assignments and colors from the reference object
modules <- GetModules(seurat_obj)
module_colors <- modules %>% 
  select(module, color) %>% 
  distinct()

# Create the module info dataframe
module_info_df <- data.frame(
  module = rownames(z_df_ref),
  moduleSize = z_df_ref$moduleSize
)
module_info_df <- inner_join(module_info_df, module_colors, by = "module")

# Generate the plots
plot_list <- PlotModulePreservationCustom(
  zsummary_qual_df = zsummary_qual_df,
  zsummary_pres_list = zsummary_pres_list,
  mod_point_size = 5,
  label_size = 5,
  module_info_df = module_info_df,
  modules_to_plot = c("TCR", "T_cell_activation", "conv_T_cytokine_signalling", "all_T_cytokine_signalling")
)
plot_list$Zsummary.qual
plot_list$G2G_ATO


# Access your plots by name from the list
qual_plot <- plot_list$Zsummary.qual
pres_plot1 <- plot_list$Jing_2D
pres_plot2 <- plot_list$Sumanaweera_3D
pres_plot3 <- plot_list$Ishiguro_2D
pres_plot4 <- plot_list$Ishiguro_3D

# Create the 2x2 grid of preservation plots
preservation_grid <- (pres_plot1 | pres_plot2) / (pres_plot3 | pres_plot4)

preservation_grid

# Combine the quality plot on top of the grid
combined_plot <- qual_plot / preservation_grid

combined_plot
preservation_grid

##Now plotting preservation of specific modules
# Extract Zsummary.pres, keeping module names from rownames
df_list <- list(
  Jing_2D = data.frame(module = rownames(EZ_T_df2), Zsummary.pres = EZ_T_df2$Zsummary.pres, dataset = "Jing_2D"),
  Sumanaweera_3D = data.frame(module = rownames(G2G_ATO_df1), Zsummary.pres = G2G_ATO_df1$Zsummary.pres, dataset = "Sumanaweera_3D"),
  Ishiguro_2D = data.frame(module = rownames(Ishi_2D_df3), Zsummary.pres = Ishi_2D_df3$Zsummary.pres, dataset = "Ishiguro_2D"),
  Ishiguro_3D = data.frame(module = rownames(Ishi_3D_df4), Zsummary.pres = Ishi_3D_df4$Zsummary.pres, dataset = "Ishiguro_3D")
  )

all_pres <- bind_rows(df_list)

modules_of_interest <- c("TCR", "conv_T_cytokine_signalling", "all_T_cytokine_signalling", "T_cell_activation")

plot_data <- all_pres %>%
  filter(module %in% modules_of_interest)

ggplot(plot_data, aes(x = module, y = Zsummary.pres, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  labs(
    title = "Comparative Module Preservation (Zsummary.pres)",
    x = "Module",
    y = "Zsummary.pres"
  ) +
  theme_classic()

ggplot(plot_data, aes(x = dataset, y = Zsummary.pres)) +
  geom_point(aes(color = dataset), size = 3) +
  geom_line(aes(group = module), linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  facet_wrap(~module, scales = "free_y") +
  labs(
    title = "Comparative Module Preservation Across Datasets",
    x = "Dataset",
    y = "Zsummary.pres"
  ) +
  theme_bw()

# Pick your module and bar colour
module_of_interest <- "TCR"
bar_colour <- "magenta"   # <-- change to your module's WGCNA color

plot_data <- all_pres %>%
  filter(module == module_of_interest)


# Manually set the desired order
desired_order <- c("Jing_2D", "Sumanaweera_3D",  "Ishiguro_2D", "Ishiguro_3D")

# Convert the 'dataset' column to a factor with the specified order
plot_data$dataset <- factor(plot_data$dataset, levels = desired_order)

ggplot(plot_data, aes(x = dataset, y = Zsummary.pres)) +
  geom_col(fill = bar_colour) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  labs(
    title = paste("Preservation of", module_of_interest),
    x = "iPSC Dataset",
    y = "Zsummary.pres"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18),   # increase x-axis label size
    axis.text.y = element_text(size = 16),   # optional: increase y-axis label size
    axis.title.x = element_text(size = 18),  # optional: increase x-axis title size
    axis.title.y = element_text(size = 18)   # optional: increase y-axis title size
  )

# T cell activation
module_of_interest <- "T_cell_activation"
bar_colour <- "turquoise"   # <-- change to your module's WGCNA color

plot_data <- all_pres %>%
  filter(module == module_of_interest)

# Convert the 'dataset' column to a factor with the specified order
plot_data$dataset <- factor(plot_data$dataset, levels = desired_order)

ggplot(plot_data, aes(x = dataset, y = Zsummary.pres)) +
  geom_col(fill = bar_colour) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  labs(
    title = paste("Preservation of", module_of_interest),
    x = "iPSC Dataset",
    y = "Zsummary.pres"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18),   # increase x-axis label size
    axis.text.y = element_text(size = 16),   # optional: increase y-axis label size
    axis.title.x = element_text(size = 18),  # optional: increase x-axis title size
    axis.title.y = element_text(size = 18)   # optional: increase y-axis title size
  )

# all t cytokine signalling
module_of_interest <- "all_T_cytokine_signalling"
bar_colour <- "pink"   # <-- change to your module's WGCNA color

plot_data <- all_pres %>%
  filter(module == module_of_interest)

# Convert the 'dataset' column to a factor with the specified order
plot_data$dataset <- factor(plot_data$dataset, levels = desired_order)

ggplot(plot_data, aes(x = dataset, y = Zsummary.pres)) +
  geom_col(fill = bar_colour) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  labs(
    title = paste("Preservation of", module_of_interest),
    x = "iPSC Dataset",
    y = "Zsummary.pres"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18),   # increase x-axis label size
    axis.text.y = element_text(size = 16),   # optional: increase y-axis label size
    axis.title.x = element_text(size = 18),  # optional: increase x-axis title size
    axis.title.y = element_text(size = 18)   # optional: increase y-axis title size
  )

# conv t cytokine signalling
module_of_interest <- "conv_T_cytokine_signalling"
bar_colour <- "blue"   # <-- change to your module's WGCNA color

plot_data <- all_pres %>%
  filter(module == module_of_interest)

# Convert the 'dataset' column to a factor with the specified order
plot_data$dataset <- factor(plot_data$dataset, levels = desired_order)

ggplot(plot_data, aes(x = dataset, y = Zsummary.pres)) +
  geom_col(fill = bar_colour) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  labs(
    title = paste("Preservation of", module_of_interest),
    x = "iPSC Dataset",
    y = "Zsummary.pres"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18),   # increase x-axis label size
    axis.text.y = element_text(size = 16),   # optional: increase y-axis label size
    axis.title.x = element_text(size = 18),  # optional: increase x-axis title size
    axis.title.y = element_text(size = 18)   # optional: increase y-axis title size
  )
