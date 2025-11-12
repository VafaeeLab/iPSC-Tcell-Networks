#in this script, we are visualising our network modules for the sake of our
#figures
library(Seurat)
library(hdWGCNA)
library(ggplot2)

setwd("~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/cont_T_analysis")

seurat_obj <- readRDS(file = "~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/annotated_reference_obj.rds")




#hub gene netowrk plot
# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

#we get weird bugs unless we change this setting
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB limit

#setting edge cutoff
cutoff <- 0.12
max <- 0.8
min <- 0.2

# hubgene network for T cell activation
T_activation_g <- CustomNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=0,
  edge_prop = 1,
  mods = "T_cell_activation",
  vertex.label.cex = 0.8,
  edge_cutoff = cutoff,
  return_graph = TRUE
)

kme_T_activation <- modules[, c("gene_name", "kME_T_cell_activation")]

T_activation_plot <- make_network_plot(T_activation_g, 
                                       kme_T_activation, 
                                       "kME_T_cell_activation", 
                                       "turquoise", 
                                       kME_min = min, 
                                       kME_max = max)


# hubgene network plot for TCR
TCR_g <- CustomNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=0,
  edge_prop = 1,
  mods = "TCR",
  vertex.label.cex = 0.8,
  edge_cutoff = cutoff,
  return_graph = TRUE
)

kME_TCR <- modules[, c("gene_name", "kME_TCR")]

TCR_plot <- make_network_plot(TCR_g, 
                              kME_TCR, 
                              "kME_TCR", 
                              "magenta",
                              kME_min = min, 
                              kME_max = max)



# hubgene network plot for conv T cytokine signalling
conv_T_cytokine_g <- CustomNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=0,
  edge_prop = 1,
  mods = "conv_T_cytokine_signalling",
  vertex.label.cex = 0.8,
  edge_cutoff = cutoff,
  return_graph = TRUE
)

kME_conv_T_cytokine <- modules[, c("gene_name", "kME_conv_T_cytokine_signalling")]

conv_T_cytokine_plot <- make_network_plot(conv_T_cytokine_g, 
                                          kME_conv_T_cytokine, 
                                          "kME_conv_T_cytokine_signalling", 
                                          "#4682B4",
                                          kME_min = min, 
                                          kME_max = max)


# hubgene network plot for all T cytokine signalling
all_T_cytokine_g <- CustomNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=0,
  edge_prop = 1,
  mods = "all_T_cytokine_signalling",
  vertex.label.cex = 0.8,
  edge_cutoff = cutoff,
  return_graph = TRUE
)

kME_all_T_cytokine <- modules[, c("gene_name", "kME_all_T_cytokine_signalling")]

all_T_cytokine_plot <- make_network_plot(all_T_cytokine_g, 
                                         kME_all_T_cytokine, 
                                         "kME_all_T_cytokine_signalling", 
                                         "pink",
                                         kME_min = min, 
                                         kME_max = max)

p1 <- TCR_plot + 
  ggtitle("TCR module (top 10 hub genes)") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- T_activation_plot + 
  ggtitle("T cell activation module (top 10 hub genes)") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- conv_T_cytokine_plot + 
  ggtitle("conventional T cytokine signalling module (top 10 hub genes)") +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- all_T_cytokine_plot + 
  ggtitle("all T cytokine signalling module (top 10 hub genes)") +
  theme(plot.title = element_text(hjust = 0.5))

(p1 + p2) / (p3 + p4)
p1
p2
p3
p4



CustomNetworkPlot <- function (seurat_obj, mods = "all", n_hubs = 6, n_other = 3, 
                                sample_edges = TRUE, edge_prop = 0.5, return_graph = FALSE, 
                                edge.alpha = 0.25, vertex.label.cex = 0.5, hub.vertex.size = 4, 
                                other.vertex.size = 1, wgcna_name = NULL, edge_cutoff = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  if (all("all" %in% mods)) {
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
  } else {
    if (!all(mods %in% unique(as.character(modules$module)))) {
      stop(paste0("Some selected modules are not found in wgcna_name: ", 
                  wgcna_name))
    }
    modules <- modules %>% subset(module %in% mods)
  }
  TOM <- GetTOM(seurat_obj, wgcna_name)
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(n_hubs) %>% 
      .$gene_name
  })
  names(hub_list) <- mods
  other_genes <- modules %>% subset(!(gene_name %in% unlist(hub_list))) %>% 
    group_by(module) %>% sample_n(n_other, replace = TRUE) %>% 
    .$gene_name %>% unique
  selected_genes <- c(unlist(hub_list), other_genes)
  selected_modules <- modules %>% subset(gene_name %in% selected_genes)
  subset_TOM <- TOM[selected_genes, selected_genes]
  selected_modules$geneset <- ifelse(selected_modules$gene_name %in% 
                                       other_genes, "other", "hub")
  selected_modules$size <- ifelse(selected_modules$geneset == 
                                    "hub", hub.vertex.size, other.vertex.size)
  selected_modules$label <- ifelse(selected_modules$geneset == 
                                     "hub", as.character(selected_modules$gene_name), "")
  selected_modules$fontcolor <- ifelse(selected_modules$color == 
                                         "black", "gray50", "black")
  
  # --- Modified section: allow manual cutoff ---
  if (is.null(edge_cutoff)) {
    edge_cutoff <- min(sapply(1:nrow(subset_TOM), function(i) {
      max(subset_TOM[i, ])
    }))
    message("Computed edge cutoff value: ", signif(edge_cutoff, 4))
  } else {
    message("User-specified edge cutoff value: ", signif(edge_cutoff, 4))
  }
  # ---------------------------------------------
  
  edge_df <- reshape2::melt(subset_TOM) %>% subset(value >= edge_cutoff)
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), function(i) {
    gene1 = as.character(edge_df[i, "Var1"])
    gene2 = as.character(edge_df[i, "Var2"])
    col1 <- modules %>% subset(gene_name == gene1) %>% .$color
    col2 <- modules %>% subset(gene_name == gene2) %>% .$color
    if (col1 == col2) col1 else "grey90"
  })
  groups <- unique(edge_df$color)
  print(groups)
  if (sample_edges) {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * edge_prop))
      cur_df[cur_sample, ]
    }))
  } else {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), wt = value)
    }))
  }
  edge_df <- temp
  edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
  edge_df$color <- sapply(1:nrow(edge_df), function(i) {
    a = edge_df$value[i]
    alpha(edge_df$color[i], alpha = a)
  })
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, 
                                     vertices = selected_modules)
  l <- igraph::layout_with_fr(g, ...)
  if (return_graph) {
    return(g)
  }
  plot(g, layout = l, edge.color = adjustcolor(igraph::E(g)$color, 
                                               alpha.f = edge.alpha), vertex.size = igraph::V(g)$size, 
       edge.curved = 0, edge.width = 0.5, vertex.color = igraph::V(g)$color, 
       vertex.frame.color = igraph::V(g)$color, vertex.label = igraph::V(g)$label, 
       vertex.label.family = "Helvetica", vertex.label.font = 3, 
       vertex.label.color = igraph::V(g)$fontcolor, vertex.label.cex = vertex.label.cex, 
       ...)
}

make_network_plot <- function(g, kMEs, kME_name, colour, kME_min = NULL, kME_max = NULL) {
  # Add kME to vertex attributes
  V(g)$kME <- kMEs[[kME_name]][match(V(g)$name, kMEs$gene_name)]
  
  # Determine range for scaling
  if (is.null(kME_min)) kME_min <- min(V(g)$kME, na.rm = TRUE)
  if (is.null(kME_max)) kME_max <- max(V(g)$kME, na.rm = TRUE)
  
  # Scale node sizes consistently using provided or calculated limits
  V(g)$size <- scales::rescale(V(g)$kME, from = c(kME_min, kME_max), to = c(1, 25))
  
  # Using ggplot2 via ggraph to plot the igraph object
  ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.2, colour = colour) +
    geom_node_point(aes(size = kME), colour = colour) +
    geom_node_text(aes(label = name), repel = TRUE, size = 9) +
    scale_size(range = c(3, 12)) +
    theme_void()
}
