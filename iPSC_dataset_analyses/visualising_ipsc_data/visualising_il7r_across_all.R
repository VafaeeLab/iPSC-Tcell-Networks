#in this script, we are visualising the expression of il7r across
#all our modules of interest

library(Seurat)
library(hdWGCNA)
library(patchwork)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

g2gato_obj$dataset <- "G2G_ATO"

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_obj$dataset <- "EZ_T"

#next ishi3d
ishi3d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_obj$dataset <- "ishi_3d"

#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_obj$dataset <- "ishi_2d"

comb_obj <- merge(g2gato_obj, c(EZ_obj, ishi3d_obj, ishi2d_obj))

comb_obj <- JoinLayers(comb_obj)

Idents(comb_obj) <- comb_obj$dataset

comb_obj <- NormalizeData(comb_obj)
comb_obj <- ScaleData(comb_obj)

VlnPlot(comb_obj, features = "IL7R", sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = "IL7R") + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = "IL7R") + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = "IL7R", reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = "IL7R", reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot,g2gato_plot,ishi2d_plot, ishi3d_plot ), ncol = 2) +
  plot_annotation(title = "IL7R expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )
EZ_obj[["TNFAIP3"]] <- 0
ishi3d_plot <- FeaturePlot(ishi3d_obj, features = "TNFAIP3") + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = "TNFAIP3") + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = "TNFAIP3", reduction = "harmony.umap") + ggtitle("Jing_2D") +
  scale_colour_gradientn(colors = c("grey"))  # force all grey
g2gato_plot <- FeaturePlot(g2gato_obj, features = "TNFAIP3", reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")



wrap_plots(list(EZ_T_plot,g2gato_plot,ishi2d_plot, ishi3d_plot ), ncol = 2) +
  plot_annotation(title = "TNFAIP3 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )
