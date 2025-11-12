#in this script, we are visualising the expression of t activation hub genes across
#all our modules of interest

library(Seurat)
library(hdWGCNA)
library(patchwork)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

g2gato_obj$dataset <- "Sumanaweera_3D"

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_obj$dataset <- "Jing_2D"

#next ishi3d
ishi3d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_obj$dataset <- "Ishiguro_3D"

#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_obj$dataset <- "Ishiguro_2D"

comb_obj <- merge(g2gato_obj, c(EZ_obj, ishi3d_obj, ishi2d_obj))

comb_obj <- JoinLayers(comb_obj)

Idents(comb_obj) <- comb_obj$dataset

# Set the desired order
Idents(comb_obj) <- factor(Idents(comb_obj),
                           levels = c("Jing_2D", "Sumanaweera_3D", "Ishiguro_2D", "Ishiguro_3D"))

comb_obj <- NormalizeData(comb_obj)
comb_obj <- ScaleData(comb_obj)





gene_of_interest <- "ACTG1"
VlnPlot(comb_obj, features = gene_of_interest)

ishi3d_obj[["ACTG1"]] <- 0
ishi2d_obj[["ACTG1"]] <- 0

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")+
  scale_colour_gradientn(colors = c("grey"))  # force all grey

ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")+
  scale_colour_gradientn(colors = c("grey"))  # force all grey

EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")

g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")


actg1_plot <- wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "ACTG1 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )


gene_of_interest <- "CNN2"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

cnn2_plot <- wrap_plots(list(EZ_T_plot, g2gato_plot, ishi3d_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "CNN2 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

wrap_plots(list(actg1_plot, cnn2_plot), nrow = 1)

gene_of_interest <- "DENND2D"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "DENND2D expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )


gene_of_interest <- "ACTB"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi3d_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "DENND2D expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#also visualising T cell activation markers
gene_of_interest <- "IL2RA"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

g2gato_obj[["IL2RA"]] <- 0

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")+
  scale_colour_gradientn(colors = c("grey"))  # force all grey

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#also visualising T cell activation markers
gene_of_interest <- "CD28"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()


ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#also visualising T cell activation markers
gene_of_interest <- "CD38"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

EZ_obj[["CD38"]] <- 0
ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")+
  scale_colour_gradientn(colors = c("grey"))  # force all grey
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "PTPRC"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "CCR7"

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "TFRC"

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "CD69"

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = paste0(gene_of_interest," expression across iPSC datasets"),
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )
