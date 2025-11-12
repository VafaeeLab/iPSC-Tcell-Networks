#in this script, we are visualising the expression of all cyt signalling genes across
#all our modules of interest

library(Seurat)
library(hdWGCNA)
library(patchwork)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated.rds")

g2gato_obj$dataset <- "Sumanaweera_3D"
g2gato_obj$surface <- "3D"

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

EZ_obj$dataset <- "Jing_2D"
EZ_obj$surface <- "2D"

#next ishi3d
ishi3d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

ishi3d_obj$dataset <- "Ishiguro_3D"
ishi3d_obj$surface <- "3D"

#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ishi2d_obj$dataset <- "Ishiguro_2D"
ishi2d_obj$surface <- "2D"


comb_obj <- merge(g2gato_obj, c(EZ_obj, ishi3d_obj, ishi2d_obj))

comb_obj <- JoinLayers(comb_obj)

Idents(comb_obj) <- comb_obj$dataset

comb_obj <- NormalizeData(comb_obj)
comb_obj <- ScaleData(comb_obj)

# Set the desired order
Idents(comb_obj) <- factor(Idents(comb_obj),
                             levels = c("Jing_2D", "Sumanaweera_3D", "Ishiguro_2D", "Ishiguro_3D"))

gene_of_interest <- "TNFAIP3"

VlnPlot(comb_obj, 
        features = gene_of_interest, 
        group.by = c("dataset"), 
        sort = TRUE)
VlnPlot(comb_obj, 
        features = c("IL7R"), 
        sort = FALSE) +
  VlnPlot(comb_obj, 
          features = c("TNFAIP3"), 
          sort = FALSE) &
  NoLegend()



ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(ishi3d_plot, g2gato_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "TNFAIP3 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#visualising other hub genes in the all t cytokine signalling module
gene_of_interest <- "ZNF331"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(ishi3d_plot, g2gato_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "ZNF331 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#visualising other hub genes in the all t cytokine signalling module
gene_of_interest <- "SLC2A3"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(ishi3d_plot, g2gato_plot, ishi2d_plot, EZ_T_plot), ncol = 2) +
  plot_annotation(title = "SLC2A3 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#visualising other hub genes in the all t cytokine signalling module
gene_of_interest <- "ICOS"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(ishi3d_plot, g2gato_plot, ishi2d_plot, EZ_T_plot), ncol = 2) +
  plot_annotation(title = "SLC2A3 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#visualising other hub genes in the all t cytokine signalling module
gene_of_interest <- "CXCR4"
VlnPlot(comb_obj, features = gene_of_interest, sort = TRUE, split.by ="")

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishi_3d")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishi_2d")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("EZ_T")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("G2G ATO")

wrap_plots(list(ishi3d_plot, g2gato_plot, ishi2d_plot, EZ_T_plot), ncol = 2) +
  plot_annotation(title = "SLC2A3 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )