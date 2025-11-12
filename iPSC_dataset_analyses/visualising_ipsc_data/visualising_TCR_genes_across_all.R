#in this script, we are visualising the expression of tcr hub genes across
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



gene_of_interest <- "CD8A"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLege

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "CD8A expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )


gene_of_interest <- "CD8B"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "CD8B expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )


gene_of_interest <- "TRAC"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "TRAC expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "TRBC1"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "TRBC1 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "TRBC2"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi3d_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "TRBC2 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )


gene_of_interest <- "CTSW"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi3d_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "CTSW expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

#now genes associated with non-conventional TCR
gene_of_interest <- "TRDC"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")

EZ_T_plot + ggtitle("TRDC expression")

wrap_plots(list(EZ_T_plot, g2gato_plot, ishi3d_plot, ishi2d_plot), ncol = 2) +
  plot_annotation(title = "TRDC expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "TRGC1"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE)

ishi3d_plot <- FeaturePlot(ishi3d_obj, features = gene_of_interest) + ggtitle("Ishiguro_3D")
ishi2d_plot <- FeaturePlot(ishi2d_obj, features = gene_of_interest) + ggtitle("Ishiguro_2D")
EZ_T_plot <- FeaturePlot(EZ_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Jing_2D")
g2gato_plot <- FeaturePlot(g2gato_obj, features = gene_of_interest, reduction = "harmony.umap") + ggtitle("Sumanaweera_3D")


wrap_plots(list(EZ_T_plot, g2gato_plot, ishi2d_plot, ishi3d_plot), ncol = 2) +
  plot_annotation(title = "TRGC1 expression across iPSC datasets",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

gene_of_interest <- "TRGC2"
VlnPlot(comb_obj, features = gene_of_interest, sort = FALSE) + NoLegend()

FeatureScatter(EZ_obj, feature1 = "TRDC", feature2 = "TRGC1", group.by = "treatment")
DimPlot(EZ_obj, group.by = "treatment", reduction = "harmony.umap") +
  FeaturePlot(EZ_obj, features = "TRDC", reduction = "harmony.umap") +
  FeaturePlot(EZ_obj, features = "TRGC1", reduction = "harmony.umap")


FeaturePlot(EZ_obj, features = "TRGC1", reduction = "harmony.umap", split.by = "treatment")
FeaturePlot(EZ_obj, features = "TRDC", reduction = "harmony.umap", split.by = "treatment")

alpha_beta_genes <- c("TRAC", "TRBC1", "TRBC2")
gamma_delta_genes <- c("TRDC", "TRGC1", "TRGC2")
genes_of_interest <- c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2")
VlnPlot(comb_obj, features = genes_of_interest)
VlnPlot(comb_obj, features = alpha_beta_genes)
VlnPlot(comb_obj, features = gamma_delta_genes)
VlnPlot(comb_obj, features = c("TRAC", "TRDC"))
FeatureScatter(EZ_obj, feature1 = "TRAC", feature2 = "TRDC", group.by = "treatment")
