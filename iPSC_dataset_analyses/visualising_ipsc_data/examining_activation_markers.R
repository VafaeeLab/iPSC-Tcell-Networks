#here we examine activation markers

library(Seurat)
library(hdWGCNA)
library(ggplot2)

seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")

ModuleFeaturePlot(seurat_obj, reduction = "harmony.umap", module_names = "T_cell_activation")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "IL2RA")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CD28")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "PECAM1")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CCR7")


FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CCL5")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "IL7R")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "SOX4")



seurat_obj <- ModuleConnectivity(seurat_obj)

modules <- GetModules(seurat_obj)

seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated_NO_TNF.rds")

seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "cell_line")

ModuleFeaturePlot(seurat_obj, reduction = "harmony.umap", module_names = "T_cell_activation")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "IL2RA")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CD28")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "PECAM1")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CCR7")


FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "CCL5")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "IL7R")

FeaturePlot(seurat_obj, reduction = "harmony.umap", features = "SOX4")

seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

ModuleFeaturePlot(seurat_obj, reduction = "umap", module_names = "T_cell_activation")

FeaturePlot(seurat_obj, reduction = "umap", features = "IL2RA")

FeaturePlot(seurat_obj, reduction = "umap", features = "CD28")

FeaturePlot(seurat_obj, reduction = "umap", features = "PECAM1")

FeaturePlot(seurat_obj, reduction = "umap", features = "CCR7")


FeaturePlot(seurat_obj, reduction = "umap", features = "CCL5")

FeaturePlot(seurat_obj, reduction = "umap", features = "IL7R")

FeaturePlot(seurat_obj, reduction = "umap", features = "SOX4")

FeaturePlot(seurat_obj, reduction = "umap", features = "PTPRC")

FeaturePlot(seurat_obj, reduction = "umap", features = "CD2")

seurat_obj <- ModuleConnectivity(seurat_obj)

modules <- GetModules(seurat_obj)

seurat_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro_ato/w6/processed/annotated_obj.rds")

FeaturePlot(seurat_obj, reduction = "umap", features = "IL2RA")

FeaturePlot(seurat_obj, reduction = "umap", features = "CD28")

FeaturePlot(seurat_obj, reduction = "umap", features = "PECAM1")

FeaturePlot(seurat_obj, reduction = "umap", features = "CCR7")


FeaturePlot(seurat_obj, reduction = "umap", features = "CCL5")

FeaturePlot(seurat_obj, reduction = "umap", features = "IL7R")

FeaturePlot(seurat_obj, reduction = "umap", features = "SOX4")

FeaturePlot(seurat_obj, reduction = "umap", features = "PTPRC")

FeaturePlot(seurat_obj, reduction = "umap", features = "CD2")

ModuleFeaturePlot(seurat_obj, module_names = c("all_T_cytokine_signalling", "conv_T_cytokine_signalling"))
