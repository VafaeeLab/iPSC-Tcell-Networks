#in this script, we are visualising the expression of cd8a and b genes
#in our metacells from reference, ishi2d, ez t and g2g ATO cells

library(Seurat)
library(hdWGCNA)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated_NO_TNF.rds")
g2gato_metacells <- GetMetacellObject(g2gato_obj)

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")
EZ_metacells <- GetMetacellObject(EZ_obj)

#then ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")
ishi2d_metacells <- GetMetacellObject(ishi2d_obj)

ref_obj <- readRDS(file = "~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/annotated_reference_obj.rds")
ref_metacells <- GetMetacellObject(ref_obj)

#checking phenotypes of metacells
unique(g2gato_metacells$maturity)
unique(EZ_metacells$phenotype)
unique(ishi2d_metacells$orig.ident) #only one phenotype
unique(ref_metacells$celltype_annotation)

#subsetting EZ and g2gato to mature T phenotypes
g2gato_metacells <- subset(g2gato_metacells, subset = maturity == "mature_T")
EZ_metacells <- subset(EZ_metacells, subset = phenotype == "mature_T")

#subsetting ref obj to just cd8+T cells
ref_metacells <- subset(ref_metacells, subset = celltype_annotation == "CD8+T")
ref_metacells <- NormalizeData(ref_metacells)

FeatureScatter(ref_metacells, feature1 = "CD8A", feature2 = "CD8B", slot = "counts", group.by = "celltype_annotation") +
  NoLegend() #+ ggtitle(label = "reference cd8+ T metacells")
FeatureScatter(EZ_metacells, feature1 = "CD8A", feature2 = "CD8B", slot = "data") + 
  NoLegend() #+ ggtitle(label = "CD8A vs CD8B expression in Jing-2D mature T metacells")
ishi2d_plot <- FeatureScatter(ishi2d_metacells, feature1 = "CD8A", feature2 = "CD8B", slot = "data") +
  NoLegend() + ggtitle(label = "Ishiguro-2D")
g2g_plot <- FeatureScatter(g2gato_metacells, feature1 = "CD8A", feature2 = "CD8B", slot = "data") +
  NoLegend() + ggtitle(label = "Sumanaweera-3D")

ishi2d_plot + g2g_plot

#making the same plot but with real cells
#subsetting EZ and g2gato to mature T phenotypes
g2gato_obj <- subset(g2gato_obj, subset = maturity == "mature_T")
EZ_obj <- subset(EZ_obj, subset = phenotype == "mature_T")

#subsetting ref obj to just cd8+T cells
ref_obj <- subset(ref_obj, subset = celltype_annotation == "CD8+T")

FeatureScatter(ref_obj, feature1 = "CD8A", feature2 = "CD8B", slot = "counts", group.by = "celltype_annotation") + ggtitle(label = "reference cd8+ T cells")
FeatureScatter(EZ_obj, feature1 = "CD8A", feature2 = "CD8B", slot = "data", group.by = "phenotype", split.by = "treatment") + ggtitle(label = "EZ T mature T cells")
FeatureScatter(ishi2d_obj, feature1 = "CD8A", feature2 = "CD8B", slot = "data") + ggtitle(label = "Ishi2d mature T cells")
FeatureScatter(g2gato_obj, feature1 = "CD8A", feature2 = "CD8B", slot = "data", group.by = "maturity") + ggtitle(label = "G2G ATO mature T cells")
