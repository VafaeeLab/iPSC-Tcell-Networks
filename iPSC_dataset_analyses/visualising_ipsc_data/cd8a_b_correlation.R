#in this script, we are visualising the expression of tcr hub genes across
#all our modules of interest

library(Seurat)
library(hdWGCNA)

#start with g2g ato
g2gato_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/g2g_ATO/ATO_integrated_NO_TNF.rds")


## Examining correlation of CD8a and b genes
g2gato_mtx <- GetDatExpr(g2gato_obj)

# 1. All pairwise gene–gene correlations (upper triangle only)
all_cors <- cor(g2gato_mtx, method = "pearson")
all_cors_vec <- all_cors[upper.tri(all_cors, diag = FALSE)]

# 2. Correlation between these two genes
g2gato_cor <- cor(g2gato_mtx[, "CD8A"], g2gato_mtx[, "CD8B"], method = "pearson")

# 3. Percentile of your observed correlation
g2gato_percentile <- ecdf(all_cors_vec)(g2gato_cor) * 100

#now EZ_T
EZ_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/EZ_T/processed/EZ_T.rds")


## Examining correlation of CD8a and b genes
EZ_mtx <- GetDatExpr(EZ_obj)

# 1. All pairwise gene–gene correlations (upper triangle only)
all_cors <- cor(EZ_mtx, method = "pearson")
all_cors_vec <- all_cors[upper.tri(all_cors, diag = FALSE)]

# 2. Correlation between these two genes
EZ_cor <- cor(EZ_mtx[, "CD8A"], EZ_mtx[, "CD8B"], method = "pearson")

# 3. Percentile of your observed correlation
EZ_percentile <- ecdf(all_cors_vec)(EZ_cor) * 100


#lastly ishi2d
ishi2d_obj <- readRDS(file = "/g/data/yr31/hb1618/iPSC_project/datasets/iPSC/ishiguro/2d/processed/annotated_obj.rds")

## Examining correlation of CD8a and b genes
ishi2d_mtx <- GetDatExpr(ishi2d_obj)

# 1. All pairwise gene–gene correlations (upper triangle only)
all_cors <- cor(ishi2d_mtx, method = "pearson")
all_cors_vec <- all_cors[upper.tri(all_cors, diag = FALSE)]

# 2. Correlation between these two genes
ishi2d_cor <- cor(ishi2d_mtx[, "CD8A"], ishi2d_mtx[, "CD8B"], method = "pearson")

# 3. Percentile of your observed correlation
ishi2d_percentile <- ecdf(all_cors_vec)(ishi2d_cor) * 100

#and reference
seurat_ref <- readRDS(file = "~/ipsc_project/R_scripts/iPSC_analysis/reference_data_analysis/annotated_reference_obj.rds")

#subset metacell mtx to CD8+ T cells
seurat_ref <- NormalizeMetacells(seurat_ref)
metacell_obj <- GetMetacellObject(seurat_ref)
metacell_obj <- subset(metacell_obj, subset = celltype_annotation == "CD8+T")
cd8_ref_mtx <- as.matrix(GetAssayData(metacell_obj))

# 1. All pairwise gene–gene correlations (upper triangle only)
all_cors <- cor(cd8_ref_mtx, method = "pearson")
all_cors_vec <- all_cors[upper.tri(all_cors, diag = FALSE)]

# 2. Correlation between these two genes
cd8_ref_cor <- cor(cd8_ref_mtx[, "CD8A"], cd8_ref_mtx[, "CD8B"], method = "pearson")

# 3. Percentile of your observed correlation
cd8_ref_percentile <- ecdf(all_cors_vec)(cd8_ref_cor) * 100