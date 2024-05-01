# The R package development
library(usethis)
library(devtools)
library(roxygen2)
# 检查
has_devel()


################################################################################
# 1. Anndata -> Seurat (V4)
################################################################################
library(scPDtools)
# library(reticulate)

# AnnData to Seurat
h5ad_file = './data/data_M_annotation.h5ad'
# anndata2seurat(h5ad_file,outFile='data_M_annotation.rds')
adata = convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
                      outFile='data/mac_sparse.rds')
adata=annd
mac_scp = mac_sparse
save(mac_scp,file="./data/mac_menstrual.rda")
