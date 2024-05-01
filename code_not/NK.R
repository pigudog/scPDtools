# library(sceasy)
library(scPDtools)
library(reticulate)
# use_condaenv('scvi-env')
# loompy <- reticulate::import('loompy')

# AnnData to Seurat
h5ad_file = './data/menstrual_NK.h5ad'
adata = convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p")
