if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

devtools::install_github("cellgeni/sceasy")

# library(sceasy)
library(scPDtools)
library(reticulate)
# use_condaenv('scvi-env')
# loompy <- reticulate::import('loompy')

# AnnData to Seurat
h5ad_file = './data/mac_sparse.h5ad'
anndata2seurat(h5ad_file,outFile='data/mac_sparse.rds')
adata = convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
                      outFile='data/mac_sparse.rds')
adata=annd
mac_scp = mac_sparse
save(mac_scp,file="./data/mac_menstrual.rda")

adata@meta.data$Phase = adata@meta.data$`Binary Stage`

CellDimPlot(
  srt = adata, group.by = "cell_type",
  pt.size = 5,
  stat.by = "Phase",
  reduction = "UMAP",
  theme_use = "theme_blank"
)
