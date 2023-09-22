if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

devtools::install_github("cellgeni/sceasy")

library(sceasy)
library(reticulate)
# use_condaenv('scvi-env')
loompy <- reticulate::import('loompy')

# AnnData to Seurat
h5ad_file = './data/mac_palantir.h5ad' 
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='data/mac_scp.rds')

save(mac_scp,file="./data/mac_menstrual.rda")
