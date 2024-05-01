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
convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
                      outFile='data/data_M_annotation.rds')


library(Seurat)
adata = readRDS("./data/data_M_annotation.rds")

save(adata,file="./data/data_M_annotation.rda")

ov_palette =c('#A499CC',
              '#E0A7C8',
              '#E069A6',
              "#f1707d",
              "#AFC2D9",
              "#6894B9",
              "#79B99D",
              "#F5D2A8",
              "#D2EBC8")

p = CellRatioPlot(object = scRNA,
                  sample.name = "Type",
                  celltype.name = "celltype",
                  fill.col = ov_palette)
ggsave(p,filename=paste0("preanalysis",'/ratioplot.pdf'),width = 5,height = 6)
