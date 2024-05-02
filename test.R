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

################################################################################
# 2. visualization
################################################################################
data("data_M_annotation")
library(dplyr)
ov_palette =c('#A499CC',
              '#E0A7C8',
              '#E069A6',
              "#f1707d",
              "#AFC2D9",
              "#6894B9",
              "#79B99D",
              "#F5D2A8",
              "#D2EBC8")

p = CellRatioPlot(object = adata,
                  sample.name = "disease_type",
                  celltype.name = "subset_celltype",flow.curve = 0.3,
                  fill.col = ov_palette)
p
ggsave(p,filename=paste0("preanalysis",'/ratioplot.pdf'),width = 5,height = 6)

adata <- RunUMAP(adata,reduction = 'scVI', dims = 1:30, verbose = FALSE)
p1 = DimPlot(adata, group.by="subset_celltype", label=F, label.size=4.5, reduction='umap',cols = ov_palette )
p1
p2 = DimPlot(adata, group.by="subset_celltype", label=F, label.size=4.5, reduction='mde_M',cols = ov_palette )
p2
ggsave(p,filename=paste0("preanalysis",'/umap_anno.pdf'),width = 7,height = 6)




