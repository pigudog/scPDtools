Sys.setenv(GITHUB_PAT = "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")
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("zhanghao-njmu/SCP")
devtools::install_local("./zhanghao-njmu-SCP-v0.5.1-0-g458ce17.tar.gz")

library(Seurat)
# library(SeuratData)
library(SeuratDisk)
# InstallData("pbmc3k")
# data("pbmc3k.final")



# #从Seurat转成h5ad
# SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat")
# Convert("pbmc3k.h5Seurat", dest = "h5ad")
#从h5ad转成Seurat
Convert("data/mac_sparse.h5ad", dest = "h5seurat", overwrite = TRUE,assay = "RNA")
seurat_obj <- LoadH5Seurat("data/mac_sparse.h5seurat",
                           assays ="RNA"
                           # ,slots='counts'
                           )


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

devtools::install_github("cellgeni/sceasy")

library(sceasy)
library(reticulate)
# use_condaenv('scvi-env')
loompy <- reticulate::import('loompy')

# AnnData to Seurat
h5ad_file = './data/mac_manual_annotation_beautify.h5ad' 
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='mac_scp.rds')

mac_scp <- readRDS("D:/program/songwriter_BI/scbest/step6_SCP/data/mac_scp.rds")
save(mac_scp,file = "mac_scp.rda")

#
devtools::install_github('JiekaiLab/dior')
library(dior)
adata = read_h5(file='data/mac_sparse.h5', target.object = 'seurat')


## manual
# 使用reticulate
library(reticulate)
# 导入scanpy
sc <- import('scanpy')

# 读入h5ad文件
scvi_h5ad <- sc$read_h5ad('data/mac_sparse.h5ad')
scvi_h5ad

## 在这里我只取需要的信息
## 1.样本信息
meta = scvi_h5ad$obs # obs储存细胞信息，对应seuratobj@meta.data
var = scvi_h5ad$var  # var储存基因信息，这是由于anndata提取的矩阵没有行列名，所以这里是为了提取基因名

## 2.scvi坐标
X_umap_scVI <- scvi_h5ad$obsm['X_umap']
rownames(X_umap_scVI) = rownames(meta)
colnames(X_umap_scVI) = c('scVI_UMAP-1', 'scVI_UMAP-2')
head(X_umap_scVI)

X_pca <- scvi_h5ad$obsm['X_pca']
rownames(X_pca) = rownames(meta)
colnames(X_pca) = paste0("PCA-",c(1:dim(X_pca)[2]))
head(X_pca)

X_scVI <- scvi_h5ad$obsm['X_scVI']
rownames(X_scVI) = rownames(meta)
colnames(X_scVI) = paste0("scVI-",c(1:dim(X_scVI)[2]))
head(X_scVI)

# 提取矩阵，这里需要行列转置，因为anndata矩阵是 cells x genes，而seurat是genes x cells
# 注意保存的是稀疏矩阵 dgRMatrix 格式，到时候要转化成matrix, library("Matrix")
filter_counts = t(as.matrix(scvi_h5ad$layers['counts']))
dim(filter_counts)

# 提取的矩阵缺少行列名
colnames(filter_counts) = rownames(meta)
# 这里的矩阵是只保留了高可变基因
# rownames(filter_counts) = rownames(var)[var$highly_variable_features]
rownames(filter_counts)= rownames(var)
filter_counts[1:5,1:3]

# 如果你希望拿到原始矩阵，你可以从这里提取: adata$raw$X (前提是你在Scanpy中保存了原始矩阵)
# 对应的，你应该从adata$raw$var提取基因名
# 如果你对anndata对象不熟悉，可以忽略这几行注释
library(Seurat)
# rownames(filter_counts) <- gsub(".", "-", rownames(filter_counts))
rownames(filter_counts) <- gsub("_", "-", rownames(filter_counts))
colnames(meta)
meta = meta[,c("SampleID","DonorID","BiopsyType","Location","Binary Stage","Stage","cell type")]
colnames(meta) = c("SampleID","DonorID","BiopsyType","Location","Binary_Stage","Stage","cell_type")
# 重新构建seurat对象
identical(colnames(filter_counts),rownames(meta))



# counts 	Either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object
# project 	 Project name for the Seurat object
# assay 	Name of the initial assay
# names.delim 	For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE-CLUSTER-CELLTYPE, set this to “-” to separate the cell name into its component parts for picking the relevant field.
scobj <- CreateSeuratObject(counts = filter_counts, 
                            project = "Mac",
                            assay = "RNA",
                            names.delim = "-",
                            meta.data = meta)
scobj

# 主义python中保存的系数矩阵是class "dgRMatrix"
# 而R中为 dgCMatrix
class(scobj@assays$RNA@data)
log1p = t(as.matrix(scvi_h5ad$X))
colnames(log1p) = rownames(meta)
rownames(log1p)= rownames(var)
scobj@assays$RNA@data = as(log1p,"dgCMatrix")

# 添加 reductions
DefaultAssay(scobj)="RNA"
scobj@reductions[['X_pca']] <- CreateDimReducObject(embeddings = X_pca,assay = "RNA",key = "PCA_")
scobj[['X_scVI']] <- CreateDimReducObject(embeddings = X_scVI,assay = "RNA",key = "SCVI_")
scobj[['X_umap_scVI']] <- CreateDimReducObject(embeddings = X_umap_scVI,assay = "RNA",key = "UMAP_")

# 添加 feature.var
scobj@assays$RNA@meta.features = var
mac_scPDtools = scobj
save(mac_scPDtools,file="data/test.rda")
