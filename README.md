# scPDtools: A Single-Cell Pipeline with R
scPDtools provides several comprehensive set of tools for visulization to assist python in single-cell analysis

The package includes the following facilities
- 
- High-quality data visualization methods.
- Differntial gene expression analysis
- Compositional analysis
- Gene set enrichment and pathway analysis
- Some R tools in inferring trajectories and Visulization

The functions in the SCP package are all developed around the [Seurat object](https://github.com/mojaveazure/seurat-object) and are compatible with other Seurat functions.
- If you only have data in the form of AnnData, you can refer to `./step0.R` for conversion

scPDtools refers to the following packages
- [SCP](https://github.com/zhanghao-njmu/SCP)
- [ClusterGVis](www.github.com/junjunlab/ClusterGVis)
- [sceasy](https://github.com/cellgeni/sceasy/)
## R version requirement
- R >= 4.2.0

## Installation in the global R environment
You can install the latest version of scPDtools from [GitHub](https://github.com/zhanghao-njmu/SCP) with:
```R
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("pigudog/scPDtools")
```

```r
################################################################################
# 1. Anndata -> Seurat (V4)
################################################################################
library(scPDtools)

# AnnData to Seurat
h5ad_file = './data/data_M_annotation.h5ad'
convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
                      outFile='data/data_M_annotation.rds')


library(Seurat)
adata = readRDS("./data/data_M_annotation.rds")

save(adata,file="./data/data_M_annotation.rda")

```

output:
```
The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
which was just loaded, will retire in October 2023.
Please refer to R-spatial evolution reports for details, especially
https://r-spatial.org/r/2023/05/15/evolution4.html.
It may be desirable to make the sf package available;
package maintainers should consider adding sf to Suggests:.
The sp package is now running under evolution status 2
     (status 2 uses the sf package in place of rgdal)
layers.`log1p` -> data; layers.`counts` -> counts
```
# Intro
# 1. Anndata -> Seurat (V4)


# Old Version
## 1. Visualization
Let's take windows running as an example
- We set the number of parallel processes to 8
```r
library(scPDtools)
library(Seurat)
library(BiocParallel)
param <- bpparam()
bpworkers(param) <- 8##设置并行的进程数
BiocParallel::bpparam()
# class: SnowParam
# bpisup: FALSE; bpnworkers: 8; bptasks: 0; bpjobname: BPJOB
# bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
# bpRNGseed: ; bptimeout: NA; bpprogressbar: FALSE
# bpexportglobals: TRUE; bpexportvariables: TRUE; bpforceGC: FALSE
# bpfallback: TRUE
# bplogdir: NA
# bpresultdir: NA
# cluster type: SOCK

# load data
data("mac_menstrual")
# mac_scp = NormalizeData(mac_scp)
mac_scp@meta.data$CellType = mac_scp@meta.data$`cell type`
mac_scp@meta.data$Binary_Stage <- mac_scp@meta.data$`Binary Stage`
Idents(mac_scp) = mac_scp@meta.data$CellType
```

We have built-in pallette called 'ov_pallette'
- To correspond to our pallette used in python
```r
######################################################
# color主要是使用RdBu和RdPu和ov_pallette
colnames(mac_scp@meta.data)
# 1.show
## normal
CellDimPlot(
  srt = mac_scp, group.by = c("CellType", "Stage"),
  pt.size = 3,palette = "ov_palette",
  reduction = "UMAP",
  theme_use = "theme_blank",
  subtitle = "celltype and stage"
)
```
![](README/Pasted%20image%2020230924114443.png)

We can set the 'stat.by' parameter to achieve the composition calculation, which can more intuitively see the cell composition in different periods
```R
## stage and celldimplot
CellDimPlot(
  srt = mac_scp, group.by = "CellType",
  pt.size = 5,palette = "ov_palette",
  stat.by = "Stage",
  reduction = "UMAP",
  theme_use = "theme_blank"
)
```
![](README/Pasted%20image%2020230923145416.png)

```r
## feature
FeatureDimPlot(
  srt = mac_scp, features = c("FOLR2", "SPARCL1", "SPP1", "FCN1"),
  reduction = "UMAP", theme_use = "theme_blank",
  palette = "RdPu",pt.size = 3

)
FeatureDimPlot(
  srt = mac_scp, features = c("FOLR2", "SPARCL1", "SPP1", "FCN1"),
  reduction = "UMAP", theme_use = "theme_blank",
  palette = "RdBu",pt.size = 3

)
```
![](README/Pasted%20image%2020230924120103.png)
```
FeatureDimPlot(
  srt = mac_scp, features = c("FOLR2", "SPARCL1", "SPP1", "FCN1"),
  compare_features = TRUE,
  label = TRUE,
  label_insitu = TRUE,
  reduction = "UMAP",
  theme_use = "theme_blank",
  pt.size = 5,
  palette = "ov_palette",
)
```
![](README/Pasted%20image%2020230924120244.png)

I prefer the 'GroupHeatmap()' function
```r
## groupheatmap!!
library(dplyr)
ht <- GroupHeatmap(
  srt = mac_scp,
  features = c(
    "FCN1", "EREG", # MONO
    "FOLR2", "SEPP1", # FOLR2
    "APOC1", "APOE","SPP1", # APOC1
    "SPARCL1", "CALD1"# SPARCL1

  ),
  group.by = c("CellType", "Stage"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Stage", "Binary_Stage"),
  cell_annotation_palette = c("Dark2", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE,
  libsize = 1
)
print(ht$plot)
```
![](README/Pasted%20image%2020230923145700.png)

## 2. differntial expression analysis
### RunDEtest
```r
# 2. findmarkers
## findallmarkers in Seurat
markers = FindAllMarkers(mac_scp)
markers$group1 = markers$cluster
mac_scp@tools$DEtest_CellType$AllMarkers_seurat = markers

## findmarker is used in RunDEtest
# you need to understand findmarkers,findallmarkes,FindConservedMarkers
# note: ident.2 means the group for compare
mac_scp <- RunDEtest(srt = mac_scp,
                          group_by = "CellType",
                          fc.threshold = 1,
                          only.pos = FALSE,
                          BPPARAM = BiocParallel::bpparam())
```

### VolcanoPlot
```r
# visualization
VolcanoPlot(srt = mac_scp, group_by = "CellType",
            test.use = "seurat",
            DE_threshold = "avg_log2FC > 0.2 & p_val_adj < 0.05",
            palette = "RdBu",
            palcolor = NULL)

VolcanoPlot(srt = mac_scp, group_by = "CellType",
                        test.use = "wilcox",
                        DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                        palette = "RdBu",
                        palcolor = NULL)
```
![](README/Pasted%20image%2020230924120754.png)
`findallmarker` seems to get fewer genes when trying, so it is recommended to use `findmarker` in `DEtest` or adjust the parameters in `findallmarkers`

### annoation features
```r
DEGs <- mac_scp@tools$DEtest_CellType$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
mac_scp <- AnnotateFeatures(mac_scp,
                                 species = "Homo_sapiens", # Homo_sapiens,Mus_musculus
                                 db = c("TF","CSPA"))
table(mac_scp@assays[["RNA"]]@meta.features[["CSPA"]])
table(mac_scp@assays[["RNA"]]@meta.features[["TF"]])
table(mac_scp@assays[["RNA"]]@meta.features[["highly_variable_genes"]])
ht <- FeatureHeatmap(
  srt = mac_scp, group.by = "CellType",
  features = DEGs$gene,
  feature_split = DEGs$group1,
  species = "Homo_sapiens",
  db = c("GO_BP", "KEGG", "WikiPathway"),
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"),c("forestgreen")),
  height = 5, width = 4,
  libsize = 1
)
print(ht$plot)
```
![](README/Pasted%20image%2020230923160111.png)

## 3. Enrichment 
```r
# enrichment
mac_scp <- RunEnrichment(
  srt = mac_scp,
  test.use = "wilcox",
  group_by = "CellType",
  db = c("MSigDB","GO_BP", "KEGG", "WikiPathway"), # c("GO_BP","KEGG", "WikiPathway", "Reactome", "PFAM", "MP","MSigDB", "MSigDB_MH")
  species = "Homo_sapiens",

  DE_threshold = "avg_log2FC > log2(0.5) & p_val_adj < 1"
)
table(mac_scp@meta.data$CellType)
EnrichmentPlot(
  srt = mac_scp,
  db = c("MSigDB"),
  group_by = "CellType",
  group_use = c("FOLR2+ Mac" ,"SPARCL1+ Mac" ),# ,   "SPP1+ Mac" 
  plot_type = "bar",
  topTerm = 10,
  padjustCutoff = 1

)
EnrichmentPlot(
  srt = mac_scp, group_by = "CellType", group_use = c("FOLR2+ Mac" ,"SPARCL1+ Mac" , "SPP1+ Mac" ),
  plot_type = "bar"
)
```
![](README/Pasted%20image%2020230923160554.png)

![](README/Pasted%20image%2020230923160625.png)

## 5. RunGSEA
```r
mac_scp <- RunGSEA(
  srt = mac_scp, group_by = "CellType", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
GSEAPlot(srt = mac_scp, group_by = "CellType", group_use = "SPP1+ Mac", id_use = "GO:0007186")
```

## 6. Dynamic heatmap
```r
# slingshot
mac_scp <- RunSlingshot(srt = mac_scp, group.by = "CellType", reduction = "UMAP")
# FeatureDimPlot(mac_scp, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
FeatureDimPlot(mac_scp, features = c("Lineage1","palantir_pseudotime","palantir_entropy"), 
               reduction = "UMAP", theme_use = "theme_blank",pt.size = 3)
```
![](README/Pasted%20image%2020230923161021.png)

```r
CellDimPlot(mac_scp, group.by = "CellType", reduction = "UMAP", 
            lineages = c("Lineage1","palantir_pseudotime","palantir_entropy"), 
            lineages_span = 1)
```
![](README/Pasted%20image%2020230923161150.png)


```r
mac_scp <- RunDynamicFeatures(srt = mac_scp, lineages = c("Lineage1","palantir_pseudotime"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = mac_scp,
  lineages = c("Lineage1","palantir_pseudotime"),
  use_fitted = TRUE, n_split = 6,
  reverse_ht = "palantir_pseudotime",
  species = "Homo_sapiens",
  db = "MSigDB",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  heatmap_palette = "RdPu",
  cell_annotation = "CellType",
  # separate_annotation = list("CellType", c("SPP1+ Mac", "FOLR2+ Mac" )),
  # separate_annotation_palette = c("Paired", "Set1"),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot)
```
![](README/Pasted%20image%2020230923163303.png)
