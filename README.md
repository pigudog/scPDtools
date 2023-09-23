# DOWNLOAD
```R
devtools::install_github("pigudog/scPDtools")
```

# 1. Visualization
我们以windows运行为例
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

我们在内部内置了名为`ov_pallette`的 pallette
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

这个显示成分的函数我很喜欢
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
![](Pasted%20image%2020230923145416.png)

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

我比较喜欢这个`GroupHeatmap()`函数
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
![](Pasted%20image%2020230923145700.png)

# 2. differntial expression analysis
## RunDEtest
`RunDEtest`
- 改代码位于`R/SCP-analysis`

我们首先要清楚`lognormalize`的概念：
- `LogNormalize`: Feature counts for each cell are divided by the total counts for that cell and multiplied by the `scale.factor` (默认是10000). This is then natural-log transformed using `log1p`
- 也就是说seurat中的常规归一化使用的也是移位对数化，而scanpy中的处理吧`scale.factor`设置为$50\times10^4$, 这也是为什么在`sc.pp.normalize`之后，count的最大值常常会增大。同时也显示了我们可以直接使用anndata数据进行我们的differential gene的查找
	- 但是不同的不同的值会使得过度离散值$\alpha$发生改变，$\alpha$描述了数据集中存在比期望更大的变异性，这也很好理解，比如两个值分别是1和3，假设一个cell的`total_umi`=10000，假设log2，那么变化之后1和2或者5.7和7.2， L越大，基因的差距变小，也可以结合$log(1+x)$函数理解
	- 移位对数是一种快速归一化技术，优于其他揭示数据集潜在结构的方法（特别是在进行主成分分析时），并且有利于方差的稳定性，以进行后续的降维和差异表达基因的识别
- $f(y)=log(\frac{y}{s}+y_0)$
	- 其中y是原始的计数，s是尺寸因子，$y_0$是伪计数，细胞的尺寸因子=>$s_c=\frac{\sum_g{y_{gc}}}{L}$
	- L等用于`scale.factor`

其次我们要搞清楚differential gene 查找的方法的区别和分别使用方式
- 单细胞数据分析在进行完细胞自聚类或者细胞类型注释后，一般需要对查到的差异基因可视化，用来显示基因和细胞群的相关性，进行后续分析。当然Seurat和scanpy本身可视化的方式有非常多，例如feature plot, violin plot, dot plot等，但是问题在于差异基因分析后，**如何快速将每个细胞簇所对应的top deg汇总**，然后再对接函数绘制成图像。  
- Seurat的操作比较简单，因为`FindMarker()`后自身生成的就是一个数据框，但scanpy的`sc.tl.rank_genes_groups()`就没有那么用户友好了

  我们先来看一下两者实现方式的区别
  ```R
library(Seurat)
library(ggplot2)
library(dplyr)

deg<-FindAllMarker(data) #首先差异基因分析获取每个细胞簇的deg
top5 <- deg %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 5) #提取top差异基因，这里n=5代表top5

mark <- unique(top5$gene) #当然，这里可以是自己选的markers，来自背景知识的细胞标志物
p=DotPlot(kc,features = marker)
p

# 上面已经绘制完成了，下面这步只是纯粹地气泡图美化（可选）
p+ggtitle('there is the title')+theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10),#x轴标识
                                      axis.text.y = element_text(size = 10),#y轴标识
                                      legend.text = element_text(size= 10),legend.title= element_text(size= 10),#设置legend
                                      plot.title = element_text(hjust = 0.5,size = 12))+#设置标题居中
                                scale_colour_gradientn(colours = viridis::viridis(20))#修改成为CNS配色
```

```python
import scanpy as sc
import pandas as pd
import numpy as np

# find all degs
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
celltype=adata.obs['leiden'].unique().tolist() #把所有细胞簇种类拿出来
deg=sc.get.rank_genes_groups_df(adata,group=celltype) #把所有细胞簇对应的deg拿出来
deg.to_csv('./spGCN_deg.CSV') #存储备份

top=deg.groupby('group') 
top5=[] #同样以top5举例
for i in range(len(celltype)): #分群提取top5
    tmp=top.get_group(str(i))
    tmp=tmp.sort_values('scores',ascending=False) #按scores排序
    top5.append(tmp['names'].head(5).tolist())
#array list 转为 list
top5=np.array(top5)
top5=np.reshape(top5,5*len(celltype),'C').tolist() 

# visualization
sc.pl.dotplot(adata, top5, groupby='leiden')
plt.savefig('./top5_dotplot.png')
```


在scp中
```r
# 2. findmarkers
## Seurat的findallmarkers
markers = FindAllMarkers(mac_scp)
markers$group1 = markers$cluster
mac_scp@tools$DEtest_CellType$AllMarkers_seurat = markers

## findmarker一个个去取
# 理解一下findmarkers,findallmarkes,FindConservedMarkers
# 注意ident.2一般为用于比较的组，none就是其他剩下的
# b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
# nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
library(stats)
mac_scp <- RunDEtest(srt = mac_scp,
                          group_by = "CellType",
                          fc.threshold = 1,
                          only.pos = FALSE,
                          BPPARAM = BiocParallel::bpparam())
```

## VolcanoPlot
`VolcanoPlot`
- 在`R/SCP-plot.R`


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

![](Pasted%20image%2020230922144615.png)

![](Pasted%20image%2020230922211814.png)

两次尝试中，都发现findallmarker得到的gene更少，所以建议使用findmarker

## annoation features
```r
DEGs <- mac_scp@tools$DEtest_CellType$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
# 这似乎是直接从网上获取信息
# "TF",
mac_scp <- AnnotateFeatures(mac_scp,
                                 species = "Homo_sapiens", # Homo_sapiens,Mus_musculus
                                 db = c("TF","CSPA"))
table(mac_scp@assays[["RNA"]]@meta.features[["CSPA"]])
table(mac_scp@assays[["RNA"]]@meta.features[["TF"]])
table(mac_scp@assays[["RNA"]]@meta.features[["highly_variable_genes"]])
# featureheatmap这个函数最好在你进行enrichment之前使用，否则可能会报错，以及db参数会自动勋章并富集
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
![](Pasted%20image%2020230923160111.png)

# 3. Enrichment 
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
![](Pasted%20image%2020230923160554.png)

![](Pasted%20image%2020230923160625.png)

# 4. Dynamic heatmap
```r
# slingshot
mac_scp <- RunSlingshot(srt = mac_scp, group.by = "CellType", reduction = "UMAP")
# FeatureDimPlot(mac_scp, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
FeatureDimPlot(mac_scp, features = c("Lineage1","palantir_pseudotime","palantir_entropy"), 
               reduction = "UMAP", theme_use = "theme_blank",pt.size = 3)
```
![](Pasted%20image%2020230923161021.png)

```r
CellDimPlot(mac_scp, group.by = "CellType", reduction = "UMAP", 
            lineages = c("Lineage1","palantir_pseudotime","palantir_entropy"), 
            lineages_span = 1)
```
![](Pasted%20image%2020230923161150.png)


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
![](Pasted%20image%2020230923163303.png)
