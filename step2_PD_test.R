library("scPDtools")
library(Seurat)
library(BiocParallel)




##查看系统存在的并行环境
registered()
# 根据上面的信息我们可以看到在linux和mac中MulticoreParam是特有的。另外两个在所有平台都是存在的

##两个数据进行顺序并列运行
fun <- function(greet, who) {
  paste(Sys.getpid(), greet, who)## Sys.getpid()获取任务ID
}
greet <- c("morning", "night")
who <- c("sun", "moon")
param <- bpparam()
original <- bpworkers(param)
bpworkers(param) <- 4##设置并行的进程数
result <- bpmapply(fun, greet, who, BPPARAM = param)
# BPPARAM = BiocParallel::bpparam()
register(MulticoreParam(workers = 8, progressbar = TRUE))
BiocParallel::bpparam()

data("mac_menstrual")
data("mac_anno")
palette_scp(palette = "ov_palette")
ov_palette =c('#A499CC','#5E4D9A',
                      '#1F577B',
                      '#A56BA7',
                      '#E0A7C8',
                      '#E069A6',
                      '#941456',
                      '#FCBC10',
                      '#EAEFC5',
                      '#01A0A7',
                      '#75C8CC',
                      '#F0D7BC',
                      '#D5B26C',
                      '#D5DA48',
                      '#B6B812',
                      '#9DC3C3',
                      '#A89C92',
                      '#FEE00C',
                      '#FEF2A1',
                      '#7CBB5F',
                      '#368650',
                      '#279AD7',
                      '#78C2ED',
                      '#866017',
                      '#9F987F',
                      '#E0DFED',
                      '#EF7B77',
                      '#F0EEF0')
# palette_list[["ov_palette"]]=ov_palette
# save(palette_list,file="~/scPDtools/SCP-main/R/sysdata.rda")

# 理解一下findmarkers,findallmarkes,FindConservedMarkers
# 注意ident.2一般为用于比较的组，none就是其他剩下的
# b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
# nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)

# mac_scp = NormalizeData(mac_scp)
mac_scp@meta.data$CellType = mac_scp@meta.data$`cell type`
mac_scp@meta.data$Binary_Stage <- mac_scp@meta.data$`Binary Stage`
Idents(mac_scp) = mac_scp@meta.data$CellType

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

## stage and celldimplot
CellDimPlot(
  srt = mac_scp, group.by = "CellType",
  pt.size = 5,palette = "ov_palette",
  stat.by = "Stage",
  reduction = "UMAP",
  theme_use = "theme_blank"
)

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


# 2. findmarkers
# ## Seurat的findallmarkers
# markers = FindAllMarkers(mac_scp)
# markers$group1 = markers$cluster
# mac_scp@tools$DEtest_CellType$AllMarkers_seurat = markers

## findmarker一个个去取
mac_scp <- RunDEtest(srt = mac_scp,
                          group_by = "CellType",
                          fc.threshold = 1,
                          only.pos = FALSE,
                          BPPARAM = BiocParallel::bpparam())
# palettes
# source("R/utils.R")
#
# VolcanoPlot(srt = mac_scp, group_by = "CellType",
#             test.use = "seurat",
#             DE_threshold = "avg_log2FC > 0.2 & p_val_adj < 0.05",
#             palette = "RdBu",
#             palcolor = NULL)

VolcanoPlot(srt = mac_scp, group_by = "CellType",
                        test.use = "wilcox",
                        DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                        palette = "RdBu",
                        palcolor = NULL)



#
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

# SankeyPlot()

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
  group_use = c("FOLR2+ Mac" ,"SPARCL1+ Mac" ,   "SPP1+ Mac" ),
  plot_type = "bar",
  topTerm = 10,
  padjustCutoff = 1

)
EnrichmentPlot(
  srt = mac_scp, group_by = "CellType", group_use = c("FOLR2+ Mac" ,"SPARCL1+ Mac" , "SPP1+ Mac" ),
  plot_type = "bar"
)

# slingshot
mac_scp <- RunSlingshot(srt = mac_scp, group.by = "leiden_res1", reduction = "UMAP")
FeatureDimPlot(mac_scp, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(mac_scp, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)

mac_scp <- RunDynamicFeatures(srt = mac_scp, lineages = c("Lineage1", "Lineage2","Lineage3"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = mac_scp,
  lineages = c("Lineage1", "Lineage2"),
  use_fitted = TRUE, n_split = 6,
  reverse_ht = "Lineage1",
  species = "Mus_musculus",
  db = "MSigDB",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  heatmap_palette = "RdPu",
  cell_annotation = "SubCellType",
  separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
  separate_annotation_palette = c("Paired", "Set1"),
  # feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot)

