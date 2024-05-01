library(scPDtools)
library(Seurat)
# AnnData to Seurat
h5ad_file = './data/combined_immune.h5ad'
scPDtools::convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
                      outFile='combined_immune.rds')
combined_immune <- readRDS("D:/program/scPDtools/combined_immune.rds")
save(combined_immune,file = 'data/combined_immune.rda')
library(SingleR)
# 人用下面, [4.2]以上R需要用到celldex::HumanPrimaryCellAtlasData()!!!!
refdata <- SingleR::HumanPrimaryCellAtlasData()
# refdata <- SingleR::DatabaseImmuneCellExpressionData()

# 鼠用下面
#refdata <- SingleR::MouseRNAseqData()
scRNA = combined_immune
library(Seurat)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$leiden
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)


scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$leiden_immune == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}


Idents(scRNA)=scRNA$celltype

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

DimPlot(scRNA,group.by = 'major_celltype',split.by = 'disease_type',cols = col_vector[20:30])
DimPlot(scRNA,group.by = 'immune_celltype',split.by = 'disease_type',cols = col_vector[20:30])

scRNA@meta.data$celltype = scRNA@meta.data$major_celltype
Idents(scRNA)=scRNA$celltype
save(scRNA,file = 'twodisease_immune.rda')

## 细胞比例图
library(reshape2)
library(ggplot2)
library(dplyr)
prop_df <- table(scRNA@meta.data$celltype,scRNA@meta.data$disease_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_color <- col_vector[1:10]

# 作图
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
prop## 细胞比例图
