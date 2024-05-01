#' Title
#'
#' @param object seurat object
#' @param sample.name the group (treat vs ctrl)
#' @param celltype.name the celltype
#' @param col.width the width of col
#' @param flow.alpha the alpha of flow
#' @param flow.curve the curve of flow
#' @param fill.col the color you filled in
#'
#' @return p
#' @export
#'
#' @examples p = CellRatioPlot(object = scRNA,sample.name = "Type",celltype.name = "celltype",fill.col = ov_palette)  ggsave(p,filename=paste0("preanalysis",'/ratioplot.pdf'),width = 5,height = 6)
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
CellRatioPlot <- function(object = NULL,
                          sample.name = NULL,
                          celltype.name = NULL,
                          col.width = 0.7,
                          flow.alpha = 0.25,
                          flow.curve = 0.1,
                          fill.col = NULL){
  # get metainfo
  meta <- object@meta.data

  # calculate percent ratio
  ratio.info <- meta %>%
    dplyr::group_by(.data[[sample.name]],.data[[celltype.name]]) %>%
    dplyr::summarise(num = n()) %>%
    dplyr::mutate(rel_num = num/sum(num))

  # color
  if(is.null(fill.col)){
    fill.col <- mycolor[1:length(unique(meta[,celltype.name]))]
  }else{
    fill.col <- fill.col
  }

  # plot
  p <-
    ggplot2::ggplot(ratio.info,
                    ggplot2::aes_string(x = sample.name,y = "rel_num")) +
    ggplot2::geom_col(ggplot2::aes_string(fill = celltype.name),
                      width = col.width) +
    ggalluvial::geom_flow(ggplot2::aes_string(stratum = celltype.name,
                                              alluvium = celltype.name,
                                              fill = celltype.name),
                          width = 0.5,
                          alpha = flow.alpha,
                          knot.pos = flow.curve) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(values = fill.col,
                               name = "Cell Type") +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black'),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black')) +
    ggplot2::xlab('') + ggplot2::ylab('Cell percent ratio')

  return(p)
}
