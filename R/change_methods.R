#' Convert between data objects
#'
#' This function converts between data format frequently used to hold single
#' cell data
#'
#' @param obj Input a scRNAobject (Maybe a anndata, Seurat object)
#' @param from Format of input object, e.g. "anndata", "seurat", "sce", "loom",
#'   etc (str)
#' @param to Format of output object (str)
#' @param main_layer Required by some formats, may be "counts", "data",
#'   "scale.data", etc (str)
#'
#' @return Output object
#' @export
convertFormat <- function(obj,
                          from = c("anndata", "seurat", "sce", "loom"),
                          to = c("anndata", "loom", "sce", "seurat", "cds"),
                          outFile = NULL,
                          main_layer = NULL, ...) {
  # match.arg is a base function
  from <- match.arg(from)
  to <- match.arg(to)
  # function parse:text
  # character vector. The text to parse. Elements are treated as if they were lines of a file. Other R objects will be coerced to character if possible.
  tryCatch(
    {
      func <- eval(parse(text = paste(from, to, sep = "2")))
    },
    error = function(e) {
      stop(paste0('Unsupported conversion from "', from, '" to "', to, '"'), call. = FALSE)
    },
    finally = {}
  )

  return(func(obj, outFile = outFile, main_layer = main_layer, ...))
}
