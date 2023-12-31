% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{loom2anndata}
\alias{loom2anndata}
\title{Convert Loom object to AnnData}
\usage{
loom2anndata(
  inFile,
  outFile = NULL,
  main_layer = c("spliced", "unspliced"),
  obs_names = "CellID",
  var_names = "Gene"
)
}
\arguments{
\item{inFile}{Path to the input Loom file on disk (str)}

\item{outFile}{Save output AnnData to this file if specified (str or NULL)}

\item{main_layer}{Assay to be converted to AnnData.X, can be "spliced" or
"unspliced", default "spliced" (str)}

\item{obs_names}{Column in cell metadata table that contains cell IDs (str)}

\item{var_names}{Column in gene metadata table that contains gene names}
}
\value{
AnnData object
}
\description{
This function converts a Loom object to an Anndata object
}
