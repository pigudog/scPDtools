% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCP-workflow.R
\name{check_srtMerge}
\alias{check_srtMerge}
\title{Check and preprocess a merged seurat object}
\usage{
check_srtMerge(
  srtMerge,
  batch = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  vars_to_regress = NULL,
  seed = 11,
  ...
)
}
\arguments{
\item{srtMerge}{A merged \code{Seurat} object with batch information.}

\item{batch}{Metadata column name containing the batch information.}

\item{do_normalization}{Whether to normalize the data. If NULL, will automatically determine.}

\item{normalization_method}{Normalization method.Can be one of "LogNormalize", "SCT".}

\item{do_HVF_finding}{Whether to find the high variable features(HVF). If NULL, will automatically determine.}

\item{HVF_source}{Source of the HVF. Can be one of "separate" and "global".}

\item{nHVF}{HVF number to use.}

\item{HVF}{Custom high variable features.}

\item{vars_to_regress}{Variables to regress out.}

\item{seed}{Set a random seed.}

\item{...}{}
}
\description{
Check and preprocess a merged seurat object
}
