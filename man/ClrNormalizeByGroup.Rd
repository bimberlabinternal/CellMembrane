% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{ClrNormalizeByGroup}
\alias{ClrNormalizeByGroup}
\title{ClrNormalizeByGroup}
\usage{
ClrNormalizeByGroup(
  seuratObj,
  groupingVar,
  assayName = "ADT",
  targetAssayName = NA,
  margin = 1,
  minCellsPerGroup = 20,
  calculatePerFeatureUCell = FALSE,
  featureInclusionList = NULL,
  featureExclusionList = NULL,
  doAsinhTransform = FALSE
)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{groupingVar}{The variable to use to partition the data}

\item{assayName}{The name of the assay}

\item{targetAssayName}{If provided, data will be saved to this assay, rather than modifying the source assay}

\item{margin}{Passed directly to NormalizeData()}

\item{minCellsPerGroup}{If provided, any group with newer than this many cells will be dropped}

\item{calculatePerFeatureUCell}{If TRUE, UCell will be run once per feature in the assay}

\item{featureInclusionList}{If provided, the input assay will be subset to just these features.}

\item{featureExclusionList}{If provided, the input assay will be subset to exclude these features.}

\item{doAsinhTransform}{If true, asinh transform will be performed on the raw counts prior to CLR}
}
\description{
This subsets the input object based on a variable (like dataset), and performs CLR normalization per-group, then combines them
}
