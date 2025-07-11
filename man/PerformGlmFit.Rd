% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseudoBulk.R
\name{PerformGlmFit}
\alias{PerformGlmFit}
\title{PerformGlmFit}
\usage{
PerformGlmFit(
  seuratObj,
  design,
  test.use = "QLF",
  assayName = "RNA",
  filterGenes = TRUE,
  legacy = FALSE,
  plotBCV = TRUE
)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{design}{The model.matrix object}

\item{test.use}{Can be either QLF or LRT. QLF runs edgeR::glmQLFTest, while LRT runs edgeR::glmLRT}

\item{assayName}{The name of the assay to use}

\item{filterGenes}{A boolean controlling whether or not to filter genes using edgeR::filterByExpr. If TRUE, genes with low counts will be filtered out.}

\item{legacy}{A passthrough variable for edgeR's glmQLF function. They recently (R 4.0) changed the default behavior, so this will break on earlier versions of R.}

\item{plotBCV}{A boolean determining if the BCV plot should be shown.}
}
\value{
An edgeR glm object
}
\description{
Creates an edgeR glm object
}
