% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_III.R
\name{ScoreCellCycle}
\alias{ScoreCellCycle}
\title{Score Cell Cycle}
\usage{
ScoreCellCycle(seuratObj, min.genes = 10, useAlternateG2M = FALSE)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{min.genes}{If less than min.genes are shared between the seurat object and the reference cell cycle genes, this method will abort.}

\item{useAlternateG2M}{If true, this will use a smaller set of G2M genes, defined from: https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv}
}
\value{
A modified Seurat object.
}
\description{
Score Cell Cycle
}
