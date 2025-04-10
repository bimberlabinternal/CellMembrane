% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeuratUtils.R
\name{SubsetSeurat}
\alias{SubsetSeurat}
\title{Subset Seurat Using String Expressions}
\usage{
SubsetSeurat(seuratObj, expressionStrings = NULL)
}
\arguments{
\item{seuratObj}{The seurat object}

\item{expressionStrings}{A vector strings, each of which will be parsed into an expression, to use in selection}
}
\description{
Subset a seurat object, selecting cells that match each expression
}
