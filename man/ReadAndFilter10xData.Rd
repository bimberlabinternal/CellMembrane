% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_III.R
\name{ReadAndFilter10xData}
\alias{ReadAndFilter10xData}
\title{Read and Filter 10X files.}
\usage{
ReadAndFilter10xData(
  dataDir,
  datasetId,
  datasetName = NULL,
  emptyDropNIters = 10000,
  emptyDropsFdrThreshold = 0.001,
  storeGeneIds = TRUE,
  emptyDropsLower = 100,
  useEmptyDropsCellRanger = FALSE,
  nExpectedCells = 8000,
  useSoupX = FALSE,
  previouslyFilteredMatrix = NULL
)
}
\arguments{
\item{dataDir}{Either the directory holding raw count data (generally the raw_feature_bc_matrix), or the parent 'outs' dir from cellranger}

\item{datasetId}{This will be used as a prefix for barcodes, and stored in metadata. Also used as the project name for the Seurat object.}

\item{datasetName}{An optional print-friendly name that will be stored in metadata}

\item{emptyDropNIters}{The number of iterations to use with PerformEmptyDrops()}

\item{emptyDropsFdrThreshold}{The FDR threshold to call cells in emptyDrops()}

\item{storeGeneIds}{If true, a map to translate geneId and name (by default rownames will use gene name)}

\item{emptyDropsLower}{Passed directly to emptyDrops(). The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.}

\item{useEmptyDropsCellRanger}{If TRUE, will use DropletUtils emptyDropsCellRanger instead of emptyDrops}

\item{nExpectedCells}{Only applied if emptyDropsCellRanger is selected. Passed to n.expected.cells argument}

\item{useSoupX}{If true, SoupX will be run against the run input data, instead of emptyDrops}

\item{previouslyFilteredMatrix}{An optional filepath to a pre-filtered count matrix in h5 format. If non-null, this file will be read instead of dataDir. Empty drops and/or soupX will be skipped.}
}
\value{
A Seurat object.
}
\description{
Reads in 10X files using Read10X and filters abberent cells using PerformEmptyDropletFiltering and returns a Seurat object.
}
