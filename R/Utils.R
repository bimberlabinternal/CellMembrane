#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table

utils::globalVariables(
  names = c('X', 'Y'),
  package = 'CellMembrane',
  add = TRUE
)

pkg.env <- new.env(parent=emptyenv());

pkg.env$RANDOM_SEED <- 1234
set.seed(pkg.env$RANDOM_SEED)

.GetCCGenes <- function(){
  # Cell cycle genes were obtained from the Seurat example (See regev_lab_cell_cycle_genes.txt)
  # and stored using use_data(internal = T) (https://github.com/r-lib/usethis and use_data)
  # cc.genes
  # g2m.genes.orig

  return(cc.genes)
}

.GetSPhaseGenes <- function(){
  return (.GetCCGenes()[1:43])
}

.GetG2MGenes <- function() {
  return(unique(c(g2m.genes.orig, .GetCCGenes()[44:97])))
}

.WriteLogMsg <- function(msg, prefixTime = TRUE, file = 'CellMembrane.log.txt') {
  if (prefixTime) {
    msg <- paste0(Sys.time(), ' ', msg)
  }

  write(msg, file = file, append = T)
}

.InferPerplexityFromSeuratObj <- function(seuratObj, perplexity = 30) {
  return(.InferPerplexity(ncol(seuratObj), perplexity))
}

.InferPerplexity <- function(sampleNumber, perplexity = 30) {
  if (sampleNumber - 1 < 3 * perplexity) {
    print(paste0('Perplexity is too large for the number of samples: ', sampleNumber))
    perplexityNew <- floor((sampleNumber - 1) / 3)
    print(paste0('lowering from ', perplexity, ' to: ', perplexityNew))
    perplexity <- perplexityNew
  }

  return(perplexity)
}

.PossiblyAddBarcodePrefix <- function(seuratObj, datasetId, datasetName = NULL) {
  if (!('BarcodePrefix' %in% names(seuratObj@meta.data))) {
    print(paste0('Adding barcode prefix: ', datasetId))
    seuratObj <- RenameCells(object = seuratObj, add.cell.id = datasetId)
    seuratObj[['BarcodePrefix']] <- c(datasetId)
    seuratObj[['DatasetId']] <- c(datasetId)
    if (!is.null(datasetName)) {
      seuratObj[['DatasetName']] <- datasetName
    }
  } else {
    print('Barcode prefix already added')
  }

  return(seuratObj)
}

#' @title Set random seed
#'
#' @description Sets the seed used for R‘s random number generator, which should be used in all internal functions
#' @param seed The random seed
#' @export
SetSeed <- function(seed) {
  pkg.env$RANDOM_SEED <- seed
  set.seed(pkg.env$RANDOM_SEED)
}

#' @title Get random seed
#'
#' @description Sets a random seed, which should be used in all internal functions
#' @param seed The random seed
#' @export
GetSeed <- function(seed) {
  return(pkg.env$RANDOM_SEED)
}