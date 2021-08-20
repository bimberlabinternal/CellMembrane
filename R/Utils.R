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

#' @title Rename a vector of genes using CD nomenclature
#'
#' @description This compares a vector of genes to CD nomenclature, based on https://www.genenames.org/data/genegroup/#!/group/471. Any gene symbol that matches will have the CD name appended to the end (i.e. KLRB1 becomes "KLRB1 (CD161)")
#' @param inputGenes A vector of genes to be aliased.
#' @export
RenameUsingCD <- function(inputGenes) {
  vals <- CellMembrane::cdGenes[c('GeneSymbol', 'PreviousSymbols', 'Synonyms')]
  vals$Aliases <- sapply(paste0(vals$PreviousSymbols, ',', vals$Synonyms), function(x){
    x <- unlist(strsplit(x, split = ','))
    x <- x[x != '']
    x <- x[grep(x, pattern = '^CD')]
    x <- sort(x)
    x <- unique(toupper(x))

    return(paste0(x, collapse = ','))
  })
  vals <- vals[!is.na(vals$Aliases) & vals$Aliases != '',]
  vals <- vals[c('GeneSymbol','Aliases')]
  vals

  toMerge <- data.frame(Name = inputGenes, SortOrder = 1:length(inputGenes))
  toMerge <- merge(toMerge, vals, by.x = 'Name', by.y = 'GeneSymbol', all.x = T)
  sel <- !is.na(toMerge$Aliases)
  toMerge$Name[sel] <- paste0(toMerge$Name[sel], ' (', toMerge$Aliases[sel], ')')
  toMerge <- arrange(toMerge, SortOrder)

  return(toMerge$Name)
}