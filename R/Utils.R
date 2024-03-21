#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table

utils::globalVariables(
  names = c('X', 'Y', 'SortOrder'),
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

#' @title Get S Phase Genes
#' @export
#' @return The default list of S phase genes
GetSPhaseGenes <- function(){
  return (.GetCCGenes()[1:43])
}

#' @title Get G2M Phase Genes
#' @param alternateGenes If true, the genes will be based on: https://hbctraining.github.io/scRNA-seq_online/lessons/cell_cycle_scoring.html
#' @export
#' @return The default list of G2M phase genes
GetG2MGenes <- function(alternateGenes = FALSE) {
  if (alternateGenes) {
    # based on: https://hbctraining.github.io/scRNA-seq_online/lessons/cell_cycle_scoring.html
    # https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv
    return(c('ANLN','ANP32E','AURKA','AURKB','BIRC5','BUB1','CBX5','CCNB2','CDC20','CDC25C','CDCA2','CDCA3','CDCA8','CDK1','CENPA','CENPE','CENPF','CKAP2','CKAP2L','CKAP5','CKS1B','CKS2','CTCF','DLGAP5','ECT2','G2E3','GAS2L3','GTSE1','HJURP','HMGB2','HMMR','JPT1','KIF11','KIF20B','KIF23','KIF2C','LBR','MKI67','NCAPD2','NDC80','NEK2','NUF2','NUSAP1','PIMREG','PSRC1','RANGAP1','SMC4','TACC3','TMPO','TOP2A','TPX2','TTK','TUBB4B','UBE2C'))
  }

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
  }

  return(seuratObj)
}


#' @title Update Gene Model
#'
#' @description Substitutes LOC genes to more common gene IDs
#' @param features a vector of features to be updated to more common gene IDS
#' @param predictions Whether or not to include lower quality/speculative updates to the gene model

.UpdateGeneModel <- function(features, predictions = T){
# TRAC = LOC710951 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJ95TNA1016&from=84561541&to=84565299)
features <- gsub("LOC710951", "TRAC", features)

# TRBC1 = LOC114677140 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJ9S5K0Z01R&from=169372796&to=169374268)
# TRBC2 maps to LOC114677140 as well
features <- gsub("LOC114677140", "TRBC1", features)
#features <- gsub("LOC114677140", ,"TRBC2", features)

#TRDC = LOC711031 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJA2XTBW01R&from=84480217&to=84482757)
features <- gsub("LOC711031", "TRDC", features)

# TRGC1 = LOC720538 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJAHHBHJ016&from=76878214&to=76881635)
features <- gsub("LOC720538", "TRGC1", features)

#TRGC2 = LOC705095 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJASX36T016&from=76909393&to=76914320)
features <- gsub("LOC705095", "TRGC2", features)

#TRDV1 = LOC720456 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=2&RID=UJB64E6Z013&from=84055100&to=84055663)
#TRDV1 also maps to TRDC (LOC710951) so use that if the variable region seems odd to include
features <- gsub("LOC720456", "TRDV1", features)

#IGKC = LOC701504 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041766.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBE65YT013&from=18130540&to=18130862)
features <- gsub("LOC701504", "IGKC", features)

#IGHG3 = LOC114679691 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBJKUU9016&from=168001652&to=168005818) 
#These are gamma-2 like instead of gamma-3 like, but they map well. 
features <- gsub("LOC114679691", "IGHG3", features)

#IGHG1 = LOC708891 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBVCMEA013&from=168070075&to=168071681)
features <- gsub("LOC708891", "IGHG1", features)

#IGHA1 = LOC720839 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJC2E4S8013&from=167910270&to=167911722)
features <- gsub("LOC720839", "IGHA1", features)

#HLA.DP1 = LOC114669810 (source: https://www.ncbi.nlm.nih.gov/nuccore/XM_028842593.1)
features <- gsub("LOC114669810", "HLA.DPA1", features)

#FCN1 = LOC712405 (source: https://www.ncbi.nlm.nih.gov/nuccore/XM_015116399.2)
features <- gsub("LOC712405", "FCN1", features)

#XCL1 = LOC100423131 (source: https://www.ncbi.nlm.nih.gov/gene/100423131)
features <- gsub("LOC100423131", "XCL1.2", features)

#CCL4 = LOC100430627 (source: https://www.ncbi.nlm.nih.gov/gene/100430627)
#LOC100426632 also maps to CCL4
features <- gsub("LOC100430627", "CCL4", features)

#CCL3L1 = LOC100426537 (source: https://www.ncbi.nlm.nih.gov/gene/100426537)
#alternative name: CCL3L3
features <- gsub("LOC100426537", "CCL3L1", features)

  #these are predicted proteins/isoforms that seem less annotated
  if(predictions){
    #TYMP = LOC716161 (source: https://www.ncbi.nlm.nih.gov/protein/1622837354)
    features <- gsub("LOC716161", "TYMP", features)
  }

return(features)
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

#' @title Update Seurat Object Gene Model From Monkey LOCs to Human Gene IDs
#'
#' @description Uses UpdateGeneModel to change LOC genes to more common human gene IDs in the RNA assay
#' @param seuratObj The Seurat Object to be updated
#' @param verbose This prints the slot and assay number, primarily for debugging once seurat updates.
#' @export

UpdateMacaqueMmul10NcbiGeneSymbols <- function(seuratObj, verbose = T){
  if (class(seuratObj)[[1]] != "Seurat"){
    stop(paste0('Please provide a Seurat Object'))
  }
  for (assay in names(seuratObj@assays)){
    for (slot in names(attributes(seuratObj@assays[[assay]]))){
      if (verbose){
        print(paste("Updating Gene Names in Assay:", assay, "Slot:", slot))
      }

      #updating common gene expression slots
      if (any(grepl(slot, c("counts", "data", "scale.data")))){
        ad <- Seurat::GetAssayData(seuratObj, assay = assay, slot = slot)
        rownames(ad) <- .UpdateGeneModel(rownames(ad))
      }

      #updating variable features (expected to be a vector)
      if (slot == "var.features"){
        attr(x = seuratObj@assays[[assay]], which = slot) <- .UpdateGeneModel(attr(x = seuratObj@assays[[assay]], which = slot))
      }

      #updating meta.features (expected to be a matrix)
      if (slot == "meta.features"){
        rownames(attr(x = seuratObj@assays[[assay]], which = slot)) <- .UpdateGeneModel(rownames(attr(x = seuratObj@assays[[assay]], which = slot)))
        #feature metadata is tracked by assay, so this matrix could be empty.
        #This is an update that seems to be CellMembrane specific
        if(dim(attr(x = seuratObj@assays[[assay]], which = slot))[[2]] != 0 & any(grepl("GeneId", names(attr(x = seuratObj@assays[[assay]], which = slot)))))
          attr(x = seuratObj@assays[[assay]], which = slot)$GeneId <- .UpdateGeneModel(attr(x = seuratObj@assays[[assay]], which = slot)$GeneId)
      }
    }
  }
  #this block attempts to update the features in PCA reductions. (expected to be two matrices)
  if (any(grepl("pca", names(seuratObj@reductions)))){
    rownames(seuratObj@reductions$pca@feature.loadings) <- .UpdateGeneModel(rownames(seuratObj@reductions$pca@feature.loadings))
    rownames(seuratObj@reductions$pca@feature.loadings.projected) <- .UpdateGeneModel(rownames(seuratObj@reductions$pca@feature.loadings.projected))
  }
  return(seuratObj)
}

.SplitIntoBatches <- function(vect, batchSize) {
  numBatches <- ceiling(length(vect) / batchSize)

  ret <- list()

  for (batchNum in 1:numBatches) {
    start0 <- (batchNum-1)*batchSize
    end <- min(length(vect), start0 + batchSize)

    ret[[batchNum]] <- vect[(start0+1):end]
  }

  return(ret)
}

#' @title Resolve NCBI Loc Genes
#'
#' @description The accepts a vector of NCBI LOC gene IDs and returns a data frame with their description and aliases, queries from rentrez
#' @param geneIds A vector of gene IDs, all of which must start with LOC
#' @param maxBatchSize The max number of IDs for each query against NCBI
#' @export
ResolveLocGenes <- function(geneIds, maxBatchSize = 100) {
  invalidIds <- geneIds[!grepl(geneIds, pattern = '^LOC')]
  if (length(invalidIds) > 0) {
    stop('All genes must start with LOC: ', paste0(invalidIds, collapse = ';'))
  }

  genesBatched <- .SplitIntoBatches(geneIds, batchSize = maxBatchSize)
  print(paste0('Total batches: ', length(genesBatched)))

  ret <- NULL
  batchNum <- 0
  for (geneBatch in genesBatched) {
    batchNum <- batchNum + 1
    print(paste0('Querying batch ', batchNum, ' of ', length(genesBatched)))

    results <- rentrez::entrez_summary(db="gene", id = gsub(geneBatch, pattern = '^LOC', replacement = ''))
    df <- do.call(rbind, lapply(results, FUN = function(x){
      return(data.frame(Name = x$name, Description = x$description, Aliases = x$otheraliases))
    }))

    rownames(df) <- paste0('LOC', rownames(df))
    df$GeneId <- rownames(df)

    if (all(is.null(ret))) {
      ret <- df
    } else {
      ret <- rbind(ret, df)
    }
  }

  # Ensure return order matches input:
  ret <- ret[geneIds,]

  return(ret)
}

#' @title ClrNormalizeByGroup
#'
#' @description This subsets the input object based on a variable (like dataset), and performs CLR normalization per-group, then combines them
#' @param seuratObj The seurat object
#' @param groupingVar The variable to use to partition the data
#' @param assayName The name of the assay
#' @param targetAssayName If provided, data will be saved to this assay, rather than modifying the source assay
#' @param margin Passed directly to NormalizeData()
#' @param minCellsPerGroup If provided, any group with newer than this many cells will be dropped
#' @param calculatePerFeatureUCell If TRUE, UCell will be run once per feature in the assay
#' @param featureInclusionList If provided, the input assay will be subset to just these features.
#' @param featureExclusionList If provided, the input assay will be subset to exclude these features.
#' @param doAsinhTransform If true, asinh transform will be performed on the raw counts prior to CLR
#' @export
ClrNormalizeByGroup <- function(seuratObj, groupingVar, assayName = 'ADT', targetAssayName = NA, margin = 1, minCellsPerGroup = 20, calculatePerFeatureUCell = FALSE, featureInclusionList = NULL, featureExclusionList = NULL, doAsinhTransform = FALSE) {
  if (!groupingVar %in% names(seuratObj@meta.data)) {
    stop(paste0('Field not found: ', groupingVar))
  }

  if (!assayName %in% names(seuratObj@assays)) {
    stop(paste0('Source assay not found: ', assayName))
  }

  if (!is.null(minCellsPerGroup) && !is.na(minCellsPerGroup)) {
    groupValues <- as.character(seuratObj[[groupingVar, drop = TRUE]])
    if (any(is.na(groupValues))) {
      stop(paste0('There were NAs for the column: ', groupingVar))
    }

    totals <- table(groupValues)
    totals <- totals[totals < minCellsPerGroup]
    toDrop <- character()
    if (length(totals) > 0) {
      for (groupName in names(totals)) {
        d <- colnames(seuratObj)[groupValues == groupName]
        if (any(is.na(d))) {
          stop(paste0('There were NAs for the column: ', groupingVar, ' with value: ', groupName))
        }
        print(paste0('Dropping group: ', groupName, ' with total cells: ', length(d)))
        toDrop <- c(toDrop, d)
      }

      if (length(toDrop) > 0) {
        cellsBefore <- ncol(seuratObj)
        print(paste0('Dropping total of ', length(toDrop), ' cells, out of ', cellsBefore))
        seuratObj <- subset(seuratObj, cells = toDrop, invert = TRUE)
        print(paste0('Remaining: ', ncol(seuratObj)))
        if (ncol(seuratObj) != (cellsBefore - length(toDrop))) {
          stop(paste0('subset did not work as expected. initial cells: ', cellsBefore, ', to drop: ', length(toDrop), ', cells remaining: ', ncol(seuratObj)))
        }
      }
    }
  }

  sourceAssay <- assayName
  if (!is.na(targetAssayName) && !is.null(targetAssayName)) {
    seuratObj@assays[[targetAssayName]] <- Seurat::CreateAssayObject(Seurat::GetAssayData(seuratObj, assay = sourceAssay, slot = 'counts'))
    seuratObj@assays[[targetAssayName]]@key <- paste0(tolower(targetAssayName), '_')
    sourceAssay <- targetAssayName
  }

  groups <- unique(seuratObj[[groupingVar, drop = TRUE]])
  normalizedMat <- NULL

  for (groupName in groups) {
    cells <- colnames(seuratObj)[!is.na(seuratObj@meta.data[[groupingVar]]) & seuratObj@meta.data[[groupingVar]] == groupName]
    print(paste0('Processing group: ', groupName, ' with ', length(cells), ' cells'))
    ad <- subset(seuratObj@assays[[sourceAssay]], cells = cells)
    if (ncol(ad) != length(cells)) {
      stop(paste0('Incorrect assay subset for group: ', groupName, '. Expected: ', length(cells), ', actual: ', ncol(ad)))
    }

    if (!all(is.null(featureInclusionList))) {
      featureInclusionList <- RIRA::ExpandGeneList(featureInclusionList)
      toKeep <- intersect(rownames(ad), featureInclusionList)
      print(paste0('Limiting to ', length(featureInclusionList), ' features, of which ', length(toKeep), ' exist in this assay'))
      if (length(toKeep) == 0) {
        stop(paste0('None of the featureInclusionList features were found in this object: ', paste0(featureInclusionList, collapse = ',')))
      }
      ad <- subset(ad, features = toKeep)
      if (nrow(ad) != length(toKeep)) {
        stop(paste0('Incorrect assay subset. Expected: ', length(toKeep), ', actual: ', nrow(ad)))
      }
      print(paste0('Total features after: ', nrow(ad)))
    }

    if (!all(is.null(featureExclusionList))){
      featureExclusionList <- RIRA::ExpandGeneList(featureExclusionList)
      toDrop <- intersect(rownames(ad), featureExclusionList)
      print(paste0('Excluding ', length(featureExclusionList), ' features(s) from the input assay, of which ', length(toDrop), ' exist in this assay'))
      if (length(toDrop) == 0) {
        stop(paste0('None of the featureExclusionList features were found in this object: ', paste0(featureExclusionList, collapse = ',')))
      }

      featuresToKeep <- rownames(ad)[!rownames(ad) %in% toDrop]
      ad <- subset(ad, features = featuresToKeep)
      if (nrow(ad) != length(featuresToKeep)) {
        stop(paste0('Incorrect assay subset. Expected: ', length(featuresToKeep), ', actual: ', nrow(ad)))
      }
      print(paste0('Total features after: ', nrow(ad)))
    }

    if (doAsinhTransform) {
      dat <- Seurat::GetAssayData(ad, slot = 'counts')
      # based on flowCore module: https://github.com/RGLab/flowCore/blob/1dee3931c7ac922052b74fcdf6ba037fe1313892/R/AllClasses.R#L5030
      # and ADTnorm: https://github.com/yezhengSTAT/ADTnorm/blob/d5d08d9cc075d1300d0ff1038ff2a3efae780b15/R/arcsinh_transform.R
      a <- 1
      b <- 1/5
      c <- 0
      dat <- asinh(a+b*dat) + c

      dat <- Seurat::NormalizeData(Seurat::as.sparse(dat), normalization.method = 'CLR', margin = margin, verbose = FALSE)
      ad <- Seurat::SetAssayData(ad, slot = 'data', new.data = dat)
    } else {
      ad <- Seurat::NormalizeData(ad, normalization.method = 'CLR', margin = margin, verbose = FALSE)
    }

    if (all(is.null(normalizedMat))) {
      normalizedMat <- Seurat::GetAssayData(ad, slot = 'data')
    } else {
      normalizedMat <- cbind(normalizedMat, Seurat::GetAssayData(ad, slot = 'data'))
    }
  }

  normalizedMat <- normalizedMat[,colnames(seuratObj)]

  assayObj <- Seurat::GetAssay(seuratObj, assay = sourceAssay)
  if (nrow(assayObj) != nrow(normalizedMat) || any(rownames(assayObj)  != rownames(normalizedMat))) {
    rawCounts <- Seurat::GetAssayData(assayObj, slot = 'counts')
    rawCounts <- rawCounts[rownames(normalizedMat),]
    assayObj <- Seurat::CreateAssayObject(counts = rawCounts)
  }

  assayObj <- Seurat::SetAssayData(assayObj, slot = 'data', new.data = as.sparse(normalizedMat))
  seuratObj[[sourceAssay]] <- NULL  # Reset this first to avoid warning
  seuratObj[[sourceAssay]] <- assayObj
  seuratObj <- ScaleData(seuratObj, verbose = FALSE, assay = sourceAssay)

  if (calculatePerFeatureUCell) {
    seuratObj <- CalculateUcellPerFeature(seuratObj, assayName = sourceAssay, columnPrefix = paste0(sourceAssay, '.'))
  }

  return(seuratObj)
}

.InferBpParam <- function(nThreads, defaultValue = NULL) {
  if (!is.null(nThreads) && nThreads > 1) {
    if (.Platform$OS.type == 'windows') {
      return(BiocParallel::SnowParam(nThreads))
    } else {
      return(BiocParallel::MulticoreParam(nThreads))
    }
  } else {
    return(defaultValue)
  }
}

# Copied from pheatmap: https://rdrr.io/cran/pheatmap/src/R/pheatmap.r
scale_mat <- function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat <- switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

scale_rows <- function(x){
  m <- apply(x, 1, mean, na.rm = T)
  s <- apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

#' @title GetAssayMetadataSlotName
#'
#' @description Returns the slotname holding the assay metadata, compatible with seurat 4 and 5
#' @param assayObj The assay object
#' @export
GetAssayMetadataSlotName <- function(assayObj) {
  if (class(assayObj)[1] == 'Assay') {
    return('meta.features')
  } else if (class(assayObj)[1] == 'Assay5') {
    return('meta.data')
  } else {
    stop(paste0('Unknown class: ', class(assayObj)[1]))
  }
}

# This is patterned after Seurat: https://github.com/satijalab/seurat/blob/41d19a8a55350bff444340d6ae7d7e03417d4173/R/utilities.R#L1455C1-L1470C4
# Seurat's Pseudobulking converts numeric columns, like cluster names, with prefixes like this
.CheckColnamesAreNumeric <- function(mat, prefix = "g") {
  col.names <- colnames(mat)
  if (any(!(grepl("^[a-zA-Z]|^\\.[^0-9]", col.names)))) {
    col.names <- ifelse(
      !(grepl("^[a-zA-Z]|^\\.[^0-9]", col.names)),
      paste0(prefix, col.names),
      col.names
    )
    colnames(mat) <- col.names
  }

  return(mat)
}

#' @title GetMsigdbGeneSet
#'
#' @description a slightly extended escape::getGeneSet wrapper to deal with one-step gene set parsing, which can't parse hierarchical gene sets (C5 + GO:BP) and non-hierarchical gene sets (hallmark, C2 itself, etc) at the same time. 
#' @param msigdbGeneSets a character vector of gene sets that is either a top level msigdb gene set name (e.g. "H" for hallmark or "C2" for curated gene sets), or a common hierarchical gene set in a large top level (e.g. C5;BP for GO:BP annotations.)
#' @return a named list of gene sets fetched by msigdbr. 
GetMsigdbGeneSet <- function(msigdbGeneSets = "H") {
  #initialize gene set
  GS <- c()
  #ensure msigdbGeneSets is formatted properly to parse. 
  if (all(!is.null(msigdbGeneSets), !is.na(msigdbGeneSets), !(length(msigdbGeneSets) == 0))) {
    #require vector. customGeneSets needs to be a list, so this could be a point of confusion
    if (!is.vector(msigdbGeneSets)) {
      stop("msigdbGeneSets is not a vector. Please coerce it to a vector of supported characters.")
    }
    ##GO:BP and other hierarchical gene sets require subcategories to be passed with them, so we should parse those individually and concatenate afterwards. I think each of these need to be individually supported.
    #check if all gene sets are non-hierarchical
    if (all(msigdbGeneSets %in% c("H", paste0("C", 1:8)))) {
      #fetch non-hierarchical gene sets
      GS <- c(GS, escape::getGeneSets(library = msigdbGeneSets))
    } else if ("GO:BP" %in% msigdbGeneSets & all(msigdbGeneSets[msigdbGeneSets != "GO:BP"] %in% c("H", paste0("C", 1:8)))) {
      #remove GO:BP from the list and fetch hierarchical gene set.
      msigdbGeneSets <- msigdbGeneSets[msigdbGeneSets != "GO:BP"]
      GS_GO_BP <- escape::getGeneSets(library = "C5", subcategory = "BP")
      #fetch non-hierarchical gene sets and concatenate gene sets
      GS <- c(GS, c(escape::getGeneSets(library = msigdbGeneSets), GS_GO_BP))
    } else {
      #if the msigdb gene set is hierarchical but unsupported, throw error: 
      unsupportedGeneSets <- msigdbGeneSets[!(msigdbGeneSets %in% c("GO:BP", "H", paste0("C", 1:8)))]
      stop(paste0(unsupportedGeneSets, " in msigdbGeneSets are unsupported. Please ensure msigdbGeneSets is any of: H, GO:BP, ", paste0(", C", 1:8), ". Please ensure your requested geneSet is listed, add the gene set to the customGeneSets argument, or contact members of the Bimber Lab to add a gene set."))
    }
  }
}
