#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

utils::globalVariables(
  names = c('Component', 'Score', 'Comp', 'GO_data', 'GeneOdds', 'BgOdds', 'pvalue', 'Count', 'p.adjust', 'Description'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title RunSDA
#'
#' @description This will run SDA on the target assay
#' @param seuratObj A Seurat object.
#' @param outputFolder The path to save results. There will be subfolders for ./rawData and ./results
#' @param numComps Passed to SDAtools::run_SDA(). 30 is a good minimum but depends on input data complexity.
#' @param minCellsExpressingFeature Can be used with perCellExpressionThreshold to drop features present in limited cells. Only features expressed >= perCellExpressionThreshold in at least minCellsExpressingFeature cells will be retained. If this value is less than one it is interpreted as a percentage of total cells. If above one it is interpreted as the min number of cells.
#' @param perCellExpressionThreshold Can be used with perCellExpressionThreshold to drop features present in limited cells. Only features expressed >= perCellExpressionThreshold in at least minCellsExpressingFeature cells will be retained.
#' @param minFeatureCount Only features where the total counts across all cells are greater than or equal to this value will be included. Setting this value to one will include all expressed genes.
#' @param featureInclusionList An optional vector of genes that will be included in SDA
#' @param featureExclusionList An optional vector of genes that will be excluded from SDA
#' @param maxFeaturesDiscarded After filtering, if fewer than this number of features remain, an error will be thrown. This is a guard against overly aggressive filtering. This can either be an integer (meaning number of features), or 0-1 (which is interpreted as a percent of the input features)
#' @param assayName The name of the assay
#' @param randomSeed Passed to SDAtools::run_SDA() set_seed
#' @param minLibrarySize Passed to dropsim::normaliseDGE() min_library_size. Only cells with library size equal or greater to this will be kept. IMPORTANT: this is applied after feature selection.
#' @param path.sda The full path to the SDA binary. By default it assumes sda_static_linux in in your $PATH
#' @param max_iter Passed directly to SDAtools::run_SDA()
#' @param nThreads Passed to SDAtools::run_SDA() num_openmp_threads
#' @param storeGoEnrichment If true, SDA_GO_Enrichment will be performed and stored in the result list
#' @export
RunSDA <- function(seuratObj, outputFolder, numComps = 50, minCellsExpressingFeature = 0.01, perCellExpressionThreshold = 2, minFeatureCount = 1, featureInclusionList = NULL, featureExclusionList = NULL, maxFeaturesDiscarded = NULL, assayName = 'RNA', randomSeed = GetSeed(), minLibrarySize = 50, path.sda = 'sda_static_linux', max_iter = 10000, nThreads = 1, storeGoEnrichment = FALSE) {
  SerObj.DGE <- seuratObj@assays[[assayName]]@counts

  n_cells <- ncol(SerObj.DGE)
  if (n_cells > 250000) {
    stop('SDA has shown to max handle ~200K cells ')
  } else if (n_cells > 150000) {
    warning('SDA has shown to max handle ~200K cells ')
  }

  print(paste0('Initial features: ', nrow(SerObj.DGE)))
  print(paste0('Initial cells: ', ncol(SerObj.DGE)))
  featuresToUse <- rownames(SerObj.DGE)

  if (!is.na(minFeatureCount) && minFeatureCount > 0) {
    numFeatures <- length(featuresToUse)
    featuresToUse <- featuresToUse[Matrix::rowSums(SerObj.DGE[featuresToUse, ]) >= minFeatureCount]
    print(paste0('After filtering to features with total counts of at least ', minFeatureCount, ': ', length(featuresToUse), ' features remain (', scales::percent(length(featuresToUse) / numFeatures),' of input)'))
    rm(numFeatures)
  }

  if (!is.na(minCellsExpressingFeature) && minCellsExpressingFeature > 0) {
    if (is.na(perCellExpressionThreshold)) {
      stop('Must provide perCellExpressionThreshold when minCellsExpressingFeature is above zero')
    }

    print('Filtering on minCellsExpressingFeature')
    if (minCellsExpressingFeature < 1) {
      minCellsExpressingFeatureRaw <- minCellsExpressingFeature
      minCellsExpressingFeature <- floor(minCellsExpressingFeatureRaw * ncol(seuratObj))
      print(paste0('Interpreting minCellsExpressingFeature as a percentage of total cells (', ncol(seuratObj), '), converted from ', minCellsExpressingFeatureRaw, ' to ', minCellsExpressingFeature))
    }

    numFeatures <- length(featuresToUse)
    numNonZeroCells <- Matrix::rowSums(SerObj.DGE >= perCellExpressionThreshold)
    featuresToUse <- names(numNonZeroCells[which(numNonZeroCells >= minCellsExpressingFeature)])
    print(paste0('After limiting to features with expression GTE ', perCellExpressionThreshold, ' in at least ', minCellsExpressingFeature, ' cells: ', length(featuresToUse), ' features remain (', scales::percent(length(featuresToUse) / numFeatures), ' of input)'))
    rm(numFeatures)
  }

  if (!all(is.null(featureInclusionList))) {
    featureInclusionList <- RIRA::ExpandGeneList(featureInclusionList)
    preExisting <- intersect(featuresToUse, featureInclusionList)
    print(paste0('Adding ', length(featureInclusionList), ' features, of which ', length(preExisting), ' are already present'))
    featuresToUse <- unique(c(featuresToUse, featureInclusionList))
    print(paste0('Total after: ', length(featuresToUse)))
  }

  if (!all(is.null(featureExclusionList))){
    featureExclusionList <- RIRA::ExpandGeneList(featureExclusionList)
    preExisting <- intersect(featuresToUse, featureExclusionList)
    print(paste0('Excluding ', length(featureExclusionList), ' features(s), of which ', length(preExisting), ' are present'))
    featuresToUse <- unique(featuresToUse[!(featuresToUse %in% featureExclusionList)])
    print(paste0('Total after: ', length(featuresToUse)))
  }

  if (length(featuresToUse) == 0) {
    stop('No features remain after filtering')
  }

  if (!is.null(maxFeaturesDiscarded)) {
    if (maxFeaturesDiscarded < 1) {
      maxFeaturesDiscarded <- maxFeaturesDiscarded * nrow(SerObj.DGE)
    }

    featsDiscarded <- nrow(SerObj.DGE) - length(featuresToUse)
    if (featsDiscarded > maxFeaturesDiscarded) {
      stop(paste0('The total number of features discarded, ', featsDiscarded, ' exceeds the threshold of ', maxFeaturesDiscarded))
    }
  }

  df <- data.frame(x = Matrix::colSums(SerObj.DGE[featuresToUse, ]))
  dens <- stats::density(df$x)
  mx <- dens$x[which.max(dens$y)]

  P1 <- ggplot(df, aes(x = x)) +
    scale_x_sqrt() +
    geom_density() +
    ggtitle('Total features per cell') +
    labs(x = 'Features/Cell', y = 'Density') +
    geom_vline(xintercept = mx, color = 'red') +
    ggtitle(paste0('Library Size: Peak = ', mx))


  if (!is.null(minLibrarySize)) {
    P1 <- P1 + geom_vline(xintercept = minLibrarySize, color = 'red')
  }

  print(P1)

  ### other methods work, perhaps we can add other options in the future
   # most cases works but can be taken as input depeding on how the density plot above looks
  print('starting dropsim normaliseDGE')
  normedDGE <- dropsim::normaliseDGE(Matrix::as.matrix(SerObj.DGE[featuresToUse, ]),
                                     center = FALSE,
                                     scale = TRUE,
                                     threshold = 10, #dont change. Any normalised values larger than this will be rounded down
                                     min_library_size = minLibrarySize, # Cells with library size below this value (for the features selected) will be dropped. See above density plot
                                     gene_subset = 1)

  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }

  if (!endsWith(outputFolder, '/')) {
    outputFolder <- paste0(outputFolder, '/')
  }

  print(paste0('Saving raw data to: ', outputFolder))
  normedDGE <- as.matrix(normedDGE)
  SDAtools::export_data(normedDGE, path = outputFolder, name = 'rawData')

  rawDataFile <- paste0(outputFolder, 'rawData')
  resultsDir <- paste0(outputFolder, 'results/')

  print(paste0('Saving results to: ', resultsDir))
  if (dir.exists(resultsDir)) {
    message(paste0('Deleting existing result folder: ', resultsDir))
    unlink(resultsDir, recursive = TRUE)
  }

  print('Running SDA')
  if (!file.exists(path.sda)) {
    x <- unname(Sys.which(path.sda))
    if (x != '') {
      print(paste0('Found SDA under PATH: ', x))
      path.sda <- x
    }
  }

  SDAtools::run_SDA(sda_location = path.sda,
          out = resultsDir,
          data = rawDataFile,
          num_comps = numComps,
          max_iter = max_iter,
          save_freq = 1000,
          set_seed = randomSeed, #TODO: consider allowing a vector of seeds for replicates?
          N = n_cells,
          eigen_parallel = (nThreads > 1),
          ignore_missing = FALSE,
          num_blocks = 8,
          num_openmp_threads = nThreads
  )

  results <- SDAtools::load_results(results_folder = resultsDir, data_path = outputFolder)
  results$CellBarcodes <- colnames(normedDGE)
  results$Features <- rownames(normedDGE)
  results <- .AddCompStats(results)
  if (storeGoEnrichment) {
    results$goEnrichment <- SDA_GO_Enrichment(results, components = 1:nrow(results$loadings[[1]]))
  }

  tryCatch({
    SDAtools::check_convergence(results)
    SDAtools::loading_distribution(results)
    SDAtools::scores_distribution(results)
  }, error = function(e){
    print(paste0('Error generating SDA first plots'))
    print(conditionMessage(e))
  })

  tryCatch({
    SDAtools::plot_scree(results)
  }, error = function(e){
    print(paste0('Error generating SDA plot_scree'))
    print(conditionMessage(e))
  })

  tryCatch({
    SDAtools::PIP_distribution(results)
  }, error = function(e){
    print(paste0('Error generating SDA PIP_distribution'))
    print(conditionMessage(e))
  })

  tryCatch({
    SDAtools::PIP_threshold_distribution(results)
  }, error = function(e){
    print(paste0('Error generating SDA PIP_threshold_distribution'))
    print(conditionMessage(e))
  })

  return(results)
}

.EnsureFeaturesInSdaResults <- function(sdaResults) {
  if (!'Features' %in% names(sdaResults) || !'CellBarcodes' %in% names(sdaResults)) {
    stop('Features and/or CellBarcodes have not been added to the sdaResults object. This should have been performed if the object was created by CellMembrane')
  }
}

Plot_CorSDA_Loadings <- function(results) {
  rownames(results$loadings[[1]]) <- paste0('SDAV', 1:nrow(results$loadings[[1]]))

  pheatmap::pheatmap(cor(t(results$loadings[[1]][,])))
}


#' @title Plot SDA Cell Scores By Feature
#'
#' @param seuratObj A Seurat object
#' @param sdaResults The result list generated by SDATools
#' @param metadataFeature The field to plot (present in the seurat object).
#' @param direction Either Pos, Neg, or Both
#' @export
Plot_SDAScoresPerFeature <- function(seuratObj, sdaResults, metadataFeature, direction = 'Both'){
  .EnsureFeaturesInSdaResults(sdaResults)

  if (!metadataFeature %in% names(seuratObj@meta.data)) {
    print(paste0('Feature not found, skipping: ', metadataFeature))
    return()
  }

  if (direction == 'Both') {
    directions <- c('Pos', 'Neg')
  } else {
    directions <- direction
  }

  for (direction in directions) {
    SDAScores <- sdaResults$scores
    MetaDF <- seuratObj[[metadataFeature, drop = FALSE]]
    MetaDF <- MetaDF[rownames(SDAScores),,drop = FALSE]

    CompsDF <- as.data.frame(lapply(levels(factor(MetaDF[,1,drop = TRUE])), function(CondX){
      apply(SDAScores[rownames(MetaDF)[which(MetaDF[,1,drop = TRUE] == CondX)], ,drop = FALSE], 2,
            function(x){
              if (direction == 'Neg'){
                round(sum(x<0)/nrow(SDAScores)*100, 2)
              } else if (direction == 'Pos'){
                round(sum(x>0)/nrow(SDAScores)*100, 2)
              } else {
                stop(paste0('Direction must be Pos or Neg: ', direction))
              }
            })
    }))

    colnames(CompsDF) <- levels(factor(MetaDF[,1]))
    CompsDF <- as.data.frame(CompsDF[naturalsort::naturalsort(rownames(CompsDF)),])
    CompsDF$SDA <- factor(rownames(CompsDF), levels=rownames(CompsDF))
    ChiT <- stats::chisq.test(CompsDF[,1:(ncol(CompsDF)-1)])
    ChiTres <- ChiT$residuals
    ChiTres[which(is.na(ChiTres))] <- 0
    ChiResSD <- round(apply(ChiTres, 1, sd),2)
    ChiResSD[which(is.na(ChiResSD))] <- 0
    ChiResSD[ChiResSD < 0.2] <- ''

    HM <- ComplexHeatmap::Heatmap((t(ChiTres)),
                                  name = direction,
                                  column_title = paste0('Direction: ', direction),
                                  cluster_columns = TRUE,
                                  cluster_rows = TRUE,
                                  col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ='RdBu')))(10),
                                  column_labels = paste0(rownames(CompsDF), ' sd_', ChiResSD)
    )

    print(HM)
  }
}

#' @title SDAToSeuratMetadata
#'
#' @description This will store the SDA scores in the seurat object metadata table
#' @param seuratObj A Seurat object.
#' @param results The SDA results list
#' @param plotComponents If true, FeaturePlots will be generated for the components
#' @export
SDAToSeuratMetadata <- function(seuratObj, results, plotComponents = TRUE){
  # Note: this adds NAs for missing cells. We could in theory change this to zeros if NAs are a problem.
  SDAScores <- results$scores[rownames(seuratObj@meta.data), ]
  colnames(SDAScores) <- paste0('SDA_', 1:ncol(SDAScores))

  for (cmp in  colnames(SDAScores)){
    seuratObj <- Seurat::AddMetaData(seuratObj, asinh(SDAScores[,cmp]), col.name = cmp)
    if (plotComponents){
      print(Seurat::FeaturePlot(seuratObj, features = cmp, order = T) & ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
    }
  }

  return(seuratObj)
}

#' @title SDAToSeuratReduction
#'
#' @description This will store the SDA results as a reduction in the seurat object
#' @param seuratObj A Seurat object.
#' @param sdaResults The SDA results list
#' @param assayName The source assay
#' @param reduction.name The name used for this reduction
#' @param reduction.key The key used for this reduction
#' @export
SDAToSeuratReduction <- function(seuratObj, sdaResults, assayName = 'RNA', reduction.name = 'sda', reduction.key = 'SDA_') {
  # rows = cells
  embeddings <- sdaResults$scores
  colnames(embeddings) <- paste0(reduction.key, 1:ncol(embeddings))

  # Note: since SDA could drop cells, add back in the missing cells with zeros
  extraCells <- setdiff(rownames(embeddings), colnames(seuratObj))
  if (length(extraCells) > 0) {
    stop(paste0('There were ', length(extraCells), ' with data in the SDA results but not present in the seurat object.  Top barcodes: ', paste0(head(extraCells), collapse = ',')))
  }

  missingCells <- setdiff(colnames(seuratObj), rownames(embeddings))
  if (length(missingCells) > 0) {
    toAdd <- matrix(rep(0, ncol(embeddings)*length(missingCells)), ncol = ncol(embeddings))
    colnames(toAdd) <- colnames(embeddings)
    rownames(toAdd) <- missingCells
    embeddings <- rbind(embeddings, toAdd)
  }

  embeddings <- embeddings[colnames(seuratObj),]

  # rows = features. This needs to use the projected slot, since it probably has fewer features than the assay
  projected <- t(sdaResults$loadings[[1]])
  colnames(projected) <- paste0(reduction.key, 1:ncol(projected))

  # See: https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette.html
  sda.reduction <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    projected = projected,
    key = reduction.key,
    assay = assayName
  )

  for (x in c('component_statistics', 'n', 'goEnrichment')) {
    if (x %in% names(sdaResults)) {
      sda.reduction@misc[[x]] <- sdaResults[[x]]
    }
  }

  seuratObj[[reduction.name]] <- sda.reduction

  return(seuratObj)
}


.AddCompStats <- function(SDAres, redoCalc = T){
  if (redoCalc) {
    SDAres$component_statistics <- NULL
  }

  if (is.null(SDAres$component_statistics) ) {
    SDAres$component_statistics <- data.frame(
      Component = 1:SDAres$n$components,
      Component_name = dimnames(SDAres$scores)[[2]],
      max_score = apply(abs(SDAres$scores),  2, max),
      max_loading = apply(abs(SDAres$loadings[[1]]), 1, max),
      mean_score = apply(abs(SDAres$scores),  2, mean),
      mean_loading = apply(abs(SDAres$loadings[[1]]), 1, mean),
      sd_score = apply(abs(SDAres$scores),  2, sd),
      sd_loading = apply(abs(SDAres$loadings[[1]]), 1, sd),
      ssqrd_score = apply(SDAres$scores,  2, function(x) sum(x^2)),
      ssqrd_loading = apply(SDAres$loadings[[1]], 1, function(x) sum(x^2))
    ) %>% dplyr::arrange(-Component)
  } else {
    print('component_statistics found, will not redo calculations. See redoCalc to repeat calculations')
  }

  SDAres$component_statistics <- data.frame(SDAres$component_statistics)

  return(SDAres)
}


# Note: this is not well tested, but should be possible. Will need work.
#' @title Run UMAP using SDA results
#'
#' @param seuratObj A Seurat object.
#' @param dimsToUse The number of dims to use.  If null, this will be inferred using FindSeuratElbow()
#' @param minDimsToUse The minimum numbers of dims to use.  If dimsToUse is provided, this will override.
#' @param reduction The name of the source reduction
#' @return A modified Seurat object.
#' @export
RunUmapOnSDA <- function(seuratObj, dimsToUse = NULL, minDimsToUse = NULL, reduction = 'sda') {
  seuratObj <- .RunUMAP(seuratObj, reduction = reduction, dimsToUse = dimsToUse, minDimsToUse = minDimsToUse, reduction.name = 'sda.umap', reduction.key = 'sdaUMAP_')
  print(DimPlot(seuratObj, reduction = 'sda.umap'))

  return(seuratObj)
}


#' @title Plot SDA Cell Scores
#'
#' @param sdaResults The result list generated by SDATools
#' @param seuratObj A Seurat object
#' @param fieldNames A vector of field names (present in the seurat object) to plot. These must be non-numeric fields.
#' @export
PlotSdaCellScores <- function(sdaResults, seuratObj, fieldNames) {
  dat <- as.data.frame(t(sdaResults$scores))
  dat$Comp <- naturalsort::naturalfactor(colnames(sdaResults$scores))
  dat <- dat %>% tidyr::pivot_longer(names_to = 'CellBarcode', cols = rownames(sdaResults$scores), values_to = 'Score')

  meta <- seuratObj@meta.data
  meta$CellBarcode <- rownames(meta)
  dat <- merge(dat, meta, by = 'CellBarcode', all.x = T, all.y = F)


  for (fieldName in fieldNames) {
    if (!fieldName %in% colnames(dat)) {
      print(paste0('Field not found, skipping: ', fieldName))
      next
    }
    else if (is.numeric(dat[[fieldName]])) {
      print(paste0('Plotting not supported for numeric fields, skipping: ', fieldName))
      next
    }

    maxBarsPerRow <- 16
    totalValues <- length(unique(dat[[fieldName]]))
    nComponents <- length(unique(dat[['Comp']]))

    nCol <- ceiling(maxBarsPerRow / totalValues)
    totalRows <- ceiling(nComponents / nCol)

    rowsPerPage <- 3
    perPage <- nCol * rowsPerPage
    totalPages <- ceiling(nComponents / perPage)

    for (i in 1:totalPages) {
      P1 <- ggplot(dat, aes_string(y = 'Score', x = fieldName)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1,alpha = 0.3) +
        egg::theme_presentation(base_size = 14) +
        ggtitle(fieldName) +
        theme(
          legend.position = 'none',
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggforce::facet_wrap_paginate(facets = . ~ Comp, scales = 'free_y', drop = FALSE, ncol = nCol, nrow = min(3, totalRows), page = i)

      print(P1)
    }
  }
}


#' @title Compute and score GO enrichment for SDA components
#'
#' @param sdaResults The result list generated by SDATools
#' @param components A vector of SDA components to inspect
#' @param orgDb The Bioconductor OrgDb (string), which is passed directed to clusterProfiler::enrichGO()
#' @param geneNumber The number of genes to inspect
#' @param direction Whether to include positive genes, negative, or both
#' @return A list mapping component name + direction to the enrichGO results
#' @export
SDA_GO_Enrichment <- function(sdaResults, components, orgDb = 'org.Hs.eg.db', geneNumber = 100, direction = 'Both'){
  ret <- list()

  if (direction == 'Both') {
    directions <- c('Pos', 'Neg')
  } else {
    directions <- direction
  }

  for (comp in components) {
    for (direction in directions) {
      print(paste0('Loading component ', comp, ', direction: ', direction))
      if (direction == 'Neg'){
        top_genes <- names(sort(sdaResults$loadings[[1]][comp,]))[1:geneNumber]
      } else if (direction == 'Pos'){
        top_genes <- names(sort(sdaResults$loadings[[1]][comp,], decreasing = T))[1:geneNumber]
      }

      gene_universe <- names(sdaResults$loadings[[1]][comp,])
      ego <- clusterProfiler::enrichGO(
                      gene = top_genes,
                      universe = gene_universe,
                      OrgDb = orgDb,
                      keyType = 'SYMBOL',
                      ont = 'BP',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 1,
                      qvalueCutoff = 1
       )

      frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))

      ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio) / frac_to_numeric(ego@result$BgRatio)
      ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, '/', fixed = TRUE), function(x){
        x <- as.numeric(x) ;x[1] / (x[2]-x[1])
      }))

      ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, '/', fixed = TRUE), function(x){
        x <- as.numeric(x) ;x[1] / (x[2]-x[1])
      }))

      ret[[paste0(comp, '-', direction)]] <- ego@result
    }
  }

  return(ret)
}


SDA_GO_VolcanoPlot <- function(x = GO_data, component = 'V5N', plotTitle = paste0('Component : ', component)){
  print(ggplot(data.frame(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size = Count)) +
    geom_point(aes(colour=p.adjust<0.05)) +
    scale_size_area() +
    ggrepel::geom_label_repel(data = data.frame(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size = 0.25), size = 3, force=2) +
    ggtitle(plotTitle) +
    xlab('Odds Ratio') +
    scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
  )
}