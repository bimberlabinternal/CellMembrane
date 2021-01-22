utils::globalVariables(
  names = c('Var1', 'value'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title Run SingleR For A Seurat Object
#'
#' @description Compute SingleR classification on a Seurat object
#' @param seuratObj A Seurat object
#' @param datasets One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, and immgen. See cellDex package for available datasets.
#' @param assay The assay in the seurat object to use
#' @param resultTableFile If provided, a table of results will be saved here
#' @param rawDataFile If provided, the complete SingleR results will be saved to this file
#' @param minFraction If provided, any labels present with fraction of this or fewer across cells will be converted to Unknown
#' @param showHeatmap If true, heatmaps will be generated showing the SingleR calls
#' @param maxCellsForHeatmap The heatmap will only be plotted if the total cells is below this number
#' @return The modified seurat object
#' @import Seurat
#' @import SingleR
#' @export
#' @importFrom scater logNormCounts
RunSingleR <- function(seuratObj = NULL, datasets = c('hpca', 'blueprint', 'dice', 'monaco'), assay = NULL, resultTableFile = NULL, rawDataFile = NULL, minFraction = 0.01, showHeatmap = TRUE, maxCellsForHeatmap = 20000){
  if (is.null(seuratObj)){
      stop("Seurat object is required")
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seuratObj)
  }

  if (length(seuratObj@assays[[assay]]@counts) == 0) {
    print('Selected assay has no count data, trying RNA')
    assay <- 'RNA'
    if (length(seuratObj@assays[[assay]]@counts) == 0) {
      warning('Unable to find counts for the seurat object, aborting SingleR')
      return(seuratObj)
    }
  }

  allFields <- c()
  completeRawData <- NULL

  for (dataset in datasets) {
    print(paste0('Adding dataset: ', dataset))
    if (dataset == 'hpca'){
        ref <- celldex::HumanPrimaryCellAtlasData()
		} else if (dataset == 'immgen') {
      ref <- celldex::ImmGenData()
    } else if (dataset == 'blueprint') {
      ref <- celldex::BlueprintEncodeData()
		} else if (dataset == 'dice') {
      ref <- celldex::DatabaseImmuneCellExpressionData()
		} else if (dataset == 'monaco') {
      ref <- celldex::MonacoImmuneData()
    } else {
      stop(paste0('unknown reference dataset: ', dataset))
    }

		#Subset genes:
    genesPresent <- intersect(rownames(seuratObj@assays[[assay]]), rownames(ref))
    ref <- ref[genesPresent,]

    seuratObjSubset <- Seurat::DietSeurat(seuratObj, assays = c(assay), counts = T)
    seuratObjSubset <- subset(seuratObj, features = genesPresent)

    Seurat::DefaultAssay(seuratObjSubset) <- assay
    print(paste0('Total genes shared with reference data: ', length(genesPresent), ' of ', nrow(seuratObj)))

    if (length(genesPresent) < 100) {
      print(paste0('Too few shared genes, skipping: ', length(genesPresent)))
      next
    }

    #Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(seuratObjSubset, assay = assay)
    sce <- scater::logNormCounts(sce)
    rm(seuratObjSubset)

    refAssay <- 'logcounts'
    if (!('logcounts' %in% names(SummarizedExperiment::assays(ref)))) {
      refAssay <- 'normcounts'
    }

    tryCatch({
      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main, assay.type.ref = refAssay, fine.tune = TRUE, prune = TRUE))
      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = rownames(pred.results), classification_type = 'Main', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
        if (is.null(completeRawData)) {
          completeRawData <- toBind
        } else {
          completeRawData <- rbind(completeRawData, toBind)
        }
      }

      if (showHeatmap) {
        cells.use <- NULL
        if (ncol(seuratObj) > maxCellsForHeatmap) {
          cells.use <- sample(1:ncol(seuratObj), size = maxCellsForHeatmap)
        }

        print(SingleR::plotScoreHeatmap(pred.results, cells.use = cells.use))
      }

      if (sum(colnames(seuratObj) != rownames(pred.results)) > 0) {
        stop('Cell barcodes did not match for all results')
      }

      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- rownames(pred.results)
      fn <- paste0(dataset, '.label')
      allFields <- c(allFields, fn)
      seuratObj[[fn]] <- toAdd

      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.ref = refAssay))
      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = rownames(pred.results), classification_type = 'Fine', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
        if (is.null(completeRawData)) {
          completeRawData <- toBind
        } else {
          completeRawData <- rbind(completeRawData, toBind)
        }
      }

      if (sum(colnames(seuratObj) != rownames(pred.results)) > 0) {
        stop('Cell barcodes did not match for all results')
      }

      if (!is.null(minFraction)){
        for (label in c(fn, fn2)) {
          l <- unlist(seuratObj[[label]])
          names(l) <- colnames(seuratObj)

          print(paste0('Filtering ', label, ' below: ', minFraction))
          d <- data.frame(table(Label = l))
          names(d) <- c('Label', 'Count')
          print(d)

          d <- d / sum(d)
          toRemove <- names(d)[d < minFraction]
          if (length(toRemove) > 0) {
            print(paste0('Will remove: ', paste0(toRemove, collapse = ', ')))
          }

          l[l %in% toRemove] <- 'Unknown'
          seuratObj[[label]] <- l

          print('After filter:')
          l <- unlist(seuratObj[[label]])
          d <- data.frame(table(Label = l))
          names(d) <- c('Label', 'Count')
          print(d)
        }
      }

      if (showHeatmap) {
        cells.use <- NULL
        if (ncol(seuratObj) > maxCellsForHeatmap) {
          cells.use <- sample(1:ncol(seuratObj), size = maxCellsForHeatmap)
        }

        print(SingleR::plotScoreHeatmap(pred.results, cells.use = cells.use))
      }

      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- rownames(pred.results)

      fn2 <- paste0(dataset, '.label.fine')
      allFields <- c(allFields, fn2)
      seuratObj[[fn2]] <- toAdd
    }, error = function(e){
      print(paste0('Error running singleR for dataset: ', dataset))
      print(conditionMessage(e))
    })

    #sanity check:
    if (length(colnames(seuratObj)) != length(rownames(pred.results))) {
      stop('SingleR did not produce results for all cells')
    }

    if (!is.null(rawDataFile)) {
      write.table(completeRawData, file = rawDataFile, sep = '\t', row.names = FALSE)
    }
  }

  print(paste0('Adding fields: ', paste0(allFields, collapse = ',')))
  df <- data.frame(CellBarcodes = rownames(pred.results))
  for (fn in allFields){
    df[fn] <- seuratObj@meta.data[[fn]]
  }
  if (!is.null(resultTableFile)){
    write.table(file = resultTableFile, df, sep = '\t', row.names = F, quote = F)
  }

  DimPlot_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)
  Tabulate_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)

  return(seuratObj)
}


DimPlot_SingleR <- function(seuratObject, plotIndividually = F, datasets = c('hpca')){
  for (dataset in datasets) {
    fn <- paste0(dataset, '.label')
    if (!(fn %in% colnames(seuratObject@meta.data))) {
      print(paste0('dataset not found: ', dataset))
      next
    }
    plots <- list(
      DimPlot(seuratObject, group.by = fn) + theme_bw() + ggtitle(paste0('SingleR Classification: ', dataset)) + theme(legend.position="bottom"),
      DimPlot(seuratObject, group.by = paste0(dataset, '.label.fine')) + theme_bw() + ggtitle(paste0('SingleR Classification, ', dataset, ' (Fine)')) + theme(legend.position="bottom")
    )

    if (plotIndividually){
      print(plots[[1]])
      print(plots[[2]])
    } else {
      print(plots[[1]] + plots[[2]] + patchwork::plot_layout(ncol = 1))
    }
  }
}


#' @import Seurat
Tabulate_SingleR <- function(seuratObject, plotIndividually = F, datasets = c('hpca', 'blueprint', 'dice', 'monaco')) {
  for (dataset in datasets) {
    fn <- paste0(dataset, '.label')
    if (!(fn %in% colnames(seuratObject@meta.data))) {
      print(paste0('dataset not found: ', dataset))
      next
    }

    plots <- list(
      ggplot(reshape2::melt(table(seuratObject@meta.data[[fn]])), aes(x=Var1, y = value, fill=Var1))  +
        geom_bar(stat="identity", position="dodge", width = 0.7) +
        # scale_fill_manual(values=col_vector) +
        theme_bw() +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle(paste0("SingleR Classification:", dataset)) +
        ylab("Number of cells"),

      ggplot(reshape2::melt(table(seuratObject@meta.data[[paste0(dataset, '.label.fine')]])), aes(x=Var1, y = value, fill=Var1)) +
        geom_bar(stat="identity", position="dodge", width = 0.7) +
        theme_bw() +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle(paste0("SingleR Classification (Fine):", dataset)) +
        ylab("Number of cells")
    )

    if (plotIndividually) {
      plot(plots[[1]])
      plot(plots[[2]])
    } else {
      print(plots[[1]] + plots[[2]] + patchwork::plot_layout(ncol = 1))
    }
  }
}
