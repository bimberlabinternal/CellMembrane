#' @title generate SingleR object
#'
#' @description Compute SingleR classification on a Seurat object
#' @param seuratObj A Seurat object
#' @param datasets One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, and immgen. See cellDex package for available datasets.
#' @param assay The assay in the seurat object to use
#' @param resultTableFile If provided, a table of results will be saved here
#' @param singlerSavePrefix If provided, the SingleR results will be saved to RDS here
#' @param minFraction If provided, any labels present with fraction of this or fewer across cells will be converted to Unknown
#' @param showHeatmap If true, heatmaps will be generated showing the SingleR calls
#' @param maxCellsForHeatmap The heatmap will only be plotted if the total cells is below this number
#' @return The modified seurat object
#' @import Seurat
#' @import SingleR
#' @export
#' @importFrom scater logNormCounts
RunSingleR <- function(seuratObj = NULL, datasets = c('hpca', 'blueprint', 'dice', 'monaco'), assay = NULL, resultTableFile = NULL, singlerSavePrefix = NULL, minFraction = 0.01, showHeatmap = TRUE, maxCellsForHeatmap = 20000){
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
    print(paste0('Total genes shared with reference data: ', length(genesPresent)))

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
      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main, method = 'single', assay.type.ref = refAssay))
      pred.results$labels[is.na(pred.results$labels)] <- 'Unknown'
      if (!is.null(singlerSavePrefix)){
        saveRDS(pred.results, file = paste0(singlerSavePrefix, '.', dataset, '.singleR.rds'))
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

      toAdd <- pred.results$labels
      names(toAdd) <- rownames(pred.results)
      fn <- paste0(dataset, '.label')
      allFields <- c(allFields, fn)
      seuratObj[[fn]] <- toAdd

      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, method = 'single', assay.type.ref = refAssay))
      pred.results$labels[is.na(pred.results$labels)] <- 'Unknown'
      if (!is.null(singlerSavePrefix)){
        saveRDS(pred.results, file = paste0(singlerSavePrefix, '.', dataset, '.singleR.fine.rds'))
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

      toAdd <- pred.results$labels
      names(toAdd) <- rownames(pred.results)

      fn2 <- paste0(dataset, '.label.fine')
      allFields <- c(allFields, fn2)
      seuratObj[[fn2]] <- toAdd
    }, error = function(e){
      print(paste0('Error running singleR for dataset: ', dataset))
    })

    #sanity check:
    if (length(colnames(seuratObj)) != length(rownames(pred.results))) {
      stop('SingleR did not produce results for all cells')
    }

    if (!is.null(minFraction)){
      for (label in c(fn, fn2)) {
        l <- unlist(seuratObj[[label]])
        names(l) <- colnames(seuratObj)

        print(paste0('Filtering ', label, ' below: ', minFraction))
        d <- table(Label = l)
        print(kableExtra::kbl(d))

        d <- d / sum(d)
        toRemove <- names(d)[d < minFraction]
        if (length(toRemove) > 0) {
            print(paste0('Will remove: ', paste0(toRemove, collapse = ', ')))
        }

        l[l %in% toRemove] <- 'Unknown'
        seuratObj[[label]] <- l

        print('After filter:')
        l <- unlist(seuratObj[[label]])
        d <- table(Label = l)
        print(kableExtra::kbl(d))
      }
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

  return(seuratObj)
}


#' @title DimPlot SingleR Class Labels
#' @description Create a Dimplot from a Seurat object with SingleR class labels
#' @param seuratObject a Seurat object, but if path given, path is prioritized.
#' @param plotIndividually If true, two separate plots will be printed.  Otherwise a single plot wil be printed with one above the other
#' @param datasets One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, and immgen. See cellDex package for available datasets.
#' @keywords Dimplot SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
DimPlot_SingleRClassLabs <- function(seuratObject, plotIndividually = F, datasets = c('hpca')){
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
      print(cowplot::plot_grid(plots[[1]], plots[[2]], ncol = 1))
    }
  }
}


#' @title Tabulate SingleR Class Labels
#' @description Tabulate SingleR class labels from a Seurat object
#' @param seuratObject a Seurat object, but if path given, path is prioritized.
#' @param plotIndividually If true, two separate plots will be printed.  Otherwise a single plot wil be printed with one above the other
#' @param datasets One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, and immgen. See cellDex package for available datasets.
#' @keywords Tabulate SingleR Classification 
#' @export
#' @import Seurat
#' @importFrom cowplot plot_grid
Tabulate_SingleRClassLabs <- function(seuratObject, plotIndividually = F, datasets = c('hpca')) {
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
        cowplot::plot_grid(plots[[1]], plots[[2]], ncol = 1)
    }
  }
}
