utils::globalVariables(
  names = c('Var1', 'value', 'Fraction'),
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
#' @param nThreads If provided, this integer value is passed to SingleR's BPPARAM argument. On windows ths is passed to BiocParallel::SnowParam(). On other OS it is passed to BiocParallel::MulticoreParam()
#' @param createConsensus If true, a pseudo-consensus field will be created from the course labels from all datasets. Labels will be simplified in an attempt to normalize into the categories of Bcells, NK/T_cells and Myeloid.
#' @return The modified seurat object
#' @importFrom pheatmap pheatmap
#' @import Seurat
#' @import SingleR
#' @export
#' @importFrom scuttle logNormCounts
RunSingleR <- function(seuratObj = NULL, datasets = c('hpca', 'blueprint', 'dice', 'monaco', 'immgen'), assay = NULL, resultTableFile = NULL, rawDataFile = NULL, minFraction = 0.01, showHeatmap = TRUE, maxCellsForHeatmap = 20000, nThreads = NULL, createConsensus = TRUE){
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
      rownames(ref) <- toupper(rownames(ref))
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
    print(paste0('Total genes shared with reference data: ', length(genesPresent), ' of ', nrow(seuratObj)))

    if (length(genesPresent) < 100) {
      print(paste0('Too few shared genes, skipping: ', length(genesPresent)))
      next
    }

    ref <- ref[genesPresent,]

    seuratObjSubset <- Seurat::DietSeurat(seuratObj, assays = c(assay), counts = TRUE)
    seuratObjSubset <- subset(seuratObjSubset, features = genesPresent)
    Seurat::DefaultAssay(seuratObjSubset) <- assay

    #Convert to SingleCellExperiment
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = GetAssayData(object = seuratObjSubset, assay = assay, slot = "counts")))
    sce <- scuttle::logNormCounts(sce)
    rm(seuratObjSubset)

    refAssay <- 'logcounts'
    if (!('logcounts' %in% names(SummarizedExperiment::assays(ref)))) {
      print('logcount not present, using normcounts as assay type')
      refAssay <- 'normcounts'
    }

    if (!is.null(nThreads)) {
      if (.Platform$OS.type == 'windows') {
        BPPARAM <- BiocParallel::SnowParam(nThreads)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(nThreads)
      }
    } else {
      BPPARAM <- BiocParallel::SerialParam()
    }

    tryCatch({
      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main, assay.type.test = 'logcounts', assay.type.ref = refAssay, fine.tune = TRUE, prune = TRUE, BPPARAM = BPPARAM))
      if (length(colnames(seuratObj)) != nrow(pred.results)) {
        stop('Length of SingleR results did not match seurat object')
      }

      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = colnames(seuratObj), classification_type = 'Main', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
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

      suppressWarnings(print(SingleR::plotDeltaDistribution(pred.results)))

      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- colnames(seuratObj)
      fn <- paste0(dataset, '.label')
      allFields <- c(allFields, fn)
      seuratObj[[fn]] <- toAdd

      tab <- table(cluster=as.character(Seurat::Idents(seuratObj)), label=unname(seuratObj[[fn]][[fn]]))
      pheatmap::pheatmap(log10(tab+10), main = dataset) # using a larger pseudo-count for smoothing.

      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = 'logcounts', assay.type.ref = refAssay, fine.tune = TRUE, prune = TRUE))
      if (length(colnames(seuratObj)) != nrow(pred.results)) {
        stop('Length of SingleR results did not match seurat object')
      }

      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = colnames(seuratObj), classification_type = 'Fine', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
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

      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- colnames(seuratObj)
      fn2 <- paste0(dataset, '.label.fine')
      allFields <- c(allFields, fn2)
      seuratObj[[fn2]] <- toAdd

      tab <- table(cluster=as.character(Seurat::Idents(seuratObj)), label=unname(seuratObj[[fn2]][[fn2]]))
      pheatmap::pheatmap(log10(tab+10), main = paste0(dataset, ': Fine Labels')) # using a larger pseudo-count for smoothing.

      if (!is.null(minFraction)){
        for (label in c(fn, fn2)) {
          print(paste0('Filtering ', label, ' below: ', minFraction))
          d <- data.frame(table(Label = unlist(seuratObj[[label]])))
          names(d) <- c('Label', 'Count')
          d$Fraction <- d$Count / sum(d$Count)

          d <- d %>% dplyr::arrange(dplyr::desc(Fraction))
          print(d)
          toRemove <- d$Label[d$Fraction < minFraction]
          if (length(toRemove) > 0) {
            print(paste0('Will remove: ', paste0(toRemove, collapse = ', ')))
          }

          if (length(toRemove) > 0) {
            l <- unlist(seuratObj[[label]])
            names(l) <- colnames(seuratObj)
            l[l %in% toRemove] <- 'Unknown'
            seuratObj[[label]] <- l
          }

          print('After filter:')
          d <- data.frame(table(Label = unlist(seuratObj[[label]])))
          names(d) <- c('Label', 'Count')
          print(d)
        }
      }
    }, error = function(e){
      print(paste0('Error running singleR for dataset: ', dataset))
      print(conditionMessage(e))
      traceback()
    })
  }

  print(paste0('Adding fields: ', paste0(allFields, collapse = ',')))
  df <- data.frame(CellBarcodes = colnames(seuratObj))
  for (fn in allFields){
    df[fn] <- seuratObj@meta.data[[fn]]
  }

  if (!is.null(rawDataFile)) {
    write.table(completeRawData, file = rawDataFile, sep = '\t', row.names = FALSE)
  }

  if (!is.null(resultTableFile)){
    write.table(file = resultTableFile, df, sep = '\t', row.names = F, quote = F)
  }

  if (length(allFields) == 0) {
    print('No singleR calls were added, this probably indicates there were errors with singleR')
  } else {
    if ('tsne' %in% names(seuratObj@reductions) || 'umap' %in% names(seuratObj@reductions)) {
      DimPlot_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)
    }

    Tabulate_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)

    # Create a pseudo consensus:
    if (createConsensus) {
      allFields <- names(seuratObj@meta.data)
      allFields <- allFields[grepl(allFields, pattern = 'label')]

      fieldsToUse <- allFields[!grepl(allFields, pattern = 'fine')]
      dat <- seuratObj@meta.data[fieldsToUse]

      for (colName in names(dat)) {
        dat[[colName]][grepl(dat[[colName]], pattern = "Myeloid|Monocyte|Macrophage", ignore.case = T)] <- "Myeloid"
        dat[[colName]][grepl(dat[[colName]], pattern = "Tcell|T_cell|T cell|T-cell|TCell|T_Cell|T Cell|T-Cell|CD8|CD4|Treg", ignore.case = T)] <- "NK/T_cell"
        dat[[colName]][grepl(dat[[colName]], pattern = "NKcell|NK_cell|NK cell|NK-cell|NKCell|NK_Cell|NK Cell|NK-Cell", ignore.case = T)] <- "NK/T_cell"
        dat[[colName]][grepl(dat[[colName]], pattern = "Bcell|B_cell|B cell|B-cell|BCell|B_Cell|B Cell|B-Cell", ignore.case = T)] <- "B_cell"
        dat[[colName]][grepl(dat[[colName]], pattern = "DC|Dendritic", ignore.case = T)] <- "DCs"
      }

      dat$SingleRConsensus <- sapply(1:nrow(dat), function(idx) {
        vals <- unlist(dat[idx, fieldsToUse, drop = T])
        vals <- unique(vals[!is.na(vals) & vals != 'Unknown'])
        if (length(vals) == 0) {
          return(NA)
        }

        return(paste0(sort(unique(vals)), collapse = ','))
      })

      seuratObj$SingleRConsensus <- dat$SingleRConsensus

      print(DimPlot(seuratObj, group.by = 'SingleRConsensus') + theme_bw() + ggtitle('SingleR Consensus') + theme(legend.position="bottom"))

      print(ggplot(reshape2::melt(table(seuratObj@meta.data$SingleRConsensus)), aes(x=Var1, y = value, fill=Var1))  +
        geom_bar(stat="identity", position="dodge", width = 0.7) +
        theme_bw() +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("SingleR Consensus") +
        ylab("Number of cells") +
        xlab("")
      )
    }
  }

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
        ggtitle(paste0("SingleR Classification: ", dataset)) +
        ylab("Number of cells") + xlab(""),

      ggplot(reshape2::melt(table(seuratObject@meta.data[[paste0(dataset, '.label.fine')]])), aes(x=Var1, y = value, fill=Var1)) +
        geom_bar(stat="identity", position="dodge", width = 0.7) +
        theme_bw() +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle(paste0("SingleR Classification (Fine): ", dataset)) +
        ylab("Number of cells") + xlab("")
    )

    if (plotIndividually) {
      plot(plots[[1]])
      plot(plots[[2]])
    } else {
      print(plots[[1]] + plots[[2]] + patchwork::plot_layout(ncol = 1))
    }
  }
}
