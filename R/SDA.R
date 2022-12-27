#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat


#' @title RunSDA
#'
#' @description This will run SDA on the target assay
#' @param seuratObj A Seurat object.
#' @param outputFolder The path to save results. There will be subfolders for ./rawData and ./results
#' @param numComps Passed to SDAtools::run_SDA(). 30 is a good minimum but depends on input data complexity.
#' @param minCellsExpressingFeature. Can be used with perCellExpressionThreshold to drop features present in limited cells. Only features detected above perCellExpressionThreshold in at least minCellsExpressingFeature will be retained. If this value is less than zero it is interpreted as a percentage of total cells. If above zero it is interpeted as the min number of cells.
#' @param perCellExpressionThreshold Can be used with perCellExpressionThreshold to drop features present in limited cells. Only features detected above perCellExpressionThreshold in at least minCellsExpressingFeature will be retained.
#' @param minFeatureCount Only features where their total counts across all cells are above this value will be included.
#' @param featureInclusionList An optional vector of genes that will be included in SDA
#' @param featureExclusionList An optional vector of genes that will be excluded from SDA
#' @param assayName The name of the assay
#' @param randomSeed Passed to SDAtools::run_SDA() set_seed
#' @param minLibrarySize Passed to dropsim::normaliseDGE() min_library_size. Only cells with library size equal or greater to this will be kept. IMPORTANT: this is applied after feature selection.
#' @param path.sda The full path to the SDA binary. By default it assumes sda_static_linux in in your $PATH
#' @param max_iter Passed directly to SDAtools::run_SDA()
#' @param nThreads Passed to SDAtools::run_SDA() num_openmp_threads
#' @export
RunSDA <- function(seuratObj, outputFolder, numComps = 50, minCellsExpressingFeature = 0.01, perCellExpressionThreshold = 0, minFeatureCount = 200, featureInclusionList = NULL, featureExclusionList = NULL, assayName = 'RNA', randomSeed = GetSeed(), minLibrarySize = 50, path.sda = "sda_static_linux", max_iter = 10000, nThreads = 1) {
  SerObj.DGE <- seuratObj@assays[[assayName]]@counts
  
  ## default gene inclusion (0.5 is basically including all detectable genes)
  featuresToUse <- rownames(SerObj.DGE)
  print(paste0('Initial features: ', length(featuresToUse)))

  if (!is.na(minCellsExpressingFeature)) {
    print('Filtering on minCellsExpressingFeature')
    if (minCellsExpressingFeature < 1) {
      minCellsExpressingFeatureRaw <- minCellsExpressingFeature
      minCellsExpressingFeature <- floor(minCellsExpressingFeatureRaw * ncol(seuratObj))

      print(paste0('Interpreting minCellsExpressingFeature as a percentage of total cells, converted from ', minCellsExpressingFeature, ' to ', minCellsExpressingFeature))
    }

    numNonZeroCells <- Matrix::rowSums(SerObj.DGE > perCellExpressionThreshold)
    featuresToUse <- names(num.cells[which(numNonZeroCells >= minCellsExpressingFeature)])
    print(paste0('After filtering to features with expression > ', perCellExpressionThreshold, ' in at least ', minCellsExpressingFeature, ' cells: ', length(featuresToUse)))
  }

  if (!is.na(minFeatureCount)) {
    featuresToUse <- rownames(SerObj.DGE)[Matrix::rowSums(SerObj.DGE) > minFeatureCount]
    print(paste0('After gene count filter of ', minFeatureCount, ': ', length(featuresToUse)), ' features remain')
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

  P1 <- ggplot(data.frame(x = sqrt(Matrix::colSums(SerObj.DGE[featuresToUse, ]))), aes(x = x)) +
    geom_density() +
    ggtitle("SQRT(Total transcript per cell)") +
    geom_vline(xintercept = 50, color = 'dodgerblue') +
    geom_vline(xintercept = 100, color = 'orange') +
    geom_vline(xintercept = 200, color = 'gold') +
    geom_vline(xintercept = 800, color = 'red') +
    geom_vline(xintercept = 1600, color = 'green') +
    ggtitle('Library Size')

  print(P1)

  n_cells <- ncol(SerObj.DGE)
  if (n_cells > 250000) {
    stop("SDA has shown to max handle ~200K cells ")
  }
  else if (n_cells > 150000) {
    warning("SDA has shown to max handle ~200K cells ")
  }

  ### other methods work, perhaps we can add other options in the future
   # most cases works but can be taken as input depeding on how the density plot above looks
  print("starting dropsim normaliseDGE")
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
  rownames(results$loadings[[1]]) <- paste0("SDAV", 1:nrow(results$loadings[[1]]))

  pheatmap::pheatmap(cor(t(results$loadings[[1]][,])))
}

Plot_SDAScoresPerFeature <- function(seuratObj, sdaResults, metadataFeature, direction = "Neg"){
  .EnsureFeaturesInSdaResults(sdaResults)

  SDAScores <- sdaResults$scores
  MetaDF <- seuratObj[[c(metadataFeature), drop = FALSE]]
  MetaDF <- MetaDF[rownames(SDAScores),,drop = FALSE]

  # EISA: why do the two round functions use different decimal precision? 1 for neg and 2 for pos??
  if (direction == "Neg"){
    CompsDF <- as.data.frame(lapply(levels(factor(MetaDF[,1])), function(CondX){
      apply(SDAScores[rownames(MetaDF)[which(MetaDF[,1] == CondX)], ], 2,
            function(x){
              round(sum(x<0)/nrow(SDAScores)*100, 1)
            })
    }))
  } else if(direction == "Pos"){
    CompsDF <- as.data.frame(lapply(levels(factor(MetaDF[,1])), function(CondX){
      apply(SDAScores[rownames(MetaDF)[which(MetaDF[,1] == CondX)], ], 2,
            function(x){
              round(sum(x>0)/nrow(SDAScores)*100, 2)
            })
    }))
  } else {
    stop("Direction needs to be Neg or Pos")
  }

  colnames(CompsDF) <- levels(factor(MetaDF[,1]))
  CompsDF <- as.data.frame(CompsDF[gtools::mixedsort(rownames(CompsDF)),])
  CompsDF$SDA <- factor(rownames(CompsDF), levels=rownames(CompsDF))
  ChiT <- chisq.test(CompsDF[,1:(ncol(CompsDF)-1)])
  ChiTres <- ChiT$residuals
  ChiTres[which(is.na(ChiTres))] <- 0
  ChiResSD <- round(apply(ChiTres, 1, sd),2)
  ChiResSD[which(is.na(ChiResSD))] <- 0
  ChiResSD[ChiResSD < 0.2] <- ""
  
  pheatmap::pheatmap((t(ChiTres)),
                     cluster_cols = TRUE, cluster_rows = TRUE,
                     color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(10),
                     labels_col = paste0(rownames(CompsDF), " sd_", ChiResSD)
  )
}

#' @title SDAToSeuratMetadata
#'
#' @description This will store the SDA scores in the seurat object metadata table
#' @param seuratObj A Seurat object.
#' @param results The SDA results list
#' @param plotComponents If true, FeaturePlots will be generated for the components
#' @export
SDAToSeuratMetadata <- function(seuratObj, results, plotComponents = TRUE){
  # EISA: so this adds NAs for missing genes, right? Is that the desired behavior?
  SDAScores <- results$scores[rownames(seuratObj@meta.data), ]
  colnames(SDAScores) <- paste0("SDA_", 1:ncol(SDAScores))

  for (cmp in  colnames(SDAScores)){
    seuratObj <- Seurat::AddMetaData(seuratObj, asinh(SDAScores[,cmp]), col.name = cmp)
    if (plotComponents){
      print(Seurat::FeaturePlot(seuratObj, features = cmp, order = T) & ggplot2::scale_colour_gradientn(colours = c("navy", "dodgerblue", "white", "gold", "red")))
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
  # EISA: does adding cells with zeros make sense? Seurat requires the matrix to have values for all cells
  # rows = cells. Note: since SDA could drop cells, add back in the missing cells with zeros
  embeddings <- sdaResults$scores
  colnames(embeddings) <- paste0(reduction.key, 1:ncol(embeddings))

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

  seuratObj[[reduction.name]] <- sda.reduction

  return(seuratObj)
}


# EISA: note, I removed maxscoreThrsh = 20, maxloadThrsh = 1, since they were not used.
.AddCompStats <- function(SDAres, sdThrsh = 0.04, redoCalc = T){
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

  # EISA: what is this designed to do? Why do we need these other variables?
  SDAres$component_statistics$Component_namev2 <- SDAres$component_statistics$Component_name
  SDAres$component_statistics$Component_name_plot <- SDAres$component_statistics$Component_name

  # EISA: what is this designed to do? This seems to be setting names to empty string?
  #SDAres$component_statistics[which(SDAres$component_statistics$sd_loading<sdThrsh),]$Component_namev2 <- rep("", length(which(SDAres$component_statistics$sd_loading<sdThrsh)))

  SDAres$component_statistics <- data.frame(SDAres$component_statistics)

  return(SDAres)
}


# EISA: this is not well tested, but should be possible. Will need work.
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

# EISA: not adding these yet due to the dependencies and they seem like something easy enough to run locally for the time being
# GO_enrichment <- function(results = NULL, component, geneNumber = 100, threshold=0.01, side="N", OrgDb = NULL){
#   # results = SDAres; component = 1; geneNumber = 100; threshold=0.01; side="N"; OrgDb = RefGenome
#   require(data.table)
#   if(side=="N"){
#     top_genes <- data.frame(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(V1)][1:geneNumber]$rn
#   }else{
#     top_genes <- data.frame(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(-V1)][1:geneNumber]$rn
#   }
#   gene_universe <- data.frame(as.matrix(results$loadings[[1]][component,]), keep.rownames = TRUE)$rn
#   print(head(top_genes))
#   ego <- enrichGO(gene = top_genes,
#                   universe = gene_universe,
#                   OrgDb = OrgDb,
#                   keyType = 'SYMBOL',
#                   ont = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff = 1,
#                   qvalueCutoff = 1)
#   ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio)/frac_to_numeric(ego@result$BgRatio)
#   ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
#   ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
#   return(ego@result)
# }
#
# frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))
#
#
# go_volcano_plot <- function(x=GO_data, component="V5N", extraTitle=""){
#   if(extraTitle=="") extraTitle = paste("Component : ", component, sep="")
#   #print(
#   ggplot(data.frame(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size = Count)) +
#     geom_point(aes(colour=p.adjust<0.05)) +
#     scale_size_area() +
#     geom_label_repel(data = data.frame(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size = 0.25), size = 3, force=2) +
#     ggtitle(paste("",extraTitle, sep="\n") ) +
#     xlab("Odds Ratio") +
#     scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
#   #)
# }
#
#
# plotGeneLocations <- function() {
#   library(biomaRt)
#   if(file.exists("./data/sda/primeseq/RIRA4_TNK_SDA1_biomaRt_gene_loc_human.rds")){
#     ens_ds="hsapiens_gene_ensembl"
#     gene_locations <- get.location(gene.symbols=colnames(results$loadings[[1]]),
#                                    data_set = ens_ds,
#                                    gene_name = "external_gene_name")
#     saveRDS(gene_locations, "./data/sda/primeseq/RIRA4_TNK_SDA1_biomaRt_gene_loc_human.rds")
#   } else {
#     gene_locations = readRDS("./data/sda/primeseq/RIRA4_TNK_SDA1_biomaRt_gene_loc_human.rds")
#   }
# }
#
# get.location <- function(gene.symbols, data_set, gene_name){
#   ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = data_set)
#   mapTab <- getBM(attributes = c(gene_name,'chromosome_name','start_position'),
#                   filters = gene_name, values = gene.symbols, mart = ensembl, uniqueRows=TRUE)
#   mapTab <- as.data.frame(mapTab)
#   setnames(mapTab, gene_name,"gene_symbol")
#   # Remove duplicate genes!!
#   # first which genes are duplicated
#   duplicates <- mapTab[duplicated(mapTab, by="gene_symbol")]$gene_symbol
#   # order by chr so patch versions to go bottom, then choose first unique by name
#   unduplicated <- unique(mapTab[gene_symbol %in% duplicates][order(chromosome_name)], by="gene_symbol")
#   # remove duplicates and replace with unique version
#   mapTab <- mapTab[!gene_symbol %in% duplicates]
#   mapTab <- rbind(mapTab, unduplicated)
#   # change all patch chr names to NA
#   mapTab[!chromosome_name %in% c(1:22,"X","Y","MT")]$chromosome_name <- NA # mice actually 19
#   mapTab[is.na(chromosome_name)]$start_position <- NA
#   return(mapTab)
# }