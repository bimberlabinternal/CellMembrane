#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat


#' @title RunSDA
#'
#' @description This will run SDA on the target assay
#' @param seuratObj A Seurat object.
#' @param outputFolder The path to save results. There will be subfolders for ./rawData and ./results
#' @param numComps Passed to SDAtools::run_SDA(). 30 is a good minimum but depends on input data complexity.
#' @param assayName The name of the assay
#' @param randomSeed Passed to SDAtools::run_SDA() set_seed
#' @param minAsinhThreshold Only features with asinh() on total counts above this value will be included.
#' @param featureInclusionList An optional vector of genes that will be included in SDA
#' @param featureExclusionList An optional vector of genes that will be excluded from SDA
#' @param minLibrarySize Passed to dropsim::normaliseDGE() min_library_size
#' @param path.sda The full path to the SDA binary. By default it assumes sda_static_linux in in your $PATH
#' @param max_iter Passed directly to SDAtools::run_SDA()
#' @param nThreads Passed to SDAtools::run_SDA() num_openmp_threads
#' @export
RunSDA <- function(seuratObj, outputFolder, numComps = 50, assayName = 'RNA', randomSeed = GetSeed(), minAsinhThreshold = 0.5, featureInclusionList = NULL, featureExclusionList = NULL, minLibrarySize = 50, path.sda = "sda_static_linux", max_iter = 10000, nThreads = 1) {
  SerObj.DGE <- seuratObj@assays[[assayName]]@counts
  
  ## default gene inclusion (0.5 is basically including all detectable genes)
  print(paste0('Initial features: ', nrow(SerObj.DGE)))
  inclusionFeats <- rownames(SerObj.DGE)[asinh(Matrix::rowSums(SerObj.DGE)) > minAsinhThreshold]
  print(paste0('After gene count filter: ', length(inclusionFeats)))

  if (!all(is.null(featureInclusionList))) {
    featureInclusionList <- RIRA::ExpandGeneList(featureInclusionList)
    preExisting <- intersect(features, featureInclusionList)
    print(paste0('Adding ', length(featureInclusionList), ' features, of which ', length(preExisting), ' are already present'))
    features <- unique(c(features, featureInclusionList))
    print(paste0('Total after: ', length(features)))
  }

  if (!all(is.null(featureExclusionList))){
    featureExclusionList <- RIRA::ExpandGeneList(featureExclusionList)
    preExisting <- intersect(features, featureExclusionList)
    print(paste0('Excluding ', length(featureExclusionList), ' features(s), of which ', length(preExisting), ' are present'))
    features <- unique(features[!(features %in% featureExclusionList)])
    print(paste0('Total after: ', length(features)))
  }

  P1 <- ggplot(data.frame(x = sqrt(Matrix::colSums(SerObj.DGE[inclusionFeats, ]))), aes(x = x)) +
    geom_density() +
    ggtitle("SQRT(Total transcript per cell)") +
    geom_vline(xintercept = 50, color = 'dodgerblue') +
    geom_vline(xintercept = 100, color = 'orange') +
    geom_vline(xintercept = 200, color = 'gold') +
    geom_vline(xintercept = 800, color = 'red') +
    geom_vline(xintercept = 1600, color = 'green')

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
  normedDGE <- dropsim::normaliseDGE(Matrix::as.matrix(SerObj.DGE[inclusionFeats, ]),
                                     center = FALSE, #dont change
                                     scale = TRUE,#dont change
                                     threshold = 10, #dont change
                                     min_library_size = minLibrarySize, #see above density plot
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

  rawDataDir <- paste0(outputFolder, 'rawData')
  print('Files after save:')
  print(list.files(outputFolder))
  print(list.files(rawDataDir))

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
          data = rawDataDir,
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

  results <- SDAtools::load_results(results_folder = resultsDir, data_path = rawDataDir)

  SDAtools::check_simulation_scores(data = normedDGE, results = results)
  SDAtools::check_convergence(results)
  SDAtools::loading_distribution(results)
  SDAtools::scores_distribution(results)
  SDAtools::plot_maximums(results)
  SDAtools::plot_scree(results)
  SDAtools::PIP_distribution(results)
  SDAtools::PIP_component_distribution(results, 2)
  SDAtools::PIP_threshold_distribution(results)

  return(results)
}
