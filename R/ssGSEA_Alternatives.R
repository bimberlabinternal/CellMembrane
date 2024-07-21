#' @title Alternative ssGSEA implementation for Seurat objects
#'
#' @description This function will run a modified version of ssGSEA. 
#' @param seuratObj A Seurat object.
#' @param geneSets A list of gene sets to run ssGSEA on. The gene sets should be named lists of gene symbols.
#' @param method The method to use for ssGSEA. The dts, wass, and ks methods from the twosamples package are supported.
#' @param groupSize The number of groups to split the data into for parallel processing.
#' @param assay The assay to use.
#' @param layer The layer to use.
#' @param outputAssayName The name of the assay to store the results in.
#' @param parallelWorkers The number of parallel workers/cores to use.
#' @param timing If TRUE, the time taken to run the function will be printed.
#' @examples 
#' \dontrun{
#' # Load the Seurat object
#' data("pbmc_small")
#' seuratObj <- pbmc_small
#' # Get geneSets
#' geneSets <- escape::getGeneSets(library = "H")
#' # Run the alternative ssGSEA
#' seuratObj <- AlternativeSsgseaSeurat(seuratObj, geneSets = geneSets)
#' }
#' @export
#' @return A Seurat object with the ssGSEA scores stored in the specified assay.


AlternativeSsgseaSeurat <- function(seuratObj = seuratObj, 
                                     geneSets, 
                                     method = "dts", 
                                     groupSize = 20000, 
                                     assay = "RNA", 
                                     layer = 'counts', 
                                     outputAssayName = "ssGSEA.alternative", 
                                     parallelWorkers = 4, 
                                     timing = FALSE) {
  # basic checks for the validity of the Seurat object
  if (is.null(seuratObj)) {
    stop('Error: seuratObj not found. Please define a Seurat Object containing a count matrix.')
  } else if (!(layer %in% SeuratObject::Layers(seuratObj))){
    stop('Error: Layer ', layer, ' not found in Seurat Object. Please check the layer name.')
  } else if (!(assay %in% SeuratObject::Assays(seuratObj))){
    stop('Error: Assay ', assay, ' not found in Seurat Object. Please check the assay name.')
  }
  
  # Check that the gene sets are named lists of vectors
  if (!is.list(geneSets)) {
    stop('Error: geneSets must be a named list of gene sets.')
  } else if (!all(sapply(geneSets, is.vector))) {
    stop('Error: Non-vector element within geneSets. Each entry within the list must be a vector of gene symbols.')
  }
  
  # harvest the count matrix
  countMatrix <- SeuratObject::GetAssayData(seuratObj, layer = layer, assay = assay)
  # create a matrix to store the enrichment scores
  scores <- matrix(0, nrow = length(geneSets), ncol = ncol(countMatrix))
  rownames(scores) <- names(geneSets)
  colnames(scores) <- colnames(countMatrix)
  
  # set up the parallel processing plan & check user's number of cpu cores.
  if (future::availableCores() == 1) {
    future::plan("sequential", workers = 1 )
  } else {
    if (future::availableCores()-1 > parallelWorkers) {
      message("Further parallelization is possible. If memory pressure/limits permit, consider increasing the parallelWorkers parameter up to a maximum of " , 
      future::availableCores()-1, ". Current parallelWorkers value is ", parallelWorkers, ".")
    }
    future::plan("multisession", workers = min(parallelWorkers, future::availableCores()-1))
  }
  
  if (timing) {
    start_time <- Sys.time()
  }
  # Create a progress handler to print a progress bar
  progressr::handlers(list(
    progressr::handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 75,
      complete = "+"
    )))
  
  progressr::with_progress({
    # Calculate the number of groups to split the data into for parallel processing. If you supply a negative number, it will default to 1 group.
    numberOfGroups <- max(ceiling(ncol(countMatrix)/groupSize),1)
    # Set the total number of progress steps
    p <- progressr::progressor(steps = length(geneSets) * numberOfGroups)
    subsetScoresList <- list()
    for (group in seq_along(1:numberOfGroups)) {
      # Subset the count matrix to a managable size for parallelization
      # the first index is a max(1, ...) to prevent a zero index
      # the second index is a min(..., ncol(countMatrix)) to prevent an index out of bounds error
      subsetMatrix <- countMatrix[, (max(1,(group-1)*groupSize+1)):min(group*groupSize, ncol(countMatrix))]
      # Run ssGSEA for each gene set
      individualGroupScoresList <- future.apply::future_lapply(names(geneSets), function(gene_set_name) {
        # Get the gene set & enforce that the gene sets exist
        gene_set <- unlist(geneSets[gene_set_name])
        gene_set <- gene_set[gene_set %in% rownames(subsetMatrix)]
        # 
        tmpScores <- .computeSsgseaScores(gene_set, subsetMatrix, method)
        names(tmpScores) <- colnames(subsetMatrix)
        # update the progress bar
        p()
        
        return(tmpScores)
      }, future.seed=1234)
      
      # Convert the list of vectors to a matrix
      individualGroupScores <- do.call(rbind, individualGroupScoresList)
      rownames(individualGroupScores) <- names(geneSets)
      #store the subset scores in a list to be combined after iterating through all groups
      subsetScoresList[[group]] <- individualGroupScores
    }
    # Combine the subset scores into a single matrix
    scores <- do.call(cbind, subsetScoresList)
  })
  if (timing) {
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  seuratObj[[outputAssayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(scores))
  return(seuratObj)
}

#' @title Compute ssGSEA scores for a gene set
#' @description This function computes the ssGSEA scores for a gene set using the twosamples package.
#' @param gene_set A vector of gene symbols to compute the ssGSEA scores for.
#' @param countMatrix A matrix of gene expression values where rows are genes and columns are samples.
#' @param method The method to use for ssGSEA. The dts, wass, and ks methods from the twosamples package are supported.
#' @return A vector of ssGSEA scores from a gene set for each sample.
#' 
.computeSsgseaScores <- function(gene_set, countMatrix, method) {
  distances <- sapply(colnames(countMatrix), function(sample) {
  ranked_genes <- rank(-countMatrix[, sample])
  genes_in_set <- gene_set
  
  # Subset the ranked genes for the gene set and all genes
  ranks_all <- sort(ranked_genes)
  ranks_set <- sort(ranked_genes[genes_in_set])
  ecdf_all <- ecdf(ranks_all)
  ecdf_set <- ecdf(ranks_set)

  ecdf_values_all <- unlist(lapply(seq(1:(nrow(countMatrix))), function(i) ecdf_all(ranks_all[i])))
  ecdf_values_set <- unlist(lapply(seq(1:(nrow(countMatrix))), function(i) ecdf_set(ranks_all[i])))
  
  # Compute a distribution distance using the twosamples package
  if (method == 'dts') {
    distance <- twosamples::dts_stat(ecdf_values_set, ecdf_values_all)
  } else if (method %in% c('wass', "wasserstein")) {
    distance <- twosamples::wass_stat(ecdf_values_all,ecdf_values_set)
  } else if (method == 'ks') {
    distance <- twosamples::ks_stat(ecdf_values_set, ecdf_values_all)
  } else {
    stop('Method: ', method, ' not supported')
  }
  
  # Return the dts distance as the enrichment score
  distance
})
  return(distances)
}
