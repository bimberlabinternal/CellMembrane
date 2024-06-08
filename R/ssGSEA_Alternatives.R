#' @title Alternative ssGSEA implementation for Seurat objects
#'
#' @description This function will run a modified version of ssGSEA. 
#' @param seuratObj A Seurat object.
#' @param geneSets A list of gene sets to run ssGSEA on. The gene sets should be named lists of gene symbols.
#' @param method The method to use for ssGSEA. The dts, wass, and ks methods from the twosamples package are supported.
#' @param groups The number of groups to split the data into for parallel processing.
#' @param assay The assay to use.
#' @param layer The layer to use.
#' @param outputAssayName The name of the assay to store the results in.
#' @examples 
#' \dontrun{
#' # Load the Seurat object
#' data("pbmc_small")
#' seuratObj <- pbmc_small
#' # Get geneSets
#' geneSets <- escape::getGeneSets(library = "H")
#' # Run the alternative ssGSEA
#' seuratObj <- .AlternativeSsgseaSeurat(seuratObj, geneSets = geneSets)
#' }


.AlternativeSsgseaSeurat <- function(seuratObj = seuratObj, 
                                     geneSets, 
                                     method = "dts", 
                                     groups = 5000, 
                                     assay = "RNA", 
                                     layer = 'counts', 
                                     outputAssayName = "ssGSEA.alternative") {
  # harvest the count matrix
  countMatrix <- SeuratObject::GetAssayData(seuratObj, layer = layer, assay = assay)
  # create a matrix to store the enrichment scores
  scores <- matrix(0, nrow = length(geneSets), ncol = ncol(countMatrix))
  rownames(scores) <- names(geneSets)
  colnames(scores) <- colnames(countMatrix)
  # set up the parallel processing plan
  future::plan("multisession", workers = parallelly::availableCores())
  # Create a progress handler to print a progress bar
  progressr::handlers("txtprogressbar")
  
  progressr::with_progress({
    # Set the total number of progress steps
    p <- progressr::progressor(steps = length(geneSets))
    # Run ssGSEA for each gene set
    scores_list <- future.apply::future_lapply(names(geneSets), function(gene_set_name) {
      # Get the gene set & enforce that the gene sets exist
      gene_set <- unlist(geneSets[gene_set_name])
      gene_set <- gene_set[gene_set %in% rownames(countMatrix)]
      # 
      scores <- sapply(colnames(countMatrix), function(sample) {
        ranked_genes <- rank(-countMatrix[, sample])
        genes_in_set <- gene_set
        
        # Subset the ranked genes for the gene set and all genes
        ranks_all <- ranked_genes
        ranks_set <- ranked_genes[genes_in_set]
        
        # Compute dts distance using the twosample package
        if (method == 'dts') {
          distance <- twosamples::dts_stat(ranks_set, ranks_all)
        } else if (method %in% c('wass', "wasserstein")) {
          distance <- twosamples::wass_stat(ranks_set, ranks_all)
        } else if (method == 'ks') {
          distance <- twosamples::ks_stat(ranks_set, ranks_all)
        } else {
          stop('Method not supported')
        }
        
        # Return the dts distance as the enrichment score
        distance
      })
      names(scores) <- colnames(countMatrix)
      # update the progress bar
      p()
      
      scores
    }, future.seed=1234)
    
    # Convert the list of vectors to a matrix
    scores <- do.call(rbind, scores_list)
    rownames(scores) <- names(geneSets)
  })
  
  seuratObj[[outputAssayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(scores))
  return(seuratObj)
}
