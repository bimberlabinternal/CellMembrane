#' @title Run CS-CORE
#'
#' @description This will run CS-CORE on the default assay
#' @param seuratObj A Seurat object.
#' @param genes_selected Passed directly to CSCORE::CSCORE()
#' @param saveFile An optional file to save an RDS of the CS-CORE results
#' @return The results after WGCNA
#' @export
RunCsCore <- function(seuratObj, genes_selected = Seurat::VariableFeatures(seuratObj), saveFile = NULL) {
  CSCORE_result <- CSCORE::CSCORE(seuratObj, genes = genes_selected)

  # Obtain CS-CORE co-expression estimates
  CSCORE_coexp <- CSCORE_result$est

  # Obtain BH-adjusted p values
  CSCORE_p <- CSCORE_result$p_value
  p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
  p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
  p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

  # Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
  CSCORE_coexp[p_matrix_BH > 0.05] <- 0

  # Compute the adjacency matrix based on the co-expression matrix
  adj <- WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 1)

  # Compute the topological overlap matrix
  TOM <- WGCNA::TOMsimilarity(adj)
  dissTOM <- 1 - TOM
  rownames(dissTOM) <- colnames(dissTOM) <- genes_selected

  # Run hierarchical clustering as in the WGCNA workflow
  hclust_dist <- stats::hclust(stats::as.dist(dissTOM), method = "average")
  memb <- dynamicTreeCut::cutreeDynamic(dendro = hclust_dist,
                                       distM = dissTOM,
                                       deepSplit = 2,
                                       pamRespectsDendro = FALSE,
                                       minClusterSize = 10)
  # For more instructions on how to tune the parameters in the WGCNA workflow,
  # please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

  names(memb) <- genes_selected
  module_list <- lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))

  if (!is.null(saveFile)) {
    saveRDS(CSCORE_result, file = saveFile)
  }

  return(module_list)
}