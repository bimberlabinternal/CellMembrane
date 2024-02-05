#' @title Run scMetabolism
#'
#' @description Runs scMetabolism on a seurat object
#' @param seuratObj A Seurat object
#' @param method Passed to scMetabolism::sc.metabolism.Seurat
#' @param doImputation Allows users to choose whether impute their data before metabolism scoring.
#' @param metabolismType Either KEGG (85 pathways) or REACTOME (82 pathways)
#' @param plotVariables A list of variables used for grouping when generating a DotPlot
#' @return The modified seurat object
#' @export
RunScMetabolism <- function(seuratObj, method = 'AUCell', doImputation = FALSE, metabolismType = 'KEGG', plotVariables = c('ClusterNames_0.2', 'ClusterNames_0.4', 'ClusterNames_0.6')) {
  seuratObj <- scMetabolism::sc.metabolism.Seurat(obj = seuratObj, method = method, imputation = doImputation, metabolism.type = metabolismType)

  print(scMetabolism::DimPlot.metabolism(obj = seuratObj, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F))

  # NOTE: the way they store these results as a non-valid assay could be a problem down the line. Consider a different method, such as @misc.
  pathwayData <- seuratObj@assays$METABOLISM$score
  seuratObj@misc[[paste0('METABOLISM.', metabolismType)]] <- pathwayData

  input.pathways <- intersect(rownames(pathwayData), c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)"))
  if (length(input.pathway) > 0) {
    for (pv in plotVariables) {
      print(scMetabolism::DotPlot.metabolism(obj = seuratObj, pathway = input.pathways, phenotype = pv, norm = "y"))
    }
  }

  return(seuratObj)
}