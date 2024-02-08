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

  # NOTE: the way they store these results as a non-valid assay could be a problem down the line. Consider a different method, such as @misc.
  pathwayData <- seuratObj@assays$METABOLISM$score
  seuratObj@misc[[paste0('METABOLISM.', metabolismType)]] <- pathwayData

  if (!'umap' %in% names(seuratObj@reductions)) {
    print('No UMAP reduction found, skipping plots')
    return(seuratObj)
  }

  # TODO: switch back to built-in function once this is resolved:  https://github.com/wu-yc/scMetabolism/pull/28
  print(.MakeDimPlot(obj = seuratObj, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F))

  input.pathways <- intersect(rownames(pathwayData), c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)"))
  if (length(input.pathways) > 0) {
    for (pv in plotVariables) {
      tryCatch({
        print(scMetabolism::DotPlot.metabolism(obj = seuratObj, pathway = input.pathways, phenotype = pv, norm = "y"))
      }, error = function(e){
        print("Error running scMetabolism::DotPlot.metabolism")
        print(conditionMessage(e))
        traceback()
      })
    }
  }

  return(seuratObj)
}

#' @import ggplot2
.MakeDimPlot <- function(obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1) {
  umap.loc<-obj@reductions$umap@cell.embeddings
  colnames(umap.loc) <- c('UMAP_1', 'UMAP_2')
  row.names(umap.loc) <- colnames(obj)
  signature_exp<-obj@assays$METABOLISM$score
  input.pathway <- pathway
  signature_ggplot<-data.frame(umap.loc, t(signature_exp[pathway,]))

  ggplot(data=signature_ggplot, aes(x=UMAP_1, y=UMAP_2, color = signature_ggplot[,3])) +  #this plot is great
    geom_point(size = size) +
    labs(color = input.pathway) +
    xlab("UMAP 1") +ylab("UMAP 2") +
    theme(aspect.ratio=1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}