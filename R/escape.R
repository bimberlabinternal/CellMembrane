#' @title Run escape ssGSEA
#'
#' @description This will run escape to calcualte ssGSEA on Hallmark gene sets
#' @param seuratObj A Seurat object.
#' @param outputAssayName The name of the assay to store results
#' @param doPlot If true, a FeaturePlot will be printed for each pathway
#' @return The seurat object with results stored in an assay
#' @export
RunEscape <- function(seuratObj, outputAssayName = "escape.ssGSEA", doPlot = FALSE) {
  GS.hallmark <- escape::getGeneSets(library = "H")

  seuratObj <- escape::runEscape(seuratObj,
                             method = "ssGSEA",
                             gene.sets = GS.hallmark,
                             groups = 5000,
                             min.size = 0,
                             new.assay.name = outputAssayName)

  seuratObj <- escape::performPCA(seuratObj,
                              assay = outputAssayName,
                              n.dim = 1:10)

  print(escape::pcaEnrichment(seuratObj,
                            dimRed = "escape.PCA",
                            x.axis = "PC1",
                            y.axis = "PC2",
                            add.percent.contribution = TRUE,
                            display.factors = TRUE,
                            number.of.factors = 10
  ))

  if (doPlot) {
    pathways <- rownames(seuratObj@assays[[outputAssayName]])
    for (fn in pathways) {
      print(suppressWarnings(Seurat::FeaturePlot(so, features = fn, min.cutoff = 'q02', max.cutoff = 'q98')))
    }
  }

  return(seuratObj)
}


