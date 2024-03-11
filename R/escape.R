#' @title Run escape ssGSEA
#'
#' @description This will run escape to calcualte ssGSEA on Hallmark gene sets
#' @param seuratObj A Seurat object.
#' @param outputAssayName The name of the assay to store results
#' @param doPlot If true, a FeaturePlot will be printed for each pathway
#' @param geneSets Which gene sets should be fetched from MSigDB and calculated. Some recommendations in increasing computation time: H (hallmark, 50 gene sets), C8 (scRNASeq cell type markers, 830 gene sets), C2 (curated pathways, 6366 gene sets), GO:BP (GO biological processes, 7658). 
#' @return The seurat object with results stored in an assay
#' @export
RunEscape <- function(seuratObj, outputAssayName = "escape.ssGSEA", doPlot = FALSE, geneSets = c("H")) {
  
  ##GO:BP and other hierarchical gene sets require subcategories to be passed with them, so we should parse those individually and concatenate afterwards. I think each of these need to be individually supported, but I haven't confirmed. 
  #check if all gene sets are non-hierarchical
  if (all(geneSets %in% c("H", paste0("C", 1:8)))) {
    GS <- escape::getGeneSets(library = geneSets)
  } else if ("GO:BP" %in% geneSets) {
    geneSets <- geneSets[geneSets != "GO:BP"]
    GS_GO_BP <- escape::getGeneSets(library = "C5", subcategory = "BP")
    GS <- c(escape::getGeneSets(library = geneSets), GS_GO_BP)
  } else {
    stop("geneSets is not one of: H, GO:BP",paste0(", C", 1:8), ". Please ensure your requested geneSet is listed, or contact members of the Bimber Lab to add a gene set.")
  }
  
  seuratObj <- escape::runEscape(seuratObj,
                             method = "ssGSEA",
                             gene.sets = GS,
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
      print(suppressWarnings(Seurat::FeaturePlot(seuratObj, features = fn, min.cutoff = 'q02', max.cutoff = 'q98')))
    }
  }

  return(seuratObj)
}


