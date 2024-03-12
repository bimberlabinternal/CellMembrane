#' @title Run escape ssGSEA
#'
#' @description This will run escape to calcualte ssGSEA on Hallmark gene sets
#' @param seuratObj A Seurat object.
#' @param outputAssayName The name of the assay to store results
#' @param doPlot If true, a FeaturePlot will be printed for each pathway
#' @param msigdbGeneSets A vector containing gene set codes specifying which gene sets should be fetched from MSigDB and calculated. Some recommendations in increasing computation time: H (hallmark, 50 gene sets), C8 (scRNASeq cell type markers, 830 gene sets), C2 (curated pathways, 6366 gene sets), GO:BP (GO biological processes, 7658). 
#' @param customGeneSets A (preferably named) list containing gene sets to be scored by escape. 
#' @return The seurat object with results stored in an assay
#' @export
RunEscape <- function(seuratObj, outputAssayName = "escape.ssGSEA", doPlot = FALSE, msigdbGeneSets = c("H"), customGeneSets = NULL) {
  #gene set vector to be populated by msigdb or custom gene sets
  GS <- c()
  
  #ensure msigdbGeneSets is formatted properly to parse. 
  if (!is.null(msigdbGeneSets) & !is.na(msigdbGeneSets) & !(length(msigdbGeneSets) == 0)) {
    #require vector. customGeneSets needs to be a list, so this could be a point of confusion
    if (!is.vector(msigdbGeneSets)) {
      stop("msigdbGeneSets is not a vector. Please coerce it to a vector of supported characters.")
    }
    ##GO:BP and other hierarchical gene sets require subcategories to be passed with them, so we should parse those individually and concatenate afterwards. I think each of these need to be individually supported, but I haven't confirmed. 
    #check if all gene sets are non-hierarchical
    if (all(msigdbGeneSets %in% c("H", paste0("C", 1:8)))) {
      #fetch non-hierarchical gene sets
      GS <- c(GS, escape::getGeneSets(library = msigdbGeneSets))
    } else if ("GO:BP" %in% msigdbGeneSets & all(msigdbGeneSets[msigdbGeneSets != "GO:BP"] %in% c("H", paste0("C", 1:8)))) {
      #remove GO:BP from the list and fetch hierarchical gene set.
      msigdbGeneSets <- msigdbGeneSets[msigdbGeneSets != "GO:BP"]
      GS_GO_BP <- escape::getGeneSets(library = "C5", subcategory = "BP")
      #fetch 
      GS <- c(GS, c(escape::getGeneSets(library = msigdbGeneSets), GS_GO_BP))
    } else {
      #if the msigdb gene set is hierarchical but unsupported, throw error: 
      unsupportedGeneSets <- msigdbGeneSets[!(msigdbGeneSets %in% c("GO:BP", "H", paste0("C", 1:8)))]
      stop(paste0(unsupportedGeneSets, " in msigdbGeneSets are unsupported. Please ensure msigdbGeneSets is any of: H, GO:BP, ", paste0(", C", 1:8), ". Please ensure your requested geneSet is listed, add the gene set to the customGeneSets argument, or contact members of the Bimber Lab to add a gene set."))
    }
  } 
  
  #parse customGeneSets
  if (!is.null(customGeneSets) & !(length(customGeneSets) == 0)) {
    #check typing (msigdbGeneSets is a vector, so this could be a point of confusion).
    if (!is.list(customGeneSets)){
      stop("customGeneSets is not a list. Please coerce it into a named list.")
    }
    #if the names are not provided in the custom gene sets, we should warn the user but we can label them here.
    if (is.null(names(customGeneSets))){
      names(customGeneSets) <- paste0("CustomGeneSet", 1:length(customGeneSets))
      warning("The customGeneSets list is unnamed. Naming the genesets with generic names such as: CustomGeneSet1, CustomGeneSet2 ...")
    }
    #concatenate the customGeneSets to the running list of genes (empty if msigdb is empty)
    GS <- c(GS, customGeneSets)
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


