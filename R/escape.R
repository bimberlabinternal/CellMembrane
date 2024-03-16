#' @title Run escape ssGSEA
#'
#' @description This will run escape to calcualte ssGSEA on Hallmark gene sets
#' @param seuratObj A Seurat object.
#' @param outputAssayName The name of the assay to store results
#' @param doPlot If true, a FeaturePlot will be printed for each pathway
#' @param doPca boolean determining if the PCA plot should be computed and printed. 
#' @param msigdbGeneSets A vector containing gene set codes specifying which gene sets should be fetched from MSigDB and calculated. Some recommendations in increasing computation time: H (hallmark, 50 gene sets), C8 (scRNASeq cell type markers, 830 gene sets), C2 (curated pathways, 6366 gene sets), GO:BP (GO biological processes, 7658). 
#' @param customGeneSets A (preferably named) list containing gene sets to be scored by escape.
#' @param assayName The name of the source assay
#' @return The seurat object with results stored in an assay
#' @export
RunEscape <- function(seuratObj, outputAssayName = "escape.ssGSEA", doPlot = FALSE, doPca = TRUE, msigdbGeneSets = c("H"), customGeneSets = NULL, assayName = 'RNA') {
  #gene set vector to be populated by msigdb or custom gene sets
  GS <- c()
  
  #retrieve msigdb gene sets
  if (any(!is.null(msigdbGeneSets), !is.na(msigdbGeneSets))) {
    #GetMsigdbGeneSet sanity checks the input, so we don't need to here.
    GS <- GetMsigdbGeneSet(msigdbGeneSets)
  }
  
  #parse customGeneSets
  if (all(!is.null(customGeneSets), !(length(customGeneSets) == 0))) {
    #check typing (msigdbGeneSets is a vector, so this could be a point of confusion).
    if (!is.list(customGeneSets)){
      stop("customGeneSets is not a list. Please coerce it into a named list.")
    }
    #check for feature existence. 
    #I think it's reasonable to expect someone might grab a new/non-msigdb human gene set and supply it to customGeneSets, so a warning rather than stop for feature existence is sufficient.
    if (!(all(unlist(customGeneSets) %in% rownames(seuratObj)))) {
      warning(paste0("Some genes in customGeneSets are not features in the Seurat object. These features are missing and will not contribute to scoring: ", 
                     paste0(unlist(customGeneSets)[!(unlist(customGeneSets) %in% rownames(seuratObj))], ", ")))
    }
    #if the names are not provided in the custom gene sets, we should warn the user but we can label them here.
    if (is.null(names(customGeneSets))){
      names(customGeneSets) <- paste0("CustomGeneSet", 1:length(customGeneSets))
      warning("The customGeneSets list is unnamed. Naming the genesets with generic names such as: CustomGeneSet1, CustomGeneSet2 ...")
    }

    if (any(grepl(names(customGeneSets), pattern = '_'))) {
      print('Converting underscore to hyphen in set names')
      toFix <- grepl(names(customGeneSets), pattern = '_')
      names(customGeneSets)[toFix] <- gsub(names(customGeneSets)[toFix], pattern = '_', replacement = '-')
    }

    #concatenate the customGeneSets to the running list of genes (empty if msigdb is empty)
    GS <- c(GS, customGeneSets)
  }
  
  seuratObj <- escape::runEscape(seuratObj,
                                 method = "ssGSEA",
                                 gene.sets = GS,
                                 groups = 5000,
                                 min.size = 0,
                                 assay = assayName,
                                 new.assay.name = outputAssayName)
  
  #you may want to only score a couple of gene sets (e.g. msigdbGeneSets = NULL), in which case you need to determine the pca parameters to ensure it will run. 
  if (doPca) {
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
  }
  
  if (doPlot) {
    pathways <- rownames(seuratObj@assays[[outputAssayName]])
    key <- seuratObj@assays[[outputAssayName]]@key
    for (fn in pathways) {
      print(suppressWarnings(Seurat::FeaturePlot(seuratObj, features = paste0(key, fn), min.cutoff = 'q02', max.cutoff = 'q98')))
    }
  }
  
  return(seuratObj)
}

