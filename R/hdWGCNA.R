#' @title Run hdWGCNA
#'
#' @description Runs hdWGCNA on a seurat object
#' @param seuratObj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment
#' @param groupName The name of the group of interest in the groupBy column
#' @param groupBy The metadata column containing the cell type info. 
#' @param moduleConnectivityGroupBy Passed directly to ModuleConnectivity() description
#' @param geneSelectionApproach The gene selection approach
#' @param fraction Fraction of cells that a gene needs to be expressed in order to be included
#' @param sampleGroupingVariables Metadata variables on which to group cells for MetacellsByGroups()
#' @param reductionName Select the dimensionality reduction to perform KNN on
#' @param soft_power Passed directly to hdWGCNA::ConstructNetwork
#' @return The modified seurat object
#' @export
RunHdWGCNA <- function(seuratObj, wgcna_name, groupName, groupBy, moduleConnectivityGroupBy, sampleGroupingVariables = c('SubjectId', 'CellType'), geneSelectionApproach = 'fraction', fraction = 0.05, reductionName = 'harmony', assayName = 'RNA', soft_power = NULL) {
  seuratObj <- hdWGCNA::SetupForWGCNA(
    seuratObj,
    gene_select = geneSelectionApproach,
    fraction = fraction,
    wgcna_name = wgcna_name
  )

  seuratObj <- hdWGCNA::MetacellsByGroups(
    seurat_obj = seuratObj,
    group.by = sampleGroupingVariables,
    reduction = reductionName, 
    k = 25, # nearest-neighbors parameter
    max_shared = 10, 
    ident.group = groupBy 
  )

  # normalize metacell expression matrix:
  seuratObj <- hdWGCNA::NormalizeMetacells(seuratObj, verbose = FALSE)

  seuratObj <- hdWGCNA::SetDatExpr(
    seuratObj,
    group_name = groupName,
    group.by = groupBy,
    assay = assayName,
    slot = 'data' 
  )

  seuratObj <- hdWGCNA::TestSoftPowers(
    seuratObj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )

  # plot the results:
  print(patchwork::wrap_plots(hdWGCNA::PlotSoftPowers(seuratObj), ncol=2))

  seuratObj <- hdWGCNA::ConstructNetwork(
    seuratObj, 
    soft_power = soft_power,
    setDatExpr = FALSE,
    overwrite_tom = TRUE,
    tom_name = groupName
  )

  print(hdWGCNA::PlotDendrogram(seuratObj, main = paste0(groupName, ' hdWGCNA Dendrogram')))

  # need to run ScaleData first or else harmony throws an error:
  seuratObj <- Seurat::ScaleData(seuratObj, features = Seurat::VariableFeatures(seuratObj), verbose = FALSE)

  # compute all MEs in the full single-cell dataset
  seuratObj <- hdWGCNA::ModuleEigengenes(
    seuratObj,
    verbose = FALSE,
    group.by.vars = sampleGroupingVariables
  )

  # compute eigengene-based connectivity (kME):
  seuratObj <- hdWGCNA::ModuleConnectivity(
    seuratObj,
    group.by = moduleConnectivityGroupBy,
    group_name = groupName
  )

  # plot genes ranked by kME for each module
  print(hdWGCNA::PlotKMEs(seuratObj, ncol=5))

  seuratObj <- hdWGCNA::ModuleExprScore(
    seuratObj,
    n_genes = 25,
    method = 'UCell'
  )

  # make a featureplot of hMEs for each module
  plot_list <- hdWGCNA::ModuleFeaturePlot(
    seuratObj,
    features = 'hMEs', # plot the hMEs
    order = TRUE # order so the points with highest hMEs are on top
  )

  # stitch together with patchwork
  print(patchwork::wrap_plots(plot_list, ncol=6))

  # plot module correlagram
  print(hdWGCNA::ModuleCorrelogram(seuratObj))

  unlink('TOM', recursive = TRUE)

  return(seuratObj)
}