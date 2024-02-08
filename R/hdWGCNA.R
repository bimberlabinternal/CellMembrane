#' @title Run hdWGCNA
#'
#' @description Runs hdWGCNA on a seurat object
#' @param seuratObj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment
#' @param groupName The name of the group of interest in the groupBy column
#' @param groupBy The metadata column containing the cell type info. 
#' @param geneSelectionApproach The gene selection approach
#' @param fraction Fraction of cells that a gene needs to be expressed in order to be included
#' @param sampleGroupingVariables Metadata variables on which to group cells for MetacellsByGroups()
#' @param reductionName Select the dimensionality reduction to perform KNN on
#' @import hdWGCNA
#' @return The modified seurat object
#' @export
RunHdWGCNA <- function(seuratObj, wgcna_name, groupName, groupBy, sampleGroupingVariables = c('SubjectId', 'CellType'), geneSelectionApproach = 'fraction', fraction = 0.05, reductionName = 'harmony', assayName = 'RNA') {
  seuratObj <- SetupForWGCNA(
    seuratObj,
    gene_select = geneSelectionApproach,
    fraction = 0.05,
    wgcna_name = wgcna_name
  )

  seuratObj <- MetacellsByGroups(
    seurat_obj = seuratObj,
    group.by = sampleGroupingVariables,
    reduction = reductionName, 
    k = 25, # nearest-neighbors parameter
    max_shared = 10, 
    ident.group = groupBy 
  )

  # normalize metacell expression matrix:
  seuratObj <- NormalizeMetacells(seuratObj, verbose = FALSE)

  seuratObj <- SetDatExpr(
    seuratObj,
    group_name = groupName,
    group.by = groupBy,
    assay = assayName,
    slot = 'data' 
  )

  seuratObj <- TestSoftPowers(
    seuratObj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )

  # plot the results:
  print(patchwork::wrap_plots(PlotSoftPowers(seuratObj), ncol=2))

  power_table <- GetPowerTable(seuratObj)
  head(power_table)

  seuratObj <- ConstructNetwork(
    seuratObj, 
    soft_power=9,
    setDatExpr=FALSE,
    tom_name = 'INH' # name of the topoligical overlap matrix written to disk
  )

  PlotDendrogram(seuratObj, main='INH hdWGCNA Dendrogram')

  # need to run ScaleData first or else harmony throws an error:
  seuratObj <- ScaleData(seuratObj, features=VariableFeatures(seuratObj))

  # compute all MEs in the full single-cell dataset
  seuratObj <- ModuleEigengenes(
    seuratObj,
    group.by.vars="Sample"
  )

  # harmonized module eigengenes:
  hMEs <- GetMEs(seuratObj)

  # module eigengenes:
  MEs <- GetMEs(seuratObj, harmonized=FALSE)

  # compute eigengene-based connectivity (kME):
  seuratObj <- ModuleConnectivity(
    seuratObj,
    group.by = 'cell_type', group_name = 'INH'
  )

  # rename the modules
  seuratObj <- ResetModuleNames(
    seuratObj,
    new_name = "INH-M"
  )

  # plot genes ranked by kME for each module
  p <- PlotKMEs(seuratObj, ncol=5)
  print(p)

  seuratObj <- ModuleExprScore(
    seuratObj,
    n_genes = 25,
    method='UCell'
  )

  # make a featureplot of hMEs for each module
  plot_list <- ModuleFeaturePlot(
    seuratObj,
    features='hMEs', # plot the hMEs
    order=TRUE # order so the points with highest hMEs are on top
  )

  # stitch together with patchwork
  wrap_plots(plot_list, ncol=6)

  # plot module correlagram
  ModuleCorrelogram(seuratObj)
}