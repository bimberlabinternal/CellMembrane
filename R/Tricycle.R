#' @title Run tricycle
#'
#' @description This will run tricycle to calculate cell cycle phase
#' @param seuratObj A Seurat object.
#' @param species The species
#' @param gname.type Passed directly to tricycle functions
#' @return The seurat object with results stored in metadata variables
#' @export
RunTricycle <- function(seuratObj, assay = NULL, species = , gname.type = 'SYMBOL') {
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seuratObj)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = GetAssayData(object = seuratObj, assay = assay, layer = "counts")))
  sce <- scuttle::logNormCounts(sce)

  sce <- tricycle::project_cycle_space(sce, species = 'human', gname.type = gname.type)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding") +
          labs(x = "Projected PC1", y = "Projected PC2") +
          ggtitle(sprintf("Projected cell cycle space (n=%d)", ncol(sce))) +
          theme_bw(base_size = 14)
  )

  sce <- tricycle::estimate_cycle_position(sce)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
          labs(x = "Projected PC1", y = "Projected PC2") +
          ggtitle(sprintf("Projected cell cycle space (n=%d), tricyclePosition", ncol(sce))) +
          theme_bw(base_size = 14)
  )

  sce <- tricycle::estimate_Schwabe_stage(sce, gname.type = gname.type, species = species)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "CCStage") +
          labs(x = "Projected PC1", y = "Projected PC2",
          title = paste0("Projected cell cycle space (n=", ncol(sce), "), CCStage")) +
          theme_bw(base_size = 14)
  )

  print(tricycle::plot_emb_circle_scale(sce, dimred = 1,
      point.size = 3.5, point.alpha = 0.9) +
      theme_bw(base_size = 14)
  )

  toAdd <- data.frame(row.names = rownames(sce@colData), tricyclePosition = sce@colData$tricyclePosition, CCStage = sce@colData$CCStage)
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = toAdd)

  print(Seurat::FeaturePlot(seuratObj, features = 'tricyclePosition'))
  print(Seurat::DimPlot(seuratObj, group.by = 'CCStage'))

  return(seuratObj)
}