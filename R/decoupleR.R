#' @title Run decoupleR
#'
#' @description This will run decoupleR to infer TF
#' @param seuratObj A Seurat object.
#' @param sourceAssay The source assay
#' @param targetAssay The target assay to store results
#' @return The seurat object
#' @export
RunDecoupleR <- function(seuratObj, sourceAssay = 'RNA', targetAssay = 'tfsulm', organism = 'human') {
  net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)

  mat <- as.matrix(Seurat::GetAssayData(seuratObj, assay = sourceAssay, layer = 'data'))

  acts <- decoupleR::run_ulm(mat = mat,
                             net = net,
                             .source = 'source',
                             .target = 'target',
                             .mor = 'mor',
                             minsize = 5
  )

  newAssay <- acts %>%
    tidyr::pivot_wider(id_cols = 'source',
                       names_from = 'condition',
                       values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

  newAssay <- Seurat::ScaleData(newAssay)
  seuratObj[[targetAssay]] <- newAssay

  PlotTfData(seuratObj, targetAssay)

  return(seuratObj)
}

PlotTfData <- function(seuratObj, assayName = 'tfsulm') {
  dat <- Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'scale.data')
  if (is.null(dat)) {
    stop(paste0('assay not found: ', assayName))
  }

  n_tfs <- 25

  # Extract activities from object as a long dataframe
  df <- t(as.matrix(dat)) %>%
    as.data.frame() %>%
    dplyr::mutate(cluster = Seurat::Idents(seuratObj)) %>%
    tidyr::pivot_longer(cols = -cluster,
                        names_to = "source",
                        values_to = "score") %>%
    dplyr::group_by(cluster, source) %>%
    dplyr::summarise(mean = mean(score))

  # Get top tfs with more variable means across clusters
  tfs <- df %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(std = stats::sd(mean)) %>%
    dplyr::arrange(-abs(std)) %>%
    head(n_tfs) %>%
    dplyr::pull(source)

  # Subset long data frame to top tfs and transform to wide matrix
  top_acts_mat <- df %>%
    dplyr::filter(source %in% tfs) %>%
    tidyr::pivot_wider(id_cols = 'cluster',
                       names_from = 'source',
                       values_from = 'mean') %>%
    tibble::column_to_rownames('cluster') %>%
    as.matrix()

  # Choose color palette
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)

  my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
                 seq(0.05, 2, length.out = floor(100 / 2)))

  pheatmap::pheatmap(mat = top_acts_mat,
                     color = colors.use,
                     border_color = "white",
                     breaks = my_breaks,
                     cellwidth = 15,
                     cellheight = 15,
                     treeheight_row = 20,
                     treeheight_col = 20)
}