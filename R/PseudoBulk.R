#' @include Utils.R
#' @import Seurat
#' @import ggplot2

utils::globalVariables(
  names = c('FDR', 'gene', 'PValue'),
  package = 'CellMembrane',
  add = TRUE
)


#' @title Pseudobulk Seurat
#'
#' @description Aggregates raw counts in the seurat object, generating a new seurat object where each same has the sum of the counts, grouped by the desired variables
#' @param seuratObj The seurat object
#' @param groupFields The set of fields on which to group
#' @param assays The assays to aggregate
#' @return An aggregated Seurat object.
#' @export
PseudobulkSeurat <- function(seuratObj, groupFields, assays = NULL) {
  if (!all(groupFields %in% names(seuratObj@meta.data))) {
    stop('All fields from groupFields must be in seuratObj@meta.data')
  }

  # TODO: perhaps filtering on saturation, min.counts or other features??
  seuratObj$KeyField <- apply(seuratObj@meta.data[,groupFields,drop = FALSE], 1, function(y){
    return(paste0(y, collapse = '_'))
  })

  Seurat::Idents(seuratObj) <- seuratObj$KeyField

  # This generates the mean() of counts. Even though we want sum(), this is a convenient way to ensure all other
  a <- Seurat::AverageExpression(seuratObj, return.seurat = T, verbose = F, slot = 'counts', assays = assays)

  metaGrouped <- unique(seuratObj@meta.data[,c('KeyField', groupFields),drop = FALSE])
  rownames(metaGrouped) <- metaGrouped$KeyField
  metaGrouped <- metaGrouped[,names(metaGrouped) != 'KeyField',drop = FALSE]
  a <- Seurat::AddMetaData(a, metaGrouped)

  totals <- as.data.frame(seuratObj@meta.data %>% dplyr::group_by(KeyField) %>% dplyr::summarise(TotalCells = n()))
  rownames(totals) <- totals$KeyField

  a <- Seurat::AddMetaData(a, totals[,'TotalCells',drop = FALSE])

  print(ggplot(a@meta.data, aes(x = TotalCells)) +
    geom_density() +
    egg::theme_presentation(base_size = 18) +
    labs(x = 'Cells/Sample', y = '# Cells') +
    ggtitle('Total Cells/Sample')
  )

  # Convert mean into sum:
  for (assayName in names(a@assays)) {
    m <- Seurat::GetAssayData(a, assay = assayName, slot = 'counts')
    m2 <- m %*% diag(a$TotalCells)
    rownames(m2) <- rownames(m)
    colnames(m2) <- colnames(m)
    a <- Seurat::SetAssayData(a, assay = assayName, slot = 'counts', new.data = Seurat::as.sparse(m2))
    a <- Seurat::NormalizeData(a, verbose = FALSE, assay = assayName)
  }

  return(a)
}


#' @title DesignModelMatrix
#'
#' @description Creates an edgeR glm object
#' @param seuratObj The seurat object
#' @param contrast_columns A vector of columns to contrast
#' @param sampleIdCol An additional column denoting the variable containing the sample (for grouping)
#' @return An edgeR glm object
#' @export
DesignModelMatrix <- function(seuratObj, contrast_columns, sampleIdCol = "cDNA_ID"){
  #Create a dummy sce@colData that unites the contrast columns into a single "group" column
  #combined columns are (by default) separated by an underscore
  colData_intermediate <- seuratObj@meta.data |>
    as.data.frame() |>
    tidyr::unite('group', contrast_columns)
  
  colData_intermediate$group <- make.names(colData_intermediate$group)
  
  #apply that group column into the original sce
  seuratObj@meta.data$group <- colData_intermediate$group
  
  # Create the sample level metadata by selecting specific columns
  experiment_information <- data.frame(seuratObj@meta.data,  row.names = NULL) %>%
    dplyr::select(sampleIdCol, "group")
  
  design <- stats::model.matrix(~ 0 + experiment_information$group) %>%
    magrittr::set_rownames(experiment_information[[sampleIdCol]]) %>%
    magrittr::set_colnames(levels(factor(experiment_information$group)))

  return(design)
}

#' @title PerformGlmFit
#'
#' @description Creates an edgeR glm object
#' @param seuratObj The seurat object
#' @param design The model.matrix object
#' @param test.use Can be either QLF or LRT. QLF runs edgeR::glmQLFTest, while LRT runs edgeR::glmLRT
#' @param assayName The name of the assay to use
#' @param minCountsPerGene Any genes with fewer than this many counts (across samples) will be dropped
#' @return An edgeR glm object
#' @export
PerformGlmFit <- function(seuratObj, design, test.use = "QLF", assayName = 'RNA', minCountsPerGene = 1){
  #convert seurat object to SingleCellExperiment for edgeR
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = seuratObj@assays[[assayName]]@counts), colData = seuratObj@meta.data)

  #filter out lowly expressed genes
  if (!is.null(minCountsPerGene)) {
    sce <- sce[rowSums(as.matrix(SingleCellExperiment::counts(sce))) > minCountsPerGene, ]
  }

  #preprocess the DGEList and fit the model
  y <- edgeR::DGEList(SingleCellExperiment::counts(sce), remove.zeros = TRUE)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  print(edgeR::plotBCV(y))

  if (test.use == "QLF"){
    fit <- edgeR::glmQLFit(y, design)
  } else if (test.use == "LRT"){
    fit <- edgeR::glmFit(y, design)
  } else {
    stop("Please supply a valid test.use argument.")
  }

  return(fit)
}

#' @title PerformDifferentialExpression
#'
#' @description This is designed to perform DE using edgeR on a pseudobulked seurat object
#' @param fit The glm fit object generated by PerformGlmFit
#' @param contrast The contrast, passed directly to edgeR::glmQLFTest or edgeR::glmLRT
#' @param contrast_name A name for the contrast
#' @param logFC_threshold The min logFC
#' @param FDR_threshold The min FDR to report
#' @param test.use Can be either QLF or LRT. QLF runs edgeR::glmQLFTest, while LRT runs edgeR::glmLRT
#' @return A list with the differential_expression results, volano ggplot object and pvalue_dist ggplot object
#' @export
PerformDifferentialExpression <- function(fit, contrast, contrast_name, logFC_threshold = 1, FDR_threshold = 0.05, test.use = "QLF"){
  #perform differential expression for the given contrast
  if (test.use == "QLF"){
    fit_test <- edgeR::glmQLFTest(fit, contrast = contrast)
  } else if (test.use == "LRT"){
    fit_test <- edgeR::glmLRT(fit, contrast = contrast)
  } else {
    stop("Please supply a valid test.use argument.")
  }

  #distill results into differential_expression list
  differential_expression <- edgeR::topTags(fit_test, n= Inf)
  differential_expression$table$gene <- rownames(differential_expression$table)

  #create
  volcano <- ggplot(differential_expression$table, aes(x= logFC, y = -log10(FDR), label = gene)) +
    #plot the upregulated genes using logFC cutoff and FDR cutoff
    geom_point(data = differential_expression$table[differential_expression$table$logFC >=logFC_threshold &
                                                      differential_expression$table$FDR <=FDR_threshold,], color = "orange")  +
    #plot the downregulated genes
    geom_point(data = differential_expression$table[differential_expression$table$logFC <=-logFC_threshold &
                                                      differential_expression$table$FDR <=FDR_threshold,], color = "cadetblue2") +
    #repel labels for the up- and down-regulated genes
    ggrepel::geom_text_repel(data = differential_expression$table[
        (differential_expression$table$logFC >= logFC_threshold & differential_expression$table$FDR < FDR_threshold) |
        (differential_expression$table$logFC <= -logFC_threshold & differential_expression$table$FDR < FDR_threshold),] ,
        max.overlaps = 10) +
    #plot the rest of the genes
    geom_point(data = differential_expression$table[(differential_expression$table$logFC > -logFC_threshold & differential_expression$table$logFC < logFC_threshold) |differential_expression$table$FDR > FDR_threshold,], color = "black") +
    #plot thresholds
    geom_hline(yintercept = -log10(FDR_threshold), color = "red", linetype="dashed")+
    geom_vline(xintercept = logFC_threshold, color = "red", linetype="dashed") +
    geom_vline(xintercept = -logFC_threshold, color = "red", linetype="dashed") +
    egg::theme_article() +
    ggtitle(contrast_name)

  print(volcano)

  pvalue_dist <- ggplot(differential_expression$table, aes( x= PValue)) + geom_histogram() + ggtitle("PValue Distribution")
  print(pvalue_dist)

  return(list(
    differential_expression = differential_expression,
    volcano = volcano,
    pvalue_dist = pvalue_dist)
  )
}

#' @title RunPairwiseContrasts
#'
#' @description A convenience function to iterate and perform pairwise contrasts using the glm object
#' @param fit The glm fit object generated by PerformGlmFit
#' @param test.use Can be either QLF or LRT. QLF runs edgeR::glmQLFTest, while LRT runs edgeR::glmLRT
#' @return A list with the differential_expression results, volano ggplot object and pvalue_dist ggplot object
#' @export
RunPairwiseContrasts <- function(fit, test.use){
  contrasts <- t(utils::combn(colnames(fit$design), m = 2))
  results <- future.apply::future_lapply(split(contrasts, 1:nrow(contrasts)) , FUN = function(x){
    up_contrast<- x[[1]]
    down_contrast <- x[[2]]
    contrast_name <- paste0(up_contrast, "-", down_contrast)
    print(paste("Performing DE for: ", contrast_name))
    contrast <- limma::makeContrasts(contrasts = contrast_name, levels = fit$design)
    result <- PerformDifferentialExpression(fit, contrast, contrast_name, logFC_threshold = 1, FDR_threshold = 0.05, test.use = test.use)
  })

  return(results)
}