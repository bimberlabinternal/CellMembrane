#' @include Utils.R
#' @import Seurat
#' @import ggplot2

utils::globalVariables(
  names = c('FDR', 'gene', 'PValue', 'KeyField', 'TotalCells', 'n_DEG', 'uniqueness', 'regulation'),
  package = 'CellMembrane',
  add = TRUE
)


#' @title Pseudobulk Seurat
#'
#' @description Aggregates raw counts in the seurat object, generating a new seurat object where each same has the sum of the counts, grouped by the desired variables
#' @param seuratObj The seurat object
#' @param groupFields The set of fields on which to group
#' @param assays The assays to aggregate
#' @param additionalFieldsToAggregate An option list of additional fields (which must be numeric). Per field, the mean per group will be computed and stored in the output object.
#' @return An aggregated Seurat object.
#' @export
PseudobulkSeurat <- function(seuratObj, groupFields, assays = NULL, additionalFieldsToAggregate = NULL) {
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

  if (!all(is.null(additionalFieldsToAggregate))) {
    for (fn in additionalFieldsToAggregate) {
      if (!fn %in% colnames(seuratObj@meta.data)) {
        stop(paste0('Missing field: ', fn))
      }

      totals <- as.data.frame(seuratObj@meta.data %>% dplyr::group_by(KeyField) %>% dplyr::summarise(Mean = mean(!!sym(fn))))
      names(totals) <- c('keyField', paste0(fn, '_mean'))
      rownames(totals) <- totals$KeyField
      a <- Seurat::AddMetaData(a, totals[,paste0(fn, '_mean'),drop = FALSE])
    }
  }

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
    tidyr::unite('group', tidyr::all_of(contrast_columns))
  
  colData_intermediate$group <- make.names(colData_intermediate$group)
  
  #apply that group column into the original sce
  seuratObj@meta.data$group <- colData_intermediate$group
  
  # Create the sample level metadata by selecting specific columns
  experiment_information <- data.frame(seuratObj@meta.data,  row.names = NULL) %>%
    dplyr::select(tidyr::all_of(c(sampleIdCol, "group")))
  
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
#' @param logFC_threshold Log fold change threshold pass-through variable to be sent to PerformDifferentialExpression.
#' @return A list with the differential_expression results, volano ggplot object and pvalue_dist ggplot object
#' @export
RunPairwiseContrasts <- function(fit, test.use, logFC_threshold = 1){
  #create a nx2 array of all possible unique pairwise combinations of contrasts
  contrasts <- t(utils::combn(colnames(fit$design), m = 2))
  results <- future.apply::future_lapply(split(contrasts, 1:nrow(contrasts)), future.seed = GetSeed(), FUN = function(x){
    up_contrast<- x[[1]]
    down_contrast <- x[[2]]
    contrast_name <- paste0(up_contrast, "-", down_contrast)
    print(paste("Performing DE for: ", contrast_name))
    contrast <- limma::makeContrasts(contrasts = contrast_name, levels = fit$design)
    result <- PerformDifferentialExpression(fit, contrast, contrast_name, logFC_threshold = logFC_threshold, FDR_threshold = 0.05, test.use = test.use)
  })
  #use contrast names to label the results list
  names(results) <- paste0(contrasts[,1],  "-", contrasts[,2])
  return(results)
}

#' @title CreateStudyWideBarPlot 
#'
#' @description Summarize the results of RunPairwiseContrasts() into a series of bar plots, showing the magnitude of differentially expressed genes over many contrasts.
#' @param pairwise_de_results The result of RunPairwiseContrasts().
#' @param pvalue_threshold The pvalue/FDR threshold used to filter the differential expression tables for significance. 
#' @param logFC_threshold Log fold change threshold used to filter the differential expression tables for significance. 
#' @param LikeVsLike Boolean controlling whether or not one of the contrasts should be held constant in the summarization (i.e., only like cell-types should be considered, since Vaccinated_T cells vs Vaccianted_Myeloid cells will be very different, but not necessarily informative.)
#' @param LikeVsLikeGroupFieldIndex The numerical index of the contrast_columns variable supplied to DesignModelMatrix() that should be held constant for like vs like comparisons.
#' @param facetLikeVsLike Boolean controlling whether or not the final plot should be faceted in a LikeVsLike contrast. If you have several grouping variables, this could be helpful.
#' @return A ggplot2 bar plot showing the magnitudes of the differentially expressed genes.
#' @export

CreateStudyWideBarPlot <- function(pairwise_de_results, pvalue_threshold = 0.05, logFC_threshold = 1, LikeVsLike = F, LikeVsLikeGroupFieldIndex= NULL, facetLikeVsLike = T){
  
  #check to see if the LikeVsLike variables are set properly
  if (!is.logical(LikeVsLike)){
    stop("Please ensure LikeVsLike is equal to either TRUE or FALSE. LikeVsLike determines whether or not you want to only keep comparisons where a specific value of the contrast is fixed between samples. For instance: the contrast Vaccinated_Tcells-Vaccinated_Myeloid could be discarded if you wanted only comparisons between T cells")
  } else if (!LikeVsLike) {
    #if not performing Like vs Like comparisons, force Like vs Like faceting to be false.
    print("LikeVsLike is set to FALSE. Forcing facetLikeVsLike to FALSE as well.")
    facetLikeVsLike <- FALSE
  }
  
  if (!is.logical(facetLikeVsLike)){
    stop("Please ensure facetLikeVsLike is equal to either TRUE or FALSE. In the case of a LikeVsLike comparison, should the unique values of the LikeVsLike contrast be faceted?")
  }
  
  if (LikeVsLike | facetLikeVsLike){
    #determine how many grouping fields are present. This splits on _ and -.  _ determines the number of grouping fields, and - splits the two contrasts. We divide by two because the assumption is that these are all pairwise contrasts.
    number_of_grouping_fields <- length(strsplit(names(pairwise_de_results)[1], split = "_|-")[[1]])/2
    if (is.null(LikeVsLikeGroupFieldIndex)){
      stop("LikeVsLikeGroupFieldIndex is NULL. Please define which groupField should be filtered to perform LikeVsLike comparisons. For instance, if you want to filter for cell type comparisons and your groupFields are c('Tissue', 'Timepoint', 'CellType', 'cDNA_ID'), then you would pass 3 as the LikeVsLikeGroupFieldIndex. Please additionally note that your sampleIdCol will be omitted from this indexing, so for the arragement: c('cDNA_ID', 'Tissue', 'Timepoint', 'CellType'), where sampleIdCol = cDNA_ID, the index should still be specified as 3.")   
    } else if (!is.numeric(LikeVsLikeGroupFieldIndex)){
      stop("LikeVsLikeGroupFieldIndex is non-numeric. Please ensure that LikeVsLikeGroupFieldIndex is the numerical index of the groupField (i.e. 3), rather than the name of the group field. ")
    } else if (LikeVsLikeGroupFieldIndex > number_of_grouping_fields){
      stop(paste0("LikeVsLikeGroupFieldIndex exceeds the total number of grouping fields possible: ", number_of_grouping_fields, ". Please ensure that you omitted the groupField that defines sampleIdCol within the DesignModelMatrix function."))
    } else {
      #If LikeVsLike is properly configured, filter out all of the non-LikeVsLike contrasts
      pairwise_de_results_filtered <- list()
      #mine out the groupingFields info from the contrast and discard contrasts where LikeVsLikeGroupFieldIndex varies.
      for (contrast in names(pairwise_de_results)){
        #isolate the contrasts
        contrasts <- strsplit(contrast, split = "-")
        positive_contrast <- contrasts[[1]][1]
        negative_contrast <- contrasts[[1]][2]
        #mine out the groupingFields info from the contrast
        positive_contrast_LikeVsLikeGroup <- strsplit(positive_contrast, split = "_")[[1]][LikeVsLikeGroupFieldIndex]
        negative_contrast_LikeVsLikeGroup <- strsplit(negative_contrast, split = "_")[[1]][LikeVsLikeGroupFieldIndex]
        #compare the LikeVsLike group and add it to the filtered list if they match
        if (positive_contrast_LikeVsLikeGroup == negative_contrast_LikeVsLikeGroup){
          pairwise_de_results_filtered[contrast] <-  pairwise_de_results[contrast]
        }
      }
      #overwrite the filtered list 
      pairwise_de_results <- pairwise_de_results_filtered
    }
  }
  
  #instantiate vector to store DEG lists during contrast iteration.
  degs <- c()
  
  #Construct Barplot 
  #iterate through contrasts (filtered contrasts if LikeVsLike = T).
  for (contrast in names(pairwise_de_results)){
    #create dataframe to store differentially expressed genes table
    deg <- pairwise_de_results[[contrast]]$differential_expression$table
    title <- pairwise_de_results[[contrast]]$volcano$labels$title
    deg <- tibble::as_tibble(deg)
    deg$contrast <- title
    
    #establish up or down regulation in DE results
    deg$regulation <- NA
    deg[deg$logFC > logFC_threshold & deg$FDR < pvalue_threshold, "regulation"] <- "upregulated"
    deg[deg$logFC < -logFC_threshold & deg$FDR < pvalue_threshold, "regulation"] <- "downregulated"
    
    #filter insignificantly differentially expressed genes
    deg <- deg |> dplyr::filter(!is.na(regulation))
    #concatenate into a cross-contrast list of differentially expressed gene results
    degs <- rbind(degs, deg)
  }

  #iterate through the concatenation of all significant & logFC filtered DEG results for uniqueness. 
  plotting_DEGs <- degs |>
    dplyr::group_by(gene) |>
    dplyr::mutate(gene_occurance = dplyr::n()) 
  
  #determine regulation-uniqueness pairs
  plotting_DEGs[plotting_DEGs$regulation == "upregulated","uniqueness"] <- "up_nonunique"
  plotting_DEGs[plotting_DEGs$regulation == "upregulated" & plotting_DEGs$gene_occurance == 1, "uniqueness"] <- "up_unique"
  plotting_DEGs[plotting_DEGs$regulation == "downregulated","uniqueness"] <- "down_nonunique"
  plotting_DEGs[plotting_DEGs$regulation == "downregulated" & plotting_DEGs$gene_occurance == 1, "uniqueness"] <- "down_unique"
  
  #geom_bar sums the n_DEG column, so fixing the sign is all that is necessary
  plotting_DEGs[plotting_DEGs$regulation == "upregulated", "n_DEG"] <- 1 
  plotting_DEGs[plotting_DEGs$regulation == "downregulated", "n_DEG"] <- -1 
  
  if (!facetLikeVsLike){
    #construct plot without faceting
    DEG_bargraph <- ggplot2::ggplot(plotting_DEGs) + 
      ggplot2::geom_bar(data = plotting_DEGs, ggplot2::aes(x = contrast, y = n_DEG, fill = uniqueness), position="stack", stat="identity") + 
      ggplot2::scale_fill_manual(values = c("cadetblue2", "blue", "orange", "red")) + 
      ggplot2::labs(fill="Unique") + 
      ggplot2::ylab("Number of DEGs")+ 
      egg::theme_article() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  } else {
    #mine the faceting grouping variables from the contrast names using LikeVsLikeGroupFieldIndex
    contrast_splits <- strsplit(plotting_DEGs$contrast, split = "-")
    for (contrast_split_index in 1:length(contrast_splits)){
      contrast_splits[[contrast_split_index]]
      plotting_DEGs[contrast_split_index,'faceting_variable_values'] <- strsplit(contrast_splits[[contrast_split_index]][1], split = "_")[[1]][LikeVsLikeGroupFieldIndex]
    }
    #Construct the faceted plot
    DEG_bargraph <- ggplot2::ggplot(plotting_DEGs) + 
      ggplot2::geom_bar(data = plotting_DEGs, ggplot2::aes(x = contrast, y = n_DEG, fill = uniqueness), position="stack", stat="identity") + 
      ggplot2::scale_fill_manual(values = c("cadetblue2", "blue", "orange", "red")) + 
      ggplot2::labs(fill="Unique") + 
      ggplot2::ylab("Number of DEGs")+ 
      egg::theme_article() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + 
      facet_wrap(~faceting_variable_values, scales = "free_x")
  }
  
  return(DEG_bargraph)
}
