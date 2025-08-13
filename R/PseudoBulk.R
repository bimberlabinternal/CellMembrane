#' @include Utils.R
#' @import Seurat
#' @import ggplot2

utils::globalVariables(
  names = c('FDR', 'gene', 'PValue', 'KeyField', 'TotalCells', 'n_DEG', 'uniqueness', 'regulation', 'contrast_name', 'sampleIdCol', 'joinedFields', 'number_of_positive_nonunique_DEGs', 'number_of_negative_nonunique_DEGs', 'number_of_positive_unique_DEGs', 'number_of_negative_unique_DEGs', "total_DEGs"),
  package = 'CellMembrane',
  add = TRUE
)


#' @title Pseudobulk Seurat
#'
#' @description Aggregates raw counts in the seurat object, generating a new seurat object where each same has the sum of the counts, grouped by the desired variables
#' @param seuratObj The seurat object
#' @param groupFields The set of fields on which to group
#' @param assayToAggregate The assay to aggregate
#' @param additionalFieldsToAggregate An option list of additional fields (which must be numeric). Per field, the mean per group will be computed and stored in the output object.
#' @param metaFieldCollapseCharacter The character to use when concatenating metadata fields together to form the sample key field
#' @param nCountRnaStratification A boolean determining whether or not automatic outlier detection of clusters with abnormal nCount_RNA should be detected.
#' @param stratificationGroupingFields The metadata fields containing groupings of cells to be specified as having abnormal nCount_RNA distributions or not. The "leftmost"/first argument that detects an outlier will break the looping and sequester cells at that value. It is recommended that this vector is arranged by a "least to most" specific hierarchy/progression.
#' @return An aggregated Seurat object.
#' @export
PseudobulkSeurat <- function(seuratObj, 
                             groupFields, 
                             assayToAggregate = Seurat::DefaultAssay(seuratObj),
                             additionalFieldsToAggregate = NULL, 
                             metaFieldCollapseCharacter = '|', 
                             nCountRnaStratification = F, 
                             stratificationGroupingFields = c("ClusterNames_0.2", "ClusterNames_0.4", "ClusterNames_0.6", "ClusterNames_0.8", "ClusterNames_1.2")) {
  if (!all(groupFields %in% names(seuratObj@meta.data))) {
    missing <- groupFields[!groupFields %in% names(seuratObj@meta.data)]
    stop(paste0('All fields from groupFields must be in seuratObj@meta.data. Missing: ', paste0(missing, collapse = ',')))
  }
  
  #QC and pseudobulk based on nCount_RNA distributions
  if (nCountRnaStratification) {
    if (any(grepl("ClusterNames", groupFields)) & any(grepl("ClusterNames", stratificationGroupingFields))){
      warning("It appears that you are pseudobulking on a ClusterNames field, and also specifying outlier removal using one or more ClusterNames fields. Please be aware that this will potentially split a low-resolution cluster groupField according to the nCount_RNA distributions at a higher resolution stratificationGroupingFields. At extreme levels, this can have unintented consequences (i.e. artificially inflating sample sizes) for downstream model fitting.")
    }
    set.seed(GetSeed())
    for (stratificationGroupingField in stratificationGroupingFields) {
      if (!is.integer(as.integer(as.character(seuratObj@meta.data[,stratificationGroupingField])))) {
        stop("Only integer-based groups are capable of being stratified by nCount_RNA. Please specify a ClusterNames-type of grouping field or convert your desired grouping field to an integer.")
      }
      #coerce to a list of dataframes containing nCount_RNA (typically) based on a clustering algorithm
      list_of_cluster_dataframes <- seuratObj@meta.data %>% 
        tibble::rownames_to_column(var = "cellbarcode") %>% 
        dplyr::select(!!sym(stratificationGroupingField), nCount_RNA) %>% 
        dplyr::group_by(!!sym(stratificationGroupingField)) %>% 
        dplyr::group_split()
      
      #determine the minimum cluster size (KL divergence requires equal sampling)
      minimum_cluster_size <- min(unlist(lapply(list_of_cluster_dataframes, FUN = nrow)))
      cluster_nCount_RNA_matrix <- matrix(nrow = minimum_cluster_size)
      #construct matrix of samples from nCount_RNA distributions
      for (cluster_dataframe in list_of_cluster_dataframes) {
        cluster_nCount_RNA_matrix <- cbind(cluster_nCount_RNA_matrix, sample(cluster_dataframe$nCount_RNA, 
                                                                             size = minimum_cluster_size, 
                                                                             replace = F))
      }
      #drop NA column formed when initializing matrix
      cluster_nCount_RNA_matrix <- cluster_nCount_RNA_matrix[,-1]
      #compute KL divergences
      KL_divergences <- flexmix::KLdiv(cluster_nCount_RNA_matrix)
      
      #perform rudimentary outlier detection on column sums of KL divergences assuming a half-normal distribution. 
      abnormal_RNA_clusters_index <- which(colSums(KL_divergences) >= mean(colSums(KL_divergences)) + 2*sd(colSums(KL_divergences)))
      
      #We use 0 as a cluster index, and which() returns integer(0) if no clusters satisfy the check, so set these to NA to prevent accidental filtering.
      #if they're not zero, then we need to zero index them (they're currently indexed as via the column indexing of KL_divergences, which starts at 1)
      if (length(abnormal_RNA_clusters_index) > 0){
        abnormal_RNA_clusters_index <- abnormal_RNA_clusters_index - 1

        metadata <- seuratObj@meta.data %>%
          dplyr::mutate(nCount_RNA_Stratification = dplyr::case_when(!!sym(stratificationGroupingField) %in% abnormal_RNA_clusters_index ~ "AbnormalRNACounts",
                                                                     !(!!sym(stratificationGroupingField) %in% abnormal_RNA_clusters_index) ~ "NormalRNACounts"))
        seuratObj <- Seurat::AddMetaData(seuratObj, metadata = metadata)
        groupFields <- c(groupFields, "nCount_RNA_Stratification")
        print(paste0("Abnormal RNA count distributions detected at grouping field ", stratificationGroupingField, ". Breaking. No further grouping fields will be evaluated."))
        break
      }
    }
  }

  # TODO: perhaps filtering on saturation, min.counts or other features??
  seuratObj$KeyField <- unname(apply(seuratObj@meta.data[,groupFields,drop = FALSE], 1, function(y){
    # NOTE: AverageExpression will convert underscores to hyphens in the sample names anyway, so proactively do this here
    return(gsub(x = paste0(make.names(y), collapse = metaFieldCollapseCharacter), replacement = '-', pattern = '_'))
  }))
  
  Seurat::Idents(seuratObj) <- seuratObj$KeyField
  
  # This generates the sum of counts
  a <- Seurat::AggregateExpression(seuratObj, group.by = 'KeyField', return.seurat = T, verbose = F, assays = assayToAggregate)
  if (class(Seurat::GetAssay(a, assay = assayToAggregate))[1] != 'Assay5') {
    print('Updating assay object to assay5')
    assay5 <- SeuratObject::CreateAssay5Object(counts = Seurat::GetAssayData(a, assay = assayToAggregate, layer = 'counts'))
    a <- Seurat::SetAssayData(a, new.data = assay5, assay = assayToAggregate, layer = 'counts')
  }
  
  metaGrouped <- unique(seuratObj@meta.data[,c('KeyField', groupFields),drop = FALSE])
  rownames(metaGrouped) <- metaGrouped$KeyField
  if (any(sort(metaGrouped$KeyField) != sort(colnames(a)))) {
    sel <- sort(metaGrouped$KeyField) != sort(colnames(a))
    x <- sort(metaGrouped$KeyField[sel])
    y <- sort(colnames(a)[sel])
    
    warning('The keyField and AverageExpression object keys to do not match. Key fields:')
    warning(paste0(x, collapse = ';'))
    warning('Seurat colnames:')
    warning(paste0(y, collapse = ';'))
    stop(paste0('The keyField and AverageExpression object keys to do not match. Metadata: ', paste0(head(x), collapse = ', '), ', matrix: ', paste0(head(y), collapse = ', ')))
  }
  
  metaGrouped <- metaGrouped[,names(metaGrouped) != 'KeyField',drop = FALSE]
  a <- Seurat::AddMetaData(a, metaGrouped)
  
  totals <- as.data.frame(seuratObj@meta.data %>% dplyr::group_by(KeyField) %>% dplyr::summarise(TotalCells = dplyr::n()))
  rownames(totals) <- totals$KeyField
  
  a <- Seurat::AddMetaData(a, totals[,'TotalCells',drop = FALSE])
  
  if (!all(is.null(additionalFieldsToAggregate))) {
    for (fn in additionalFieldsToAggregate) {
      if (!fn %in% colnames(seuratObj@meta.data)) {
        stop(paste0('Missing field: ', fn))
      }
      
      totals <- as.data.frame(seuratObj@meta.data %>% dplyr::group_by(KeyField) %>% dplyr::summarise(Mean = mean(!!sym(fn))))
      names(totals) <- c('KeyField', paste0(fn, '_mean'))
      rownames(totals) <- totals$KeyField
      totals <- totals[,names(totals) != 'KeyField',drop = FALSE]
      a <- Seurat::AddMetaData(a, totals[,paste0(fn, '_mean'),drop = FALSE])
    }
  }
  
  print(ggplot(a@meta.data, aes(x = TotalCells)) +
          geom_density() +
          egg::theme_presentation(base_size = 18) +
          labs(x = 'Cells/Sample', y = '# Cells') +
          ggtitle('Total Cells/Sample')
  )
  
  # Normalize and store pct.expression for each assay:
  for (assayName in names(a@assays)) {
    a <- Seurat::NormalizeData(a, verbose = FALSE, assay = assayName)

    # Makes a new matrix with the expression percentage per gene per KeyField.
    counts <- Seurat::GetAssayData(seuratObj, assay = assayName, layer = "counts")

    percentages <- NULL
    for (keyfield in colnames(a)) {
      if (sum(seuratObj$KeyField == keyfield) == 0) {
        stop(paste0('There are no cells in the seuratObj where KeyField equals: ', keyfield))
      }

      pcts <- counts[,rownames(seuratObj@meta.data)[seuratObj$KeyField == keyfield],drop = FALSE]
      nCells <- ncol(pcts)
      pcts <- apply(pcts, MARGIN = 1, FUN = function(x) {
        return(sum(x > 0))
      })
      
      pcts <- matrix(pcts, nrow = length(pcts))
      pcts <- pcts / nCells
      rownames(pcts) <- rownames(counts)
      colnames(pcts) <- keyfield

      percentages <- cbind(percentages, pcts)
    }

    if (any(colnames(percentages) != colnames(a))) {
      stop('The columns on the pct.expression object do not match the parent seurat object')
    }

    if (any(rownames(percentages) != rownames(counts))) {
      stop('The rows on the pct.expression object do not match the parent seurat object')
    }

    # Adds percentages as a new assay.
    SeuratObject::LayerData(a, assay = assayName, layer = 'pct.expression') <- Seurat::as.sparse(percentages)
  }

  return(a)
}


#TODO: add support to extend this to propeller (pseudobulk the object, then optionally discard.)

#' @title DesignModelMatrix
#'
#' @description This creates a model matrix that groups the samples according to the experimental groups used downstream in statistical analyses.
#' @param seuratObj The seurat object
#' @param contrast_columns A vector of columns to contrast
#' @param sampleIdCol An additional column denoting the variable containing the sample (for grouping)
#' @return A model matrix (samples x groups data frame.)
#' @export
DesignModelMatrix <- function(seuratObj, contrast_columns, sampleIdCol = "cDNA_ID"){
  #Create a dummy sce@colData that unites the contrast columns into a single "group" column
  #combined columns are (by default) separated by an underscore
  
  colData_intermediate <- seuratObj@meta.data |>
    as.data.frame() |>
    select(all_of(contrast_columns)) |>
    mutate(dplyr::across(everything(), \(x) gsub(x,pattern = "_", replacement = "."))) |>
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
  #stash metadata-specific values in the design matrix for retrieval downstream. 
  attr(design, "contrast_columns") <- contrast_columns
  attr(design, "sampleIdCol") <- sampleIdCol
  return(design)
}

#' @title PerformGlmFit
#'
#' @description Creates an edgeR glm object
#' @param seuratObj The seurat object
#' @param design The model.matrix object
#' @param test.use Can be either QLF or LRT. QLF runs edgeR::glmQLFTest, while LRT runs edgeR::glmLRT
#' @param assayName The name of the assay to use
#' @param filterGenes A boolean controlling whether or not to filter genes using edgeR::filterByExpr. If TRUE, genes with low counts will be filtered out.
#' @param legacy A passthrough variable for edgeR's glmQLF function. They recently (R 4.0) changed the default behavior, so this will break on earlier versions of R. 
#' @param plotBCV A boolean determining if the BCV plot should be shown.
#' @return An edgeR glm object
#' @export
PerformGlmFit <- function(seuratObj, design, test.use = "QLF", assayName = 'RNA', filterGenes = TRUE, legacy = FALSE, plotBCV = TRUE){
  #convert seurat object to SingleCellExperiment for edgeR
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts')), colData = seuratObj@meta.data)
  
  #preprocess the DGEList and fit the model
  y <- edgeR::DGEList(SingleCellExperiment::counts(sce), remove.zeros = TRUE)
  if (filterGenes) {
    geneFilter <- edgeR::filterByExpr(y, design = design)
    message(paste0("Number of genes passing filterByExpr filter: ", sum(geneFilter), ". Number of genes with any counts: ", nrow(y), " (", round((sum(geneFilter)/nrow(y)),2) * 100 , "% of detectable genes kept)"))
    y <- y[geneFilter,]
  }
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  if (plotBCV){
    print(edgeR::plotBCV(y))
  }

  
  if (test.use == "QLF"){
    fit <- edgeR::glmQLFit(y, design, legacy = legacy)
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
#' @param showPlots Boolean determining if the volcano plots and p value dsitributions should be shown. 
#' @return A list with the differential_expression results, volcano ggplot object and pvalue_dist ggplot object
#' @export

PerformDifferentialExpression <- function(fit, contrast, contrast_name, logFC_threshold = 1, FDR_threshold = 0.05, test.use = "QLF", showPlots = T){

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
  
  #create volcano plot
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
  
  
  
  pvalue_dist <- ggplot(differential_expression$table, aes(x= PValue)) + 
    geom_histogram() + 
    ggtitle("PValue Distribution") + 
    egg::theme_article()
  
  if (showPlots) {
  print(pvalue_dist)
  print(volcano)
  }
  
  return(list(
    differential_expression = differential_expression,
    volcano = volcano,
    pvalue_dist = pvalue_dist)
  )
}

#' @title FilterPseudobulkContrasts
#'
#' @description This is designed to accept a study design defined by a series of logical gates applied to metadata fields, test them for equivalence to some criterion, and then filter pairwise contrasts accordingly.
#' @param logicList The list that defines the study design. Please see examples for the setup and use. This list be a list of lists where each sub-list has three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor. The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.
#' @param design a design/model matrix returned by DesignModelMatrix().
#' @param useRequireIdenticalLogic Whether or not to apply require_identical logic using requireIdenticalFields. It is possible to use logic gates to keep metadata fields constant, but using this saves time & effort/iteration.
#' @param requireIdenticalFields The metadata columns of the SeuratObj that you wish to keep constant. For example: defining "Tissue" in this vector would filter all contrasts comparing Liver to PBMC, as their tissue is not identical.
#' @param filteredContrastsOutputFile The file to write the list of filtered contrasts to (one pair of samples per row, separated into their metadata fields).
#' @return A dataframe of pairwise contrasts, with column names indicating the directionality (positive vs negative) and associated metadata column (inferred by contrast_columns).
#' @examples
#' \dontrun{
#' #Set up design matrix
#' design <- DesignModelMatrix(pseudobulked_seuratObj, 
#'                             contrast_columns = c('Challenge', 'Tissue', 'cell_type', 'SampleType'), 
#'                             sampleIdCol = 'cDNA_ID')
#' 
#' 
#' #Form the study design. This should be a list of gates (each specified as a list itself) with three elements. 
#' #The first element defines the metadata column or field that will be used for filtering. 
#' #The second defines the logic gate (AND -> any, OR -> any, xor, xnor, nand, nor). 
#' #The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.
#' logicList <- list(
#'                     list('Challenge', 'xor', 'Mock.challenged'), 
#'                     list('SampleType', 'all', 'Necropsy')
#'             )
#'
#' 
#' #Finally, enumerate all possible contrasts, then filter and return them according to the study design defined by logicList.
#' #If you want to perform "Like-vs-Like" contrasts (i.e. T cell vs T cell within a cell_type metadata field), specify those columns as required to be identical and set useRequireIdenticalLogic to TRUE.
#' filtered_contrasts <- FilterPseudobulkContrasts(logicList = logicList, 
#'                                                 design = design, 
#'                                                 useRequireIdenticalLogic = T, 
#'                                                 requireIdenticalFields = c('Tissue','cell_type', 'SampleType'), 
#'                                                 filteredContrastsOutputFile = './filtered_contrasts.tsv')
#'                                                 
#' }
#' @details
#' 
#' | \strong{Filtered} | \emph{Reason}   | \strong{positive_contrast_Challenge} | \emph{positive_contrast_Tissue} | \strong{positive_contrast_cell_type} | \emph{positive_contrast_SampleType} | \strong{negative_contrast_Challenge} | \emph{negative_contrast_Tissue} | \strong{negative_contrast_cell_type} | \emph{negative_contrast_SampleType} |
#' | :-------------: |:-------------:| :----------------------------:| :----------------------------:| :----------------------------:| :----------------------------:| :----------------------------:| :----------------------------:| :----------------------------:| :----------------------------:|
#' | yes   | non-identical cell_types | Mtb | Spleen| \emph{\strong{Myeloid}} | Necropsy | Mock.challenged | Spleen | \emph{\strong{Bcell}} | Necropsy |
#' | | | | | | | | | | | 
#' | | | | | | | | | | | 
#' | no   | N/A      |   Mtb | Lung | T.NK | Necropsy | Mock.challenged | Lung | T.NK | Necropsy |
#' | | | | | | | | | | | 
#' | | | | | | | | | | | 
#' | yes   | fails Challenge xor gate      |   \emph{\strong{Mock.challenged}} | Lung | T.NK | Necropsy | \emph{\strong{Mock.challenged}} | Lung | T.NK | Necropsy |
#' | | | | | | | | | | | 
#' | | | | | | | | | | | 
#' | yes   | fails SampleType AND gate      |   Mock.challenged | MesLN | Bcell | \emph{\strong{Baseline}} | Mock.challenged | MesLN | Bcell | \emph{\strong{Necropsy}} |
#' 
#' @export
#' @md

FilterPseudobulkContrasts <- function(logicList = NULL, design = NULL, useRequireIdenticalLogic = T, requireIdenticalFields = NULL, filteredContrastsOutputFile = './filtered_contrasts.tsv'){
  #check design matrix.
  if (is.null(design)){
    stop("Please define a design matrix. The design matrix is intended to be returned by DesignModelMatrix().")
  }
  #ensure design matrix is from DesignModelMatrix() (or at least has contrast_columns in its attributes.)
  if (is.null(attr(design, "contrast_columns"))){
    stop("'contrast_columns' is missing from the design model matrix attributes. Please use CellMembrane's DesignModelMatrix to create the design matrix, or stash the contrast_columns vector into the design matrix using \n  attr(design, 'contrast_columns') <- contrast_columns ")
  }
  contrast_columns <- attr(design, "contrast_columns")
  #ensure logicList is properly defined.
  if (is.null(logicList)){
    stop("Please supply a list that defines the logic gates used for filtering the contrasts. This list should have three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (typeof(logicList) != "list"){
    stop("Please ensure that logicList is a list composed of lists. Each sub-list should define a gate that will filter possible contrasts. This list should have three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (!(all(lengths(logicList) == 3))){
    stop("The lengths of all of the elements (i.e. logic gates) within logicList are not equal to 3. Please ensure there are exactly three entries in each element of logicList. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (!(all(sapply(logicList, "[[", 1) %in% contrast_columns))){
    stop("There are elements (logic gates) within logicList whose first element is not in the contrast_columns found within the design matrix. contrast_columns is the vector supplied to DesignModelMatrix() upstream of this function and holds the metadata field names of the pseudobulked Seurat object that comprise the contrasts. If you're seeing this, something has gone awry (like using an old design matrix or a typo in the logicList), or you're trying to supply your own design/model matrix, which is currently not supported.")
  } else if (any(sapply(logicList, "[[", 2) %in% c("and", "AND", "or", "OR"))){
    stop("Error in one of the logic gate specifications. For an AND gate, please use the function name 'all'. For an OR gate, please use the function name 'any'.")
  } else if (!all(tolower(sapply(logicList, "[[", 2)) %in% c("any", "all", "xor", "nand", "nor", "xnor"))){
    stop("Error in one of the logic gate specifications. Please use one of: 'any', 'all', 'xor', 'nand', 'nor', 'xnor' to specify your logic gate.")
  }
  
  #check require_identical fields arguments.
  if (!is.logical(useRequireIdenticalLogic)){
    stop("Please set useRequireIdenticalLogic to TRUE or FALSE.")
  }
  if (!is.null(requireIdenticalFields) & !useRequireIdenticalLogic){
    warning("useRequireIdenticalLogic is set to FALSE, but requireIdenticalFields are supplied. These will not be used for filtering.")
  } else if (useRequireIdenticalLogic & is.null(requireIdenticalFields)){
    stop("useRequireIdenticalLogic is TRUE, but no requireIdenticalFields have been supplied.")
  } else if (!all(requireIdenticalFields %in% contrast_columns)){
    stop("Metadata columns requested in requireIdenticalFields do not appear in the supplied design matrix's contrast_columns attribute. Please ensure that the contrast_columns vector supplied to DesignModelMatrix() is correct and that all elements in the requireIdenticalFields vector are metadata column names whose values you intended to keep constant across all pairwise contrasts.")
  }
  
  #Sanitize the logicList's criterion field according to the make.names() usage in the DesignModelMatrix (where it will compare the str_split(column) values)
  #The metadata's name field is fetched from a raw value stored in attr(design, "contrast_columns"), so it doesn't need to be fixed. 
  for (logicListIndex in seq_along(logicList)){
    logicList[[logicListIndex]][[3]] <- gsub("_", ".", make.names(logicList[[logicListIndex]]))[[3]]
  }
  
  #create a nx2 array of all possible unique pairwise combinations of contrasts from the design matrix.
  contrasts <- t(utils::combn(colnames(design), m = 2))
  
  #initialize contrast indices vector (to be used to subset the total contrast list after filtering).
  filtered_contrast_indices <- c()
  #iterate through all of the contrasts, filter those that don't satisfy the logic gates/requireIdenticalFields logic.
  for (row_index in seq_len(nrow(contrasts))){
    #define the relevant metadata fields for each side of the contrast.
    positive_contrast <- unlist(strsplit(contrasts[row_index,1], split = "_"))
    names(positive_contrast) <- contrast_columns
    negative_contrast <- unlist(strsplit(contrasts[row_index,2], split = "_"))
    names(negative_contrast) <- contrast_columns
    
    #define iterator to track how many logic gates were passed for the contrast pair.
    gates_satisfied <- 0
    for (logic_gate_index in seq_along(logicList)){
      #define variables to set up the logic gate from logicList.
      field_to_check <- logicList[[logic_gate_index]][[1]]
      gate <- logicList[[logic_gate_index]][[2]]
      criterion <- logicList[[logic_gate_index]][[3]]
      #set up logic gate and check if the contrast satisfies the gate.
      if (get(gate)(
        positive_contrast[field_to_check] == criterion,
        negative_contrast[field_to_check] == criterion)){
        #if true, mark the gate as satisfied.
        gates_satisfied <- gates_satisfied + 1
        #it is possible to fully define your study design within logic_gate without using the requireIdenticalFields logic as a shortcut, so we check if it's used or not.
        if (useRequireIdenticalLogic == T){
          #check if the contrast passed all of the logic gates, if so, test for require_identical fields.
          if (gates_satisfied == length(logicList)){
            #reset counter for matching require_identical fields.
            requireIdenticalFields_satisfied <- 0
            for (require_identical in requireIdenticalFields){
              #keep track of how many require_identical fields are equal to each other.
              if (positive_contrast[require_identical] == negative_contrast[require_identical]){
                requireIdenticalFields_satisfied <- requireIdenticalFields_satisfied + 1
              }
              #if all of the require_identical fields are require_identical keep the contrast.
              if (requireIdenticalFields_satisfied == length(requireIdenticalFields)){
                filtered_contrast_indices <- c(filtered_contrast_indices, row_index)
              }
            }
          }
        } else{
          #if all of the gates are satisfied and the require_identical logic isn't used, keep the contrast.
          if (gates_satisfied == length(logicList)){
            filtered_contrast_indices <- c(filtered_contrast_indices, row_index)
          }
        }
      }
    }
  }
  #filter the contrasts using the indices that satisfy the logic gates/require_identical fields.
  #if no filtering was actually performed, filtered_contrast_indices will be NULL/empty, so we don't need to filter the contrasts.
  if (!is.null(filtered_contrast_indices)) {
    contrasts <- contrasts[filtered_contrast_indices,]
  } else {
    warning("No filtering in FilterPseudobulkContrasts() was performed.")
  }
  
  #if the contrasts matrix only has one row, it will be automatically coerced into a vector. We need to ensure it remains a matrix.
  if (is.vector(contrasts)) {
    contrasts <- matrix(contrasts, ncol = 2)
  }
  
  if (is.null(contrasts)) {
    stop("All contrasts were filtered. Please ensure the first entry of logicList's lists is a field supplied to the contrast_columns argument of DesignModelMatrix(). Please also check that the third entry of the logicList's lists reacts predictably with the make.names() function.")
  } else if (nrow(contrasts) == 0) {
    stop("All contrasts were filtered. Please ensure the first entry of logicList's lists is a field supplied to the contrast_columns argument of DesignModelMatrix(). Please also check that the third entry of the logicList's lists reacts predictably with the make.names() function.")
  }
  
  #initialize a dataframe to store the filtered contrasts.
  filteredContrastsDataframe <- data.frame()
  #iterate through the filtered contrasts and coerce the matrix into a data frame. 
  for (contrast_index in seq_len(nrow(contrasts))){
    #split the contrasts by _ and coerce into a dataframe. 
    temporary_data_frame <- t(as.data.frame(unlist(strsplit(contrasts[contrast_index,], split = "_"))))
    #assign the column names of the dataframe via the ordering supplied by contrast_columns.
    colnames(temporary_data_frame) <- c(paste0("positive_contrast_", contrast_columns), 
                                        paste0("negative_contrast_", contrast_columns))
    rownames(temporary_data_frame) <- contrast_index
    filteredContrastsDataframe <- rbind(filteredContrastsDataframe, temporary_data_frame)
  }
  #write the filtered output file in "data.frame" format, with column names and no row indices. 
  write.table(filteredContrastsDataframe, file = filteredContrastsOutputFile, row.names = F, col.names = T)
  return(filteredContrastsDataframe)
}

#logic gates to be used with FilterPseudobulkContrasts()
nand <- function(x,y){ !(x&y) } 
nor <- function(x,y){ !(x|y) } 
xnor <- function(x,y){(x == y)}

#' @title RunFilteredContrasts
#'
#' @description This is the companion function to FilterPseudobulkContrasts. FilterPseudobulkContrasts determines the valuable contrasts relative to your metadata variable of interest, and this function performs the model fitting and collecting the differential expression results.
#' @param seuratObj A pseudobulked Seurat object that was used for DesignModelMatrix() and FilterPseudobulkContrasts(). 
#' @param filteredContrastsFile A file containing the output of FilterPseudobulkContrasts. This can be omitted if you use the dataframe returned by FilterPseudobulkContrasts insetad.
#' @param filteredContrastsDataframe A dataframe containing the output of FilterPseudobulkContrasts. This can be omitted if you read the results of FilterPseudobulkContrasts from the file written by the function instead. 
#' @param design a design/model matrix returned by DesignModelMatrix().
#' @param test.use A passthrough argument specifying the statistical test to be run by PerformGlmFit and PerformDifferentialExpression. 
#' @param logFC_threshold A passthrough argument specifying the log fold change threshold to be used by PerformDifferentialExpression (for plotting only). 
#' @param FDR_threshold A passthrough argument specifying the FDR threshold to be used by PerformDifferentialExpression (for plotting only).
#' @param filterGenes A passthrough argument specifying whether or not PerformGlmFit should filter genes using edgeR::filterByExpr.
#' @param assayName A passthrough argument specifying in which assay the counts are held. 
#' @param showPlots A passthrough argument specifying whether or not PerformDifferentialExpression should show the volcano plots.
#' @return A list of dataframes containing the differentially expressed genes in each contrast supplied by one of the filteredGenes arguments. 
#' @export

RunFilteredContrasts <- function(seuratObj, filteredContrastsFile = NULL, filteredContrastsDataframe = NULL, design, test.use, logFC_threshold = 1, filterGenes = TRUE, FDR_threshold = 0.05, assayName = "RNA", showPlots = FALSE){
  if (is.null(filteredContrastsFile) & is.null(filteredContrastsDataframe)){
    stop("Please define either filteredContrastsDataframe or filteredContrastsFile. Both of these are output by FilterPseudobulkContrasts() and are equivalent. filteredContrastsDataframe is returned by FilterPseudobulkContrasts(), and filteredContrastsFile is written to disk (by default: ./filtered_contrasts.tsv)")
  } else if (!is.null(filteredContrastsFile) & !is.null(filteredContrastsDataframe)){
    warning("Both filteredContrastsFile and filteredContrastsDataframe are defined. Only one of these is necessary, so filteredContrastsDataframe will be used as the list of filtered contrasts and the filteredContrastsFile argument will be ignored.")
    contrasts <- filteredContrastsDataframe
  } else if (!is.null(filteredContrastsDataframe) ){
    contrasts <- filteredContrastsDataframe
  } else if (!is.null(filteredContrastsFile)){
    contrasts <- read.table()
  }
  
  if (!("counts" %in% SeuratObject::Layers(seuratObj))) {
    stop("Could not find the 'counts' layer in the supplied Seurat Object.")
  }
  
  #Check for extraneous assays and remove them
  if (!all(Assays(seuratObj) %in% assayName)) {
    extraAssays <- Assays(seuratObj)[!Assays(seuratObj) %in% assayName]
    for (extraAssay in extraAssays) {
      seuratObj[[extraAssay]] <- NULL
    }
  }
  
  results <- future.apply::future_lapply(split(contrasts, seq_len(nrow(contrasts))), future.seed = GetSeed(), FUN = function(x){
    #initialize two Seurat objects, one to be subset according to the positive contrasts, and one to be subset according to the negative contrast. These will be merged downstream.
    seuratObj.positive.contrast <- seuratObj
    seuratObj.negative.contrast <- seuratObj
    positive_contrast <- NULL
    negative_contrast <- NULL
    #construct the positive and negative contrast by accessing each value within the filtered contrast data frame (row-wise) by adding the appropriate prefix to the metadata field name (i.e. the "contrast column")
    print(paste0("Contrast columns: ", attr(design,"contrast_columns")))
    for (contrast_column in attr(design, "contrast_columns")){
      #check if the contrast column in the parent Seurat object needs sanitizing before populating seuratObj.positive.contrast and seuratObj.negative.contrast downstream.
      if (!all(identical(seuratObj@meta.data[,contrast_column], .RemoveSpecialCharacters(seuratObj@meta.data[,contrast_column])))) {
        print("Converting metadata columns to a make.names() format. Hyphens, spaces, underscores, and other non-alphanumeric characters will be converted to periods. Factor levels will be retained.")
        seuratObj@meta.data[,contrast_column] <- .RemoveSpecialCharacters(seuratObj@meta.data[,contrast_column])
      }
      #if the contrasts have only just been initialized, don't use an underscore delimiter when concatenating.
      if (is.null(positive_contrast)){
        positive_contrast <- x[,paste0("positive_contrast_", contrast_column)]
        negative_contrast <- x[,paste0("negative_contrast_", contrast_column)]
      } else {
        positive_contrast <- paste0(positive_contrast, "_", x[,paste0("positive_contrast_", contrast_column)])
        negative_contrast <- paste0(negative_contrast, "_", x[,paste0("negative_contrast_", contrast_column)])
      }
      #establish a whitelist of values from the contrast_column to subset the seurat object
      #valid_values_from_contrast <- c(x[,paste0("positive_contrast_", contrast_column)], x[,paste0("negative_contrast_", contrast_column)])
      contrast_name <- paste0(positive_contrast, "-", negative_contrast)
      print(paste("Subsetting Cells for: ", contrast_name))
      #the behavior of make.names() in DesignModelMatrix calls make.names() on the whole contrast rather than it's constituent parts, which yields divergent behavior (4 -> X4) when the value starts with a number. However, this is necessary to correct for spaces/hypens in metadata variables. 
      #TODO: ensure the first contrast column in DesignModelMatrix is non-numeric.
      print("Filtering cells...")
      if(!any(grepl("^[0-9]",seuratObj.positive.contrast@meta.data[,contrast_column])) | !any(grepl("^[0-9]",seuratObj.negative.contrast@meta.data[,contrast_column]))){
        #check for factor ordering on the positive contrast Seurat object and make.names() if necessary
        seuratObj.positive.contrast@meta.data[,contrast_column] <- .RemoveSpecialCharacters(seuratObj.positive.contrast@meta.data[,contrast_column])
        #check for factor ordering on the negative contrast Seurat object and make.names() if necessary
        seuratObj.negative.contrast@meta.data[,contrast_column] <- .RemoveSpecialCharacters(seuratObj.negative.contrast@meta.data[,contrast_column])
      }
      seuratObj.contrast <-tryCatch(
        {
          #Iteratively subset SeuratObj.contrast so we fit the glm on just the samples going into the contrast.
          seuratObj.positive.contrast <- seuratObj.positive.contrast |> subset(subset = !!rlang::sym(contrast_column) %in% gsub("_",".", x[,paste0("positive_contrast_", contrast_column)]))
          seuratObj.negative.contrast <- seuratObj.negative.contrast |> subset(subset = !!rlang::sym(contrast_column) %in% gsub("_",".", x[,paste0("negative_contrast_", contrast_column)]))
          #this merge could be more efficient (specifically, this merge could be performed just once instead of once per contrast_column), but it offers a convenient check-in during the subset and offers the most correct data to eBayes().
          seuratObj.contrast <- merge(seuratObj.positive.contrast, seuratObj.negative.contrast)
          if (HasSplitLayers(seuratObj.contrast)) {
            seuratObj.contrast <- MergeSplitLayers(seuratObj.contrast)
          }
          print(paste(c("Remaining Cells:", nrow(seuratObj.contrast@meta.data))))
          print(paste("colsums:", table(seuratObj.contrast@meta.data[,contrast_column])))
          seuratObj.contrast
        }, error = function(e){
          print(paste0("Error subsetting in contrast: ", contrast_name, ". Column:",  contrast_column))
          print(paste("Error:", e))
          return(NULL)
        }
      )
    }
    
    if (!is.null(seuratObj.contrast)){
      print(paste("Performing DE for: ", contrast_name))
      #After filtering, each side of the contrast needs to have more than 1 sample or the estimated dispersion will be NA and the GLM fitting will fail. 
      #This only tests the last contrast column for >1 samples, but this should be sufficient. If a contrast yielded a subset with no samples remaining, it should error in the trycatch and seuratObj should be NULL. 
      if (all(table(seuratObj.contrast@meta.data[,contrast_column]) > 1)) {
        filtered_design_matrix <- DesignModelMatrix(seuratObj.contrast,
                                                    contrast_columns = attr(design, "contrast_columns"),
                                                    sampleIdCol = attr(design, "sampleIdCol"))
        
        fit <- PerformGlmFit(seuratObj.contrast, design = filtered_design_matrix, test.use = test.use, assayName = assayName, filterGenes = filterGenes)
        #format contrast using limma
        contrast <- limma::makeContrasts(contrasts = contrast_name, levels = colnames(filtered_design_matrix))
        result <- PerformDifferentialExpression(fit, contrast, contrast_name, logFC_threshold = logFC_threshold, FDR_threshold = FDR_threshold, test.use = test.use, showPlots = showPlots)
        
        #calculate the percentage of expression for each gene in each side of the contrast
        percentages.positive.contrast <- SeuratObject::GetAssayData(seuratObj.positive.contrast, layer = 'pct.expression', assay = assayName)
        percentages.positive.contrast <- as.matrix(percentages.positive.contrast) %*% seuratObj.positive.contrast@meta.data$TotalCells / sum(seuratObj.positive.contrast@meta.data$TotalCells)
        percentages.negative.contrast <- SeuratObject::GetAssayData(seuratObj.negative.contrast, layer = 'pct.expression', assay = assayName)
        percentages.negative.contrast <- as.matrix(percentages.negative.contrast) %*% seuratObj.negative.contrast@meta.data$TotalCells / sum(seuratObj.negative.contrast@meta.data$TotalCells)
        
        #add the percentage of expression within each contrast (pct.1 is the positive contrast, pct.2 is the negative contrast) to the result
        result$differential_expression$table$pct.1 <- percentages.positive.contrast[result$differential_expression$table$gene,]
        result$differential_expression$table$pct.2 <- percentages.negative.contrast[result$differential_expression$table$gene,]
        
        #if no DEGs are returned, then spoof the table with a "null DEG".
        if (nrow(result$differential_expression$table)==0){
          print(paste0("empty DE results for contrast:", contrast_name))
          result <- list()
          result$differential_expression$table <- data.frame(logFC = 0,
                                     logCPM = 0,
                                     `F` = 1,
                                     PValue = 1,
                                     pct.1 = 0, 
                                     pct.2 = 0,
                                     FDR = 1,
                                     gene = "none")
        }
      } else {
        print("Not enough samples to estimate dispersion/fit GLM")
        result <- list()
        result$differential_expression$table <- data.frame(logFC = 0,
                                   logCPM = 0,
                                   `F` = 1,
                                   PValue = 1,
                                   pct.1 = 0, 
                                   pct.2 = 0,
                                   FDR = 1,
                                   gene = "none")
      }
    } else {
      print("Filtering eliminated all samples.")
      contrast_name <- paste0(positive_contrast, "-", negative_contrast)
      result <- list()
      result$differential_expression$table <- data.frame(logFC = 0,
                                 logCPM = 0,
                                 `F` = 1,
                                 PValue = 1,
                                 pct.1 = 0, 
                                 pct.2 = 0,
                                 FDR = 1,
                                 gene = "none")
    }
    result$differential_expression$table$contrast_name <- contrast_name
    result$differential_expression$table$positive_contrast <- positive_contrast
    result$differential_expression$table$negative_contrast <- negative_contrast
    result$differential_expression$table <- cbind(result$differential_expression$table, x)
    return(result$differential_expression$table)
  })
  return(results)
}
#' @title RemoveSpecialCharacters
#' 
#' @description This function is used to remove special characters from metadata values. This ensures that the seuratObj metadata columns are formatted to match entries in the design matrix. 
#' @param vector_of_metadata_values A vector of metadata values to be sanitized
#' @return A vector of metadata values with special characters removed.

.RemoveSpecialCharacters <- function(vector_of_metadata_values) { 
  if (is.factor(vector_of_metadata_values)) {
    vector_of_metadata_values  <- forcats::fct_relabel(vector_of_metadata_values, ~gsub("_", ".", make.names(.)))
  } else {
    vector_of_metadata_values <- gsub("_", ".", make.names(vector_of_metadata_values))
  }
  return(vector_of_metadata_values)
}

#' @title PseudobulkingBarPlot
#'
#' @description This is a plotting function downstream of RunFilteredContrasts that will create a  bar plot with the results of each contrast ordered by the magnitude of differentially expressed genes. 
#' @param filteredContrastsResults A list of dataframes returned by RunFilteredContrasts.
#' @param metadataFilterList An optional list of lists specifying further filtering to be performed. These lists must follow the format: list( list("filterDirection", "filterField", "filterValue" )) where: filterDirection determines whether the Positive, Negative, or Both sides of the contrasts should be filtered, filterField determines which metadata field (originally supplied to groupFields in PseudobulkSeurat and contrast_columns in DesignModelMatrix), filterValue corresponds to which values of the filterField should be filtered. This is useful if you passed parallel hypotheses (e.g. what genes are differentially expressed in each tissue?) in the requireIdenticalFields argument of RunFilteredContrasts, but want to plot the results of only one tissue at a time. 
#' @param title Title for the bar plot. 
#' @param logFC_threshold A passthrough argument specifying the log fold change threshold to be used by .addRegulationInformationAndFilterDEGs to filter genes and determine regulation direction. 
#' @param FDR_threshold A passthrough argument specifying the FDR threshold to be used by .addRegulationInformationAndFilterDEGs to filter genes and determine regulation direction.
#' @param swapContrastDirectionality A boolean determining if the contrast directionality should be swapped. This is useful if you want the "control" condition in your contrasts to appear in the opposite directionality of the default. 
#' @return A list containing the filtered dataframe used for plotting and the bar plot itself. 
#' @export

PseudobulkingBarPlot <- function(filteredContrastsResults, metadataFilterList = NULL, title = "Please Title The Bar Plot", logFC_threshold = 1, FDR_threshold = 0.05, swapContrastDirectionality = FALSE) {
  
  if (!is.list(filteredContrastsResults)) {
    stop("filteredContrastsResults is not a list. Please ensure filteredContrastsResults is a list of dataframes returned by RunFilteredContrasts().")
  } else if (length(filteredContrastsResults) == 0) {
    stop("filteredContrastsResults is empty. Please ensure filteredContrastsResults is a list of dataframes returned by RunFilteredContrasts().")
  } else if (!all(sapply(filteredContrastsResults, function(x) class(x) == "data.frame"))) {
    stop("filteredContrastsResults is not a list of dataframes. Please ensure filteredContrastsResults is a list of dataframes returned by RunFilteredContrasts().")
  } 
  
  if (!is.null(metadataFilterList)) {
    if(!is.list(metadataFilterList)){
      stop("metadataFilterList is not a list. Please ensure metadataFilterList includes a contrast directionality ('Positive', 'Negative', or 'Both'), a metadata variable that was included in groupFields, and a whitelist value that is a member of the metadata variable. Proper formatting should look like: \nmetadataFilterList = list(list('Positive', 'Tissue', 'Lung'), list('Negative', 'Timepoint', 'Baseline')).")
    } else if (length(metadataFilterList) == 0) {
      stop("metadataFilterList is empty. Please ensure metadataFilterList includes a contrast directionality ('Positive', 'Negative', or 'Both'), a metadata variable that was included in groupFields, and a whitelist value that is a member of the metadata variable. Proper formatting should look like: \nmetadataFilterList = list(list('Positive', 'Tissue', 'Lung'), list('Negative', 'Timepoint', 'Baseline')).")
    } else if (length(unlist(metadataFilterList)) %% 3 != 0) {
      stop("Not all lists in metadataFilterList are of length 3! Please ensure metadataFilterList includes a contrast directionality ('Positive', 'Negative', or 'Both'), a metadata variable that was included in groupFields, and a whitelist value that is a member of the metadata variable. Proper formatting should look like: \nmetadataFilterList = list(list('Positive', 'Tissue', 'Lung'), list('Negative', 'Timepoint', 'Baseline')).")
    }
  }
  
  #convert from list of dataframes to one large data frame (note: row names will be duplicated when genes appear more than once).
  filteredContrastsResults <- dplyr::bind_rows(filteredContrastsResults)
  #tag the genes as either up or down regulated
  filteredContrastsResults <- .addRegulationInformationAndFilterDEGs(filteredContrastsResults, logFC_threshold = logFC_threshold, FDR_threshold = FDR_threshold)
  
  #check for DEGs. If there are none, raise an error. 
  if (all(unique(filteredContrastsResults$n_DEG %in% c(0)))) {
    stop("All of the genes in all of the contrasts failed to pass the FDR and logFC thresholds. You can consider adjusting the logFC_threshold and FDR_threshold arguments, but this is a reasonable result (i.e. no DEGs) when you are comparing very similar samples.")
  }
  
  #Further filter contrasts associated with a list of vectors (metadataFilterList).
  if (!is.null(metadataFilterList)) {
    if (!is.list(metadataFilterList)) {
      stop("metadataFilterList is improperly formatted. Please ensure metadataFilterList includes a contrast directionality ('Positive', 'Negative', or 'Both'), a metadata variable that was included in groupFields, and a whitelist value that is a member of the metadata variable. Proper formatting should look like: metadataFilterList = list(c('Positive', 'Tissue', 'Lung')).")
    } else if (length(unlist(metadataFilterList)) %% 3 != 0) {
      stop("Not all vectors in metadataFilterList are of length 3! Please ensure metadataFilterList includes a contrast directionality ('Positive', 'Negative', or 'Both'), a metadata variable that was included in groupFields, and a whitelist value that is a member of the metadata variable. Proper formatting should look like: metadataFilterList = list(c('Positive', 'Tissue', 'Lung')).")
    } 
    #parse the filterList for the metadata field and whitelist value 
    for (filter in metadataFilterList) {
      #define filters
      filterDirection <- filter[1]
      filterField <- filter[2]
      filterValue <- filter[3]
      
      #Make sure the filters are valid 
      if (! (filterDirection %in% c("Positive", "Negative", "Both"))) {
        stop(paste0("Error: Invalid argument ", filterDirection, " supplied in metadataFilterList. Please ensure the first argument in each vector of metadataFilterList specifies a direction that is either 'Positive' for filtering only the positive side of each contrast, 'Negative' for filtering the negative side, or 'Both' for filtering values from either side of each contrast."))
      } else if (!any(grepl(filterField), colnames(filteredContrastsResults))) {
        stop(paste0("Error: No contrast columns associated with ", filterField, " were found in the column names of filteredContrastsResults. Please ensure your metadata column was both supplied to groupFields in PseudobulkSeurat and designated as a contrastColumn in DesignModelMatrix."))
      } else if (!any(grepl(filterValue, filteredContrastsResults[,paste0("positive_contrast_", filterField)]) | 
                    grepl(filterValue, filteredContrastsResults[,paste0("negative_contrast_", filterField)]))){
        warning(paste0(filterValue, " was not found in the columns associated with ", filterField,". No filtering was performed associated with the following filter: ", filter))
      }
      
      #Apply post-RunFilteredContrasts filtering rules to the contrasts
      if (filterDirection == "Positive" | filterDirection == "Both") {
        filteredContrastsResults <- filteredContrastsResults %>% 
          dplyr::filter(grepl(filterValue, paste0("positive_contrast_", filterField)))
      } 
      if (filterDirection == "Negative" | filterDirection == "Both") {
        filteredContrastsResults <- filteredContrastsResults %>% 
          dplyr::filter(grepl(filterValue, paste0("negative_contrast_", filterField)))
      }
      
    }
  }
  
  if (swapContrastDirectionality){
    filteredContrastsResults <- .swapContrast(filteredContrastsResults)
  }
  
  #Define order of magnitude DEG groupings
  filteredContrastsResults <- filteredContrastsResults %>% 
    #filter genes from failed model fits
    dplyr::filter(gene != "none") %>% 
    dplyr::filter(!is.na(uniqueness)) %>% 
    #compute the number of DEGs for each group of DEGs (unique vs up/down regulated)
    dplyr::group_by(contrast_name, across(starts_with("positive_contrast")), across(starts_with("negative_contrast"))) %>% 
    reframe(number_of_positive_nonunique_DEGs = sum(n_DEG[uniqueness %in% c("up_nonunique")]), 
            number_of_positive_unique_DEGs = sum(n_DEG[uniqueness %in% c("up_unique")]),
            number_of_negative_nonunique_DEGs = sum(n_DEG[uniqueness %in% c("down_nonunique")]), 
            number_of_negative_unique_DEGs = sum(n_DEG[uniqueness %in% c("down_unique")]), 
            ) %>%     
    group_by(contrast_name) %>% 
    mutate(total_DEGs = max(number_of_positive_nonunique_DEGs) + 
             max(number_of_positive_unique_DEGs) + 
             max(abs(number_of_negative_nonunique_DEGs)) + 
             max(abs(number_of_negative_unique_DEGs))) %>%
    unique.data.frame() %>% 
    dplyr::mutate(DEG_Magnitude = dplyr::case_when(
      total_DEGs > 1000 ~ "1000+ DEGs",
      total_DEGs <= 1000 & total_DEGs >= 100 ~ "1000-100 DEGs",
      total_DEGs <= 100 & total_DEGs >= 10 ~ "100-10 DEGs",
      total_DEGs < 10 ~ "<10 DEGs"
    ))
  
  filteredContrastsResults$DEG_Magnitude <- factor(filteredContrastsResults$DEG_Magnitude, levels = c("1000+ DEGs", "1000-100 DEGs", "100-10 DEGs", "<10 DEGs"))
  

  
  #plot the bar graph
  #note: geom_bar doesn't play nicely with y axis transformations, so scales = free_y and log transforming the y axis does not work. 
    bargraph <- ggplot2::ggplot(filteredContrastsResults) + 
      ggplot2::geom_col(ggplot2::aes(x = stats::reorder(contrast_name, -abs(total_DEGs)), y = number_of_positive_nonunique_DEGs, fill = "up_nonunique"), position="stack") + 
      ggplot2::geom_col(ggplot2::aes(x = stats::reorder(contrast_name, -abs(total_DEGs)), y = number_of_negative_nonunique_DEGs, fill = "down_nonunique"), position="stack") + 
      #plot downregulated genes
      ggplot2::geom_col(ggplot2::aes(x = stats::reorder(contrast_name, -abs(total_DEGs)), y = number_of_positive_unique_DEGs, fill = "up_unique"), position="stack") + 
      ggplot2::geom_col(ggplot2::aes(x = stats::reorder(contrast_name, -abs(total_DEGs)), y = number_of_negative_unique_DEGs, fill = "down_unique"), position="stack", ) + 
      #define fixed color scheme
      ggplot2::scale_fill_manual(values = c("up_unique" = "red",
                                            "up_nonunique" = "orange", 
                                            "down_unique" = "blue", 
                                            "down_nonunique" = "cadetblue2"),
                                 breaks = c("up_unique", "up_nonunique", "down_unique", "down_nonunique"),
                                 aesthetics = "fill") + 
      ggplot2::labs(fill="Gene Uniqueness") + 
      ggplot2::ylab("Number of DEGs")+ 
      egg::theme_article() + 
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) + 
      ggplot2::ggtitle(title) + 
      ggplot2::xlab("Differential Expression Contrasts") + 
      ggplot2::facet_grid(~DEG_Magnitude, scales = "free_x")
  
  print(bargraph)
  
  return(list(dataframe = filteredContrastsResults, barPlot = bargraph))
}

#' @title .addRegulationInformationAndFilterDEGs
#' 
#' @description a helper function for PseudobulkingBarPlot to tag genes as up- or down-regulated and check for uniqueness of differential gene expression across the contrasts. 
#' @param tibble a (presumably large) tibble from dplyr containing the results of RunFilteredPseudobulkingContrasts
#' @param logFC_threshold The minimum value for a gene's absolute value of log fold change to be considered differentially expressed. 
#' @param FDR_threshold The maximum value for a gene's FDR to be considered differentially expressed. 
#' @return a tibble with the uniqueness and direction of regulation for each gene stored in tibble$regulation

.addRegulationInformationAndFilterDEGs <- function(tibble, logFC_threshold = 1, FDR_threshold = 0.05){
  tibble <- dplyr::as_tibble(tibble)
  tibble <- tibble %>% 
    dplyr::group_by(contrast_name, gene) %>% 
    dplyr::mutate(regulation = dplyr::case_when( logFC > logFC_threshold & FDR < FDR_threshold ~ "upregulated",
                                   logFC < -logFC_threshold & FDR < FDR_threshold ~ "downregulated", 
                                   abs(logFC) < logFC_threshold | FDR > FDR_threshold ~ "noDEG"))
  
  tibble <- tibble |>
    dplyr::group_by(gene) |>
    dplyr::mutate(gene_occurance = dplyr::n()) 
  
  #determine regulation-uniqueness pairs
  tibble[tibble$regulation == "upregulated","uniqueness"] <- "up_nonunique"
  tibble[tibble$regulation == "upregulated" & tibble$gene_occurance == 1, "uniqueness"] <- "up_unique"
  tibble[tibble$regulation == "downregulated","uniqueness"] <- "down_nonunique"
  tibble[tibble$regulation == "downregulated" & tibble$gene_occurance == 1, "uniqueness"] <- "down_unique"
  
  #geom_bar sums the n_DEG column, so fixing the sign is all that is necessary
  tibble[tibble$regulation == "upregulated", "n_DEG"] <- 1 
  tibble[tibble$regulation == "downregulated", "n_DEG"] <- -1 
  tibble[tibble$regulation == "noDEG", "n_DEG"] <- 0
  return(tibble)
}

#' @title PseudobulkingDEHeatmap
#'
#' @description This is a plotting function that formats the FilterPseudobulkingContrasts logic to plot a heatmap of log fold changes for a given gene space. It is recommended that you filter this gene space through some kind of feature selection, but full transcriptome heatmaps are possible.  
#' @param seuratObj A pseudobulked Seurat object.
#' @param geneSpace An optional vector of gene names that specify which genes should populate the heatmap. It is recommended that you perform some kind of feature selection upstream of this to select genes for plotting. If NULL, all rows in the assay object will be used.
#' @param contrastField The primary grouping variable to display (on top of the heatmap). 
#' @param negativeContrastValue The value of contrastField to be treated as "downregulated". 
#' @param positiveContrastValue An optional variable to define a specific positive contrast value. While negativeContrastValue determines the log fold changes, this argument operates primarily as a filtering variable to eliminate groups and show a particular value of the contrast field.
#' @param subgroupingVariable If there is a second relevant grouping variable, you can supply a second metadata column to group the data by. 
#' @param assayName the name of the assay in the seurat object storing the count matrix. 
#' @param showRowNames a passthrough variable for ComplexHeatmap controlling if the gene names should be shown or not in the heatmap. 
#' @param sampleIdCol The metadata column denoting the variable containing the sample identifier (for grouping). 
#' @param filterGenes A passthrough variable for PerformGlmFit, used to determine if genes should be filtered using edgeR::filterByExpr.
#' @param subsetExpression An optional string containing an expression to subset the Seurat object. This is useful for selecting an exact subpopulation to in which to show DEGs. Please note that for string-based metadata fields, you will need to mix single and double quotes to ensure your expression is properly parsed. For instance, note the double quotes around the word unvax in this expression: subsetExpression = 'vaccine_cohort == "unvax"'. 
#' @return A list containing the filtered dataframe used for plotting and the heatmap plot itself. 
#' @export

PseudobulkingDEHeatmap <- function(seuratObj, geneSpace = NULL, contrastField = NULL, negativeContrastValue = NULL, positiveContrastValue = NULL, subgroupingVariable = NULL, showRowNames = FALSE, assayName = "RNA", sampleIdCol = 'cDNA_ID', filterGenes = FALSE, subsetExpression = NULL) {
  #sanity check arguments
  if (is.null(contrastField)) {
    stop("Please define a contrastField. This is a metadata variable (supplied to groupFields during PseudobulkSeurat()) that will be displayed on the top of the heatmap.")
  }
  #check negativeContrastValue & positiveCOntrastValue
  if (is.null(negativeContrastValue)) {
    stop("Please define a negativeContrastValue. This is the value of contrastField that will be treated as 'downregulated' in log fold changes shown in the heatmap.")
  } else if (!(negativeContrastValue %in% seuratObj@meta.data[,contrastField])) {
    stop(paste0("Error: could not find negativeContrastValue: ", negativeContrastValue, " in the metadata field ", contrastField, " within the Seurat Object. Please ensure the negativeContrastValue is a member of the ",  contrastField, " metadata field."))
  } else if (!is.null(positiveContrastValue)) {
    if (!(positiveContrastValue %in% seuratObj@meta.data[,contrastField])) {
      stop(paste0("Error: could not find positiveContrastValue: ", positiveContrastValue, " in the metadata field ", contrastField, " within the Seurat Object. Please ensure the positiveContrastValue is a member of the ",  contrastField, " metadata field."))
    }
  }
  #check subgroupingVariable
  if (!is.null(subgroupingVariable)) {
    if (!(subgroupingVariable %in% colnames(seuratObj@meta.data))) {
      stop(paste0("Error: could not find subgroupingVariable: ", subgroupingVariable, " in the metadata fields of the Seurat Object. Please ensure the subgroupingVariable: ", subgroupingVariable, " is a member of the metadata and was passed as a groupField to PseudobulkSeurat()."))
    }
  }
  #check geneSpace
  if (length(geneSpace) == 1) {
    stop("Error: geneSpace must contain more than one gene. Please ensure geneSpace is a vector of gene names that you would like to plot in the heatmap.")
  } else if (is.null(geneSpace)) {
    stop("Error: please define a geneSpace. This should be a vector of gene names within the pseudobulked Seurat object.")
  } else if (all(!(geneSpace %in% rownames(seuratObj)))) {
    stop("Error: no genes in geneSpace were found in the supplied seurat object. Please ensure geneSpace is a vector of character gene names that you would like to plot in the heatmap.")
  } else if (any(!(geneSpace %in% rownames(seuratObj)))) {
    warning(paste0("Genes: ", geneSpace[!(geneSpace %in% rownames(seuratObj))], " not found in supplied Seurat Object! These genes will be omitted."))
  } else if (sum(geneSpace %in% rownames(seuratObj)) <= 1) {
    stop(paste0('Only ', sum(geneSpace %in% rownames(seuratObj)), ' gene(s) overlapped between geneSpace (length ', length(geneSpace),') and the assay (size ', length(rownames(seuratObj)),'). With Seurat V5 there must be more than one feature'))
  }
  
  if (!is.null(subsetExpression)) { 
    if (!rlang::is_string(subsetExpression)) {
      stop("Error: subsetExpression must be a string containing the expression you wish to pass to the subset function. An example of a valid subsetExpression is: subsetExpression = 'ClusterNames_0.2 == 0'. For string matching, you may mix single and double quotes to ensure your expression is properly parsed. For more details, please see the documentation ?PseudobulkingDEHeatmap.")
    } else {
      seuratObj <- subset(seuratObj, cells = eval(parse(text = paste0("Seurat::WhichCells(seuratObj, expression = ", subsetExpression, ")"))))
    }
      }
  
  #Sanitize negativeContrastValue and positiveContrastValue, since the user is unlikely to know that values need to be compatible with a post-make.names() call to the variables from the design matrix. 
  if (negativeContrastValue != .RemoveSpecialCharacters(negativeContrastValue)) { 
    negativeContrastValue <- .RemoveSpecialCharacters(negativeContrastValue)
  } 
  if (!is.null(positiveContrastValue)) {
    if (positiveContrastValue != .RemoveSpecialCharacters(positiveContrastValue)) { 
      positiveContrastValue <- .RemoveSpecialCharacters(positiveContrastValue) 
    } 
  }
  
  #parse the contrastField, contrastValues arguments, and sampleIdCol to construct the model matrix for performing the desired contrast for the heatmap.
  design <- DesignModelMatrix(seuratObj, contrast_columns = c(contrastField, subgroupingVariable), sampleIdCol = sampleIdCol)
  #similarly, parse these arguments for setting up a logicList
  logicList <- list(list(contrastField, "xor", negativeContrastValue))
  #positiveContrastValue is mostly a filtering tool since we need to establish what our 'control' is via negativeContrastValue.  
  if (!is.null(positiveContrastValue)) {
    logicList[[2]] <- list(contrastField, 'xor', positiveContrastValue)
  }
  
  if (!is.null(subgroupingVariable)) {
    #use subgroupingVariable information to infer arguments for FilterPseudobulkContrasts
    useRequireIdenticalLogic <- TRUE
    requireIdenticalFields <- subgroupingVariable
    filteredContrasts <- FilterPseudobulkContrasts(logicList = logicList, 
                                                   design = design, 
                                                   useRequireIdenticalLogic = useRequireIdenticalLogic, 
                                                   requireIdenticalFields = requireIdenticalFields, 
                                                   filteredContrastsOutputFile = tempfile())
  } else {
    #use subgroupingVariable information to infer arguments for FilterPseudobulkContrasts
    useRequireIdenticalLogic <- FALSE
    filteredContrasts <- FilterPseudobulkContrasts(logicList = logicList, 
                                                   design = design, 
                                                   useRequireIdenticalLogic = useRequireIdenticalLogic, 
                                                   filteredContrastsOutputFile = tempfile())
  }
  
  #we can run RunFilteredContrasts without any differential expression-based filtering because our gene filtering should occur when we pass-in our geneSpace variable. We're just populating log fold changes here. 
  lfc_results <- RunFilteredContrasts(seuratObj = seuratObj, 
                                               filteredContrastsDataframe = filteredContrasts, 
                                               design = design,
                                               test.use = "QLF", 
                                               logFC_threshold = 0,
                                               FDR_threshold = 1.1,
                                               filterGenes = filterGenes, 
                                               assayName = assayName)
  #bind RunFilteredContrasts into a dataframe
  lfc_results <- dplyr::bind_rows(lfc_results)
  
  #swap directionality of contrasts such that negativeContrastValue is always negative if necessary. Contrasts are populated alphabetically upstream, so it's easier to just ensure ordering here. 
  lfc_proper_negative_contrast <- lfc_results %>% 
    dplyr::filter(!!sym(paste0("negative_contrast_",contrastField)) == negativeContrastValue)
  lfc_swapped_negative_contrast <- lfc_results %>% 
    dplyr::filter(!!sym(paste0("negative_contrast_",contrastField)) != negativeContrastValue) %>% 
    .swapContrast()
  lfc_results <- rbind(lfc_proper_negative_contrast, lfc_swapped_negative_contrast)
  
  #filter the matrix to include just genes passed to geneSpace
  lfc_results <- lfc_results %>% filter(gene %in% geneSpace)
  
  #pivot the long-form tibble into a matrix for complexheatmap
  #if we have a second grouping variable (subgroupingVariable), we include it in the else statement. 
  if (is.null(subgroupingVariable)) {
    lfc_results_wide <- lfc_results %>% 
      as.data.frame() %>% 
      dplyr::select(all_of(c("gene", "logFC", paste0("positive_contrast_",contrastField)))) %>% 
      tidyr::pivot_wider(values_from = logFC, names_from = !!sym(paste0("positive_contrast_",contrastField)))
  } else {
    lfc_results_wide <- lfc_results %>% 
      tidyr::unite(joinedFields, c(paste0("positive_contrast_",contrastField), paste0("positive_contrast_",subgroupingVariable)), remove = FALSE) %>% 
      as.data.frame() %>% 
      dplyr::select(all_of(c("gene", "logFC", "joinedFields"))) %>% 
      tidyr::pivot_wider(values_from = logFC, names_from = joinedFields )
  }
  #format the tibble from dplyr into a numeric matrix for ComplexHeatmap
  heatmap_matrix <- as.data.frame(lfc_results_wide)
  gene_vector <- heatmap_matrix[,"gene"] 
  contrast_values_vector <- colnames(heatmap_matrix)[colnames(heatmap_matrix) != "gene"] #if the heatmap matrix ends up being a single column after we strip away the gene column, it will be converted to a vector and lose its name.
  heatmap_matrix <- heatmap_matrix[,!colnames(heatmap_matrix) %in% "gene"]
  heatmap_matrix <- as.matrix(sapply(heatmap_matrix, as.numeric)) 
  rownames(heatmap_matrix) <- gene_vector 
  colnames(heatmap_matrix) <- contrast_values_vector #reassign the contrast values to the column names in case heatmap_matrix became a vector in the above lines. 
  
  #If a contrast "failed" due to not enough samples to fit a GLM in the contrast, it will return an NA LFC and Pvalue. Impute this to zero. 
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  #If there's no second grouping variable, compute a heatmap with just the grouping labels. 
  #else, compute a heatmap with the top labels as positive values of contrastField, while the bottom labels are entries of subgroupingVariable.
  if (is.null(subgroupingVariable)) {
    
    split <- colnames(heatmap_matrix)
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill =  rep(x = "white", length(colnames(heatmap_matrix)))), 
                       labels = colnames(heatmap_matrix)))
    #force a symmetric and zero-centered color scale
    heatmap_extreme_value <- max(abs(min(heatmap_matrix)), abs(max(heatmap_matrix)))
    col_fun <- circlize::colorRamp2(c(-abs(heatmap_extreme_value), 0, abs(heatmap_extreme_value)), c("dodgerblue3", "white", "red"))
    
    heatmap <- ComplexHeatmap::Heatmap(heatmap_matrix,
                                       name = "Log Fold Changes", 
                                       col = col_fun,
                                       show_row_names = showRowNames, 
                                       cluster_columns = FALSE, 
                                       top_annotation = top_annotation, 
                                       show_column_names = FALSE, 
                                       column_split = split, 
                                       column_title = NULL, 
                                       border = TRUE)
    
  } else {
    labels_df <- do.call(args = strsplit(colnames(heatmap_matrix), split = "_"), rbind)
    split <- labels_df[,1]
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = rep(x = "white", length(unique(labels_df[,1])))), 
                       labels = unique(labels_df[,1])))
    #force a symmetric and zero-centered color scale
    heatmap_extreme_value <- max(abs(min(heatmap_matrix)), abs(max(heatmap_matrix)))
    col_fun <- circlize::colorRamp2(c(-abs(heatmap_extreme_value), 0, abs(heatmap_extreme_value)), c("dodgerblue3", "white", "red"))
    
    heatmap <- ComplexHeatmap::Heatmap(heatmap_matrix,
                            column_labels = labels_df[,2],
                            col = col_fun,
                            name = "Log Fold Changes", 
                            column_names_rot = 0,
                            column_names_centered = T,
                            show_row_names = showRowNames, 
                            cluster_columns = FALSE, 
                            top_annotation = top_annotation, 
                            column_split = split, 
                            column_title = NULL, 
                            border = TRUE)
  }
  
  return(list(heatmap = heatmap, matrix = heatmap_matrix))
}


#' @title .swapContrast
#' 
#' @description a helper function for PseudobulkingDEHeatmap to swap the "positive" and "negative" sides of the contrasts such that the desired negative contrast value is correctly designated as negative. 
#' @param tibble a large tibble from dplyr containing the results of RunFilteredPseudobulkingContrasts called within PseudobulkingDEHeatmap.
#' @return a tibble with the positive and negative contrast directions swapped.

.swapContrast <- function(tibble){
  swapped_tibble <- tibble
  #swap positive contrast columns to negative
  swapped_tibble[,grepl("negative_contrast", colnames(tibble))] <- tibble[,grepl("positive_contrast", colnames(tibble))]
  #swap negative contrasts to positive
  swapped_tibble[,grepl("positive_contrast", colnames(tibble))] <- tibble[,grepl("negative_contrast", colnames(tibble))]
  #swap direction of log fold change
  swapped_tibble$logFC <- -swapped_tibble$logFC
  #swap the contrast names
  split_contrast_names <- strsplit(swapped_tibble$contrast_name, split = '-')
  swapped_tibble$contrast_name <- sapply(split_contrast_names, FUN = function(x) { paste0(x[2], "-", x[1])})
  #swap the percentage expression values
  swapped_tibble$pct.1 <- tibble$pct.2
  swapped_tibble$pct.2 <- tibble$pct.1
  return(swapped_tibble)
}

#' @title FitRegularizedClassificationGlm 
#'
#' @description Treating gene expression like a classification problem, this function trains a penalized model to classify a metadata feature. 
#' @param seuratObj a Seurat object
#' @param metadataVariableForClassification The metadata feature to be classified. If non-binary, then multinomial regression will automatically be performed. 
#' @param rescale The feature selection will optimize for "heatmap-interpretable genes" so the features are intended to be scaled. If TRUE, this will rescale the variable features.
#' @param numberOfVariableFeatures A parameter to select how many features should be selected as variable for scaling, by default, all genes will be used. 
#' @param assay Seurat Object's assay
#' @param layer layer within the Seurat object assay. Recommended to be "scale.data".
#' @param devianceCutoff Tolerance for model error when deciding how much regularization should be performed. 1 = no tolerance for error, 0 = intercept only, no genes used for prediction.
#' @param split the option to provide a previous model's training/testing set. This is necessary if you're performing multiple iterations of model fitting. 
#' @param returnModelAndSplits A boolean option to return a list containing the fitted model and training/testing splits in addition to the useful features. 
#' @param excludeVariableGenes A boolean option to exclude variable genes from the analysis. Only the Variable Gene Exclusion (v2) list from RIRA is supported.
#' @return A vector of genes useful for classification and, optionally, the model and training/testing sets.
#' @import mlr3verse
#' @importFrom mlr3 as_task_classif lrn partition
#' @export

FitRegularizedClassificationGlm <- function(seuratObj,
                                            metadataVariableForClassification = NULL,
                                            rescale = TRUE,
                                            numberOfVariableFeatures = 3000,
                                            assay = "RNA",
                                            layer = "scale.data",
                                            devianceCutoff = 0.8,
                                            split = NULL,
                                            returnModelAndSplits = F, 
                                            excludeVariableGenes = T) {
  #sanity check arguments
  if (is.null(metadataVariableForClassification)){
    stop("Please supply a column of the seurat object's metadata to classify.")
  }
  if (!(metadataVariableForClassification %in% colnames(seuratObj@meta.data))){
    stop("Supplied metadataVariableForClassification not found in the seurat object's metadata. Please ensure your metadata column is spelled correctly and exists in seuratObj@meta.data.")
  }

  set.seed(GetSeed())
 
  #convert the scale.data matrix to include a labeled classification column
  target_labeled_data <- .PrepareDataForClassificationGlm(seuratObj, 
                                                          assay = assay, 
                                                          layer = layer, 
                                                          rescale = rescale,
                                                          numberOfVariableFeatures = numberOfVariableFeatures, 
                                                          excludeVariableGenes = excludeVariableGenes, 
                                                          metadataVariableForClassification = metadataVariableForClassification)

  #set up task
  task_metadata_classification <- mlr3::as_task_classif(target_labeled_data,
                                                        target = metadataVariableForClassification,
                                                        id = "metadata_classification")
  #if a training set wasn't provided, create one, otherwise use the supplied split.
  if (is.null(split)){
    split <- mlr3::partition(task_metadata_classification)
  }

  #cv_glmnet to parameter scan regularization
  learner <- mlr3::lrn("classif.cv_glmnet")
  learner$train(task = task_metadata_classification, row_ids = split$train)
  #find %deviance over parameter scan
  deviance_vector <-
    round(1 - (stats::deviance(learner$model$glmnet.fit) / learner$model$glmnet.fit$nulldev), 2)
  #select minimum lambda value above deviance cutoff
  deviance_cutoff_index <- min(which(deviance_vector > devianceCutoff))

  #iterate through the classes of beta coefficients in the case of multinomial regression and collect useful features (genes) from each class
  classification_features <- c()
  if(learner$model$call$family == "binomial"){
    classification_features <- learner$selected_features()
  } else {
    for (class in names(learner$model$glmnet.fit$beta)) {
      #apply deviance cutoff to select lambda value
      class_weights_vector <-
        learner$model$glmnet.fit$beta[[class]][, deviance_cutoff_index]
      #harvest genes with non-zero beta coefficients
      classification_features_for_class <-
        names(class_weights_vector[abs(class_weights_vector) > 0])
      #store selected genes
      classification_features <- c(classification_features,
                                   classification_features_for_class)
    }

    #put the dashes back in the feature names, delete the "leadingNumber" prefix, and replace "star" with asterisks.
    classification_features <- gsub("dash", "-", classification_features)
    classification_features <- gsub("^leadingNumber", "", classification_features)
    classification_features <- gsub("star", "*", classification_features)
  }

  #return either a vector of genes or both a model and vector of genes.
  if (!returnModelAndSplits) {
    return(classification_features)
  } else {
    return(list(
      classification_features = classification_features,
      model = learner$model,
      split = split
    ))
  }
}

#' @title .PrepareDataForClassificationGlm
#' @description a helper function for FitClassificationGlm to remove special characters from the gene names.
#' @param seuratObj Input seurat object used for training a glmnet model. 
#' @param assay Seurat Object's assay
#' @param layer layer within the Seurat object assay. Recommended to be "scale.data".
#' @param rescale The feature selection will optimize for "heatmap-interpretable genes" so the features are intended to be scaled. If TRUE, this will rescale the variable features.
#' @param numberOfVariableFeatures A parameter to select how many features should be selected as variable for scaling, by default, all genes will be used.
#' @param excludeVariableGenes A boolean option to exclude variable genes from the analysis. Only the Variable Gene Exclusion list (v2) from RIRA is supported. 
#' @param metadataVariableForClassification The metadata feature to be classified. If non-binary, then multinomial regression will automatically be performed.
#' 
#' @return a matrix with the special characters removed from the gene names.

.PrepareDataForClassificationGlm <- function(seuratObj, 
                                             assay = "RNA", 
                                             layer = 'scale.data', 
                                             rescale = TRUE, 
                                             numberOfVariableFeatures = NULL,  
                                             excludeVariableGenes = T, 
                                             metadataVariableForClassification = NULL) {
  # rescale the input data (in case of an upstream subset since it was last rescaled).
  
  if (excludeVariableGenes) {
    exclusionList <- NULL
  } else {
    exclusionList <- RIRA::GetGeneSet("VariableGenes_Exclusion.2")
  }
  
  if (rescale) {
    if (is.null(numberOfVariableFeatures)) {
      seuratObj <- CellMembrane::NormalizeAndScale(seuratObj,
                                                   nVariableFeatures = length(rownames(seuratObj)),
                                                   variableGenesBlacklist = exclusionList,
                                                   scoreCellCycle = F, 
                                                   featuresToRegress = NULL
      )
    } else {
      seuratObj <- CellMembrane::NormalizeAndScale(seuratObj,
                                                   nVariableFeatures = numberOfVariableFeatures,
                                                   variableGenesBlacklist = exclusionList,
                                                   scoreCellCycle = F, 
                                                   featuresToRegress = NULL
      )
    }
  }
  
  target_labeled_data <-
    #merge the seuratObj's requested layer (converted to dense just in case a non-scale.data layer was used)
    merge(
      Matrix::t(as.matrix(Seurat::GetAssayData(
        seuratObj,
        assay = assay,
        layer = layer
      ))),
      seuratObj@meta.data |>
        dplyr::select(dplyr::all_of(metadataVariableForClassification)),
      by = 0) |>
    #drop rownames column
    dplyr::select(-dplyr::all_of("Row.names"))
  
  #fix gene names like MAMU-A or UGT2B9*2 with an easily replaceable and uniquely mapping "dash", "leadingNumber", or "star".
  
  # Perform a test to ensure we wont get conflicts:
  for (token in c("dash", "leadingNumber", "star")) {
    if (any(grepl(colnames(target_labeled_data), pattern = token))) {
      matches <- grep(colnames(target_labeled_data), pattern = token, value = TRUE)
      stop(paste0('The input feature names contained the unexpected pattern: ', token, '. Feature(s) were: ', paste0(matches, collapse = ', ')))
    }
  }
  colnames(target_labeled_data) <-
    gsub("-", "dash", colnames(target_labeled_data))
  colnames(target_labeled_data) <-
    gsub("^[0-9]", "leadingNumber", colnames(target_labeled_data))
  colnames(target_labeled_data) <- gsub("\\*", "star", colnames(target_labeled_data))
  
  return(target_labeled_data)
}