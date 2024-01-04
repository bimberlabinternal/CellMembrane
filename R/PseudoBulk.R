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
#' @param metaFieldCollapseCharacter The character to use when concatenating metadata fields together to form the sample key field
#' @param nCountRnaStratification A boolean determining whether or not automatic outlier detection of clusters with abnormal nCount_RNA should be detected.
#' @param stratificationGroupingFields The metadata field containing groupings of cells to be specified as having abnormal nCount_RNA distributions or not. 
#' @return An aggregated Seurat object.
#' @export
PseudobulkSeurat <- function(seuratObj, 
                             groupFields, 
                             assays = NULL, 
                             additionalFieldsToAggregate = NULL, 
                             metaFieldCollapseCharacter = '|', 
                             nCountRnaStratification = F, 
                             stratificationGroupingFields = c("ClusterNames_0.2", "ClusterNames_0.4", "ClusterNames_0.6", "ClusterNames_0.8", "ClusterNames1.2")) {
  if (!all(groupFields %in% names(seuratObj@meta.data))) {
    stop('All fields from groupFields must be in seuratObj@meta.data')
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
      if (length(abnormal_RNA_clusters_index) == 0){
        abnormal_RNA_clusters_index <- NA
      } else {
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
    return(paste0(make.names(y), collapse = metaFieldCollapseCharacter))
  }))
  
  Seurat::Idents(seuratObj) <- seuratObj$KeyField
  
  # This generates the mean() of counts. Even though we want sum(), this is a convenient way to ensure all other
  a <- Seurat::AverageExpression(seuratObj, group.by = 'KeyField', return.seurat = T, verbose = F, slot = 'counts', assays = assays)
  
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
    stop('The keyField and AverageExpression object keys to do not match')
  }
  
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


#TODO: add support to extend this to propeller (pseudobulk the object, then optionally discard.)
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

#' @title FilterPseudobulkContrasts
#'
#' @description This is designed to accept a study design defined by a series of logical gates applied to metadata fields, test them for equivalence to some criterion, and then filter pairwise contrasts accordingly.
#' @param logic_list The list that defines the study design. Please see examples for the setup and use. This list be a list of lists where each sub-list has three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor. The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.
#' @param design a design/model matrix returned by DesignModelMatrix().
#' @param use_require_identical_logic Whether or not to apply require_identical logic using require_identical_fields. It is possible to use logic gates to keep metadata fields constant, but using this saves time & effort/iteration.
#' @param require_identical_fields The metadata columns of the SeuratObj that you wish to keep constant. For example: defining "Tissue" in this vector would filter all contrasts comparing Liver to PBMC, as their tissue is not identical.
#' @param filtered_contrasts_output_file The file to write the list of filtered contrasts to (one pair of samples per row, separated into their metadata fields).
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
#' logic_list <- list(
#'                     list('Challenge', 'xor', 'Mock.challenged'), 
#'                     list('SampleType', 'all', 'Necropsy')
#'             )
#'
#' 
#' #Finally, enumerate all possible contrasts, then filter and return them according to the study design defined by logic_list.
#' #If you want to perform "Like-vs-Like" contrasts (i.e. T cell vs T cell within a cell_type metadata field), specify those columns as required to be identical and set use_require_identical_logic to TRUE.
#' filtered_contrasts <- FilterPseudobulkContrasts(logic_list = logic_list, 
#'                                                 design = design, 
#'                                                 use_require_identical_logic = T, 
#'                                                 require_identical_fields = c('Tissue','cell_type', 'SampleType'), 
#'                                                 filtered_contrasts_output_file = './filtered_contrasts.tsv')
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

FilterPseudobulkContrasts <- function(logic_list = NULL, design = NULL, use_require_identical_logic = T, require_identical_fields = NULL, filtered_contrasts_output_file = './filtered_contrasts.tsv'){
  #check design matrix.
  if (is.null(design)){
    stop("Please define a design matrix. The design matrix is intended to be returned by DesignModelMatrix().")
  }
  #ensure design matrix is from DesignModelMatrix() (or at least has contrast_columns in its attributes.)
  if (is.null(attr(design, "contrast_columns"))){
    stop("'contrast_columns' is missing from the design model matrix attributes. Please use CellMembrane's DesignModelMatrix to create the design matrix, or stash the contrast_columns vector into the design matrix using \n  attr(design, 'contrast_columns') <- contrast_columns ")
  }
  contrast_columns <- attr(design, "contrast_columns")
  #ensure logic_list is properly defined.
  if (is.null(logic_list)){
    stop("Please supply a list that defines the logic gates used for filtering the contrasts. This list should have three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (typeof(logic_list) != "list"){
    stop("Please ensure that logic_list is a list composed of lists. Each sub-list should define a gate that will filter possible contrasts. This list should have three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (!(all(lengths(logic_list) == 3))){
    stop("The lengths of all of the elements (i.e. logic gates) within logic_list are not equal to 3. Please ensure there are exactly three entries in each element of logic_list. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor). The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.")
  } else if (!(all(sapply(logic_list, "[[", 1) %in% contrast_columns))){
    stop("There are elements (logic gates) within logic_list whose first element is not in the contrast_columns found within the design matrix. contrast_columns is the vector supplied to DesignModelMatrix() upstream of this function and holds the metadata field names of the pseudobulked Seurat object that comprise the contrasts. If you're seeing this, something has gone awry (like using an old design matrix or a typo in the logic_list), or you're trying to supply your own design/model matrix, which is currently not supported.")
  } else if (any(sapply(logic_list, "[[", 2) %in% c("and", "AND", "or", "OR"))){
    stop("Error in one of the logic gate specifications. For an AND gate, please use the function name 'all'. For an OR gate, please use the function name 'any'.")
  } else if (!all(tolower(sapply(logic_list, "[[", 2)) %in% c("any", "all", "xor", "nand", "nor", "xnor"))){
    stop("Error in one of the logic gate specifications. Please use one of: 'any', 'all', 'xor', 'nand', 'nor', 'xnor' to specify your logic gate.")
  }
  
  #check require_identical fields arguments.
  if (!is.logical(use_require_identical_logic)){
    stop("Please set use_require_identical_logic to TRUE or FALSE.")
  }
  if (!is.null(require_identical_fields) & !use_require_identical_logic){
    warning("use_require_identical_logic is set to FALSE, but require_identical_fields are supplied. These will not be used for filtering.")
  } else if (use_require_identical_logic & is.null(require_identical_fields)){
    stop("use_require_identical_logic is TRUE, but no require_identical_fields have been supplied.")
  } else if (!all(require_identical_fields %in% contrast_columns)){
    stop("Metadata columns requested in require_identical_fields do not appear in the supplied design matrix's contrast_columns attribute. Please ensure that the contrast_columns vector supplied to DesignModelMatrix() is correct and that all elements in the require_identical_fields vector are metadata column names whose values you intended to keep constant across all pairwise contrasts.")
  }
  
  #create a nx2 array of all possible unique pairwise combinations of contrasts from the design matrix.
  contrasts <- t(utils::combn(colnames(design), m = 2))
  
  #initialize contrast indices vector (to be used to subset the total contrast list after filtering).
  filtered_contrast_indices <- c()
  #iterate through all of the contrasts, filter those that don't satisfy the logic gates/require_identical_fields logic.
  for (row_index in 1:nrow(contrasts)){
    #define the relevant metadata fields for each side of the contrast.
    positive_contrast <- unlist(strsplit(contrasts[row_index,1], split = "_"))
    names(positive_contrast) <- contrast_columns
    negative_contrast <- unlist(strsplit(contrasts[row_index,2], split = "_"))
    names(negative_contrast) <- contrast_columns
    
    #define iterator to track how many logic gates were passed for the contrast pair.
    gates_satisfied <- 0
    for (logic_gate_index in 1:length(logic_list)){
      #define variables to set up the logic gate from logic_list.
      field_to_check <- logic_list[[logic_gate_index]][[1]]
      gate <- logic_list[[logic_gate_index]][[2]]
      criterion <- logic_list[[logic_gate_index]][[3]]
      #set up logic gate and check if the contrast satisfies the gate.
      if (get(gate)(
        positive_contrast[field_to_check] == criterion,
        negative_contrast[field_to_check] == criterion)){
        #if true, mark the gate as satisfied.
        gates_satisfied <- gates_satisfied + 1
        #it is possible to fully define your study design within logic_gate without using the require_identical_fields logic as a shortcut, so we check if it's used or not.
        if (use_require_identical_logic == T){
          #check if the contrast passed all of the logic gates, if so, test for require_identical fields.
          if (gates_satisfied == length(logic_list)){
            #reset counter for matching require_identical fields.
            require_identical_fields_satisfied <- 0
            for (require_identical in require_identical_fields){
              #keep track of how many require_identical fields are equal to each other.
              if (positive_contrast[require_identical] == negative_contrast[require_identical]){
                require_identical_fields_satisfied <- require_identical_fields_satisfied + 1
              }
              #if all of the require_identical fields are require_identical keep the contrast.
              if (require_identical_fields_satisfied == length(require_identical_fields)){
                filtered_contrast_indices <- c(filtered_contrast_indices, row_index)
              }
            }
          }
        } else{
          #if all of the gates are satisfied and the require_identical logic isn't used, keep the contrast.
          if (gates_satisfied == length(logic_list)){
            filtered_contrast_indices <- c(filtered_contrast_indices, row_index)
          }
        }
      }
    }
  }
  #filter the contrasts using the indices that satisfy the logic gates/require_identical fields.
  contrasts <- contrasts[filtered_contrast_indices,]
  
  #initialize a dataframe to store the filtered contrasts.
  filtered_contrasts_dataframe <- data.frame()
  #iterate through the filtered contrasts and coerce the matrix into a data frame. 
  for(contrast_index in 1:nrow(contrasts)){
    #split the contrasts by _ and coerce into a dataframe. 
    temporary_data_frame <- t(as.data.frame(unlist(strsplit(contrasts[contrast_index,], split = "_"))))
    #assign the column names of the dataframe via the ordering supplied by contrast_columns.
    colnames(temporary_data_frame) <- c(paste0("positive_contrast_", contrast_columns), 
                                        paste0("negative_contrast_", contrast_columns))
    rownames(temporary_data_frame) <- contrast_index
    filtered_contrasts_dataframe <- rbind(filtered_contrasts_dataframe, temporary_data_frame)
  }
  #write the filtered output file in "data.frame" format, with column names and no row indices. 
  write.table(filtered_contrasts_dataframe, file = filtered_contrasts_output_file, row.names = F, col.names = T)
  return(filtered_contrasts_dataframe)
}

#logic gates to be used with FilterPseudobulkContrasts()
nand <- function(x,y){ !(x&y) } 
nor <- function(x,y){ !(x|y) } 
xnor <- function(x,y){(x == y)}

#' @title FitRegularizedClassificationGlm 
#'
#' @description Treating gene expression like a classification problem, this function trains a penalized model to classify a metadata feature. 
#' @param seuratObj a Seurat object
#' @param metadataVariableForClassification The metadata feature to be classified. If non-binary, then multinomial regression will automatically be performed. 
#' @param rescale The feature selection will optimize for "heatmap-interpretable genes" so the features are intended to be scaled. If TRUE, this will rescale the variable features.
#' @param numberOfVariableFeatures A parameter to select how many features should be selected as variable for scaling, by default, all genes will be used. 
#' @param assay Seurat Object's assay
#' @param slot Slot within the Seurat object assay. Recommended to be "scale.data".
#' @param devianceCutoff Tolerance for model error when deciding how much regularization should be performed. 1 = no tolerance for error, 0 = intercept only, no genes used for prediction.
#' @param split the option to provide a previous model's training/testing set. This is necessary if you're performing multiple iterations of model fitting. 
#' @param returnModelAndSplits A boolean option to return a list containing the fitted model and training/testing splits in addition to the useful features. 
#' @return A vector of genes useful for classification and, optionally, the model and training/testing sets.
#' @import mlr3verse
#' @importFrom mlr3 as_task_classif lrn partition
#' @export

FitRegularizedClassificationGlm <- function(seuratObj,
                                            metadataVariableForClassification = NULL,
                                            rescale = TRUE,
                                            numberOfVariableFeatures = 3000,
                                            assay = "RNA",
                                            slot = "scale.data",
                                            devianceCutoff = 0.8,
                                            split = NULL, 
                                            returnModelAndSplits = F) {
  #sanity check arguments
  if (is.null(metadataVariableForClassification)){
    stop("Please supply a column of the seurat object's metadata to classify.")
  }
  if (!(metadataVariableForClassification %in% colnames(seuratObj@meta.data))){
    stop("Supplied metadataVariableForClassification not found in the seurat object's metadata. Please ensure your metadata column is spelled correctly and exists in seuratObj@meta.data.")
  }
  # rescale the input data (in case of an upstream subset since it was last rescaled).
  if (rescale) {
    if (is.null(numberOfVariableFeatures)) {
      seuratObj <- CellMembrane::NormalizeAndScale(seuratObj, 
                                                   nVariableFeatures = length(rownames(seuratObj)),
                                                   variableGenesBlacklist = RIRA::GetGeneSet("VariableGenes_Exclusion.2"),
                                                   scoreCellCycle=F
      )
    } else {
      seuratObj <- CellMembrane::NormalizeAndScale(seuratObj, 
                                                   nVariableFeatures = numberOfVariableFeatures,
                                                   variableGenesBlacklist = RIRA::GetGeneSet("VariableGenes_Exclusion.2"),
                                                   scoreCellCycle=F
      )
    }
  }
  #convert the scale.data matrix to include a labeled classification column
  target_labeled_data <-
    #merge the seuratObj's requested slot (converted to dense just in case a non-scale.data slot was used)
    merge(
      Matrix::t(as.matrix(Seurat::GetAssayData(
        seuratObj,
        assay = assay,
        slot = slot
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
  
  ##set up task
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
