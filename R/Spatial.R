#' @import Seurat ggplot2 Matrix
#' @importFrom compositions geometricmean geometricmeanRow 
#' @importFrom RUVSeq RUVg
#' @importFrom matrixStats rowMedians

#' @title Q3_Normalization
#'
#' @description a normalization method that uses the 75th percentile to find size factors for data normalization.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param targetAssayName The assay name that the normalized data should be stored in. 
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the Q3 assay. 
#' @export

Q3_Normalization <- function(seuratObj, assay = "RNA", targetAssayName = "Q3"){
  #from here: https://www.dlongwood.com/wp-content/uploads/2021/05/WP_MK2833_Introduction_to_GeoMx_Normalization_CTA_R5-1.pdf page 7
  # first, divide all the genes per AOI by their respective Q3 count, 
  #ii) second, multiply all the genes in all AOIs by a constant, defined as the geometric mean of Q3 counts for all AOIs.
  
  #calculate column (ROI) Q3s
  counts <- as.matrix(seuratObj[[assay]]@counts)
  Q3_Sample_Counts <- apply(counts,MARGIN = 2,FUN = quantile,prob=c(.75))
  #Ensure no downstream division by zero errors
  if(any(Q3_Sample_Counts == 0)){
    stop("Error: 75th percentile is zero for some samples. Q3 is either too low to normalize by, or the data is too sparse.")
  }
  #divide count matrix by vector of Q3s
  Q3_Divided_Counts <- counts/Q3_Sample_Counts
  #multiply Q3 divided counts by geometric mean of Q3s (scalar)
  Q3_Final_Counts <- Q3_Divided_Counts * compositions::geometricmean(Q3_Sample_Counts)
  
  #stash into seurat object under new assay
  seuratObj[[targetAssayName]] <- Seurat::CreateAssayObject(methods::as(Q3_Final_Counts, "dgCMatrix"))
  
  return(seuratObj)
}

#' @title RUVg_Housekeeping_Normalization
#'
#' @description a normalization method that uses residuals from a fitted model to Remove Unwanted Variation by leveraging housekeeping genes (https://www.nature.com/articles/nbt.2931).
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param targetAssayName The assay name that the normalized data should be stored in. 
#' @param k the number of factors of unwanted variation. See section "Choice of number of factors of unwanted variation." in paper cited in the description.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the RUVg assay. 
#' @export

RUVg_Housekeeping_Normalization <- function(seuratObj, assay = "RNA", targetAssayName = "RUVg",  k = 1){
  housekeeping_genes <- c("ABCF1", "ACTB", "ALAS1", "B2M", "CLTC", "G6PD", "GAPDH", 
                          "GUSB", "HPRT1", "LDHA", "PGK1", "POLR1B", "POLR2A",
                          "RPL19", "RPLP0", "SDHA", "TBP", "TUBB")
  counts <- as.matrix(seuratObj[[assay]]@counts)
  rounded_counts <- round(as.matrix(counts))
  RUVg <- RUVSeq::RUVg(rounded_counts,
                       housekeeping_genes[housekeeping_genes %in% rownames(seuratObj)], 
                       k = k)
  seuratObj[[targetAssayName]] <- Seurat::CreateAssayObject(methods::as(RUVg$normalizedCounts, "dgCMatrix"))
  return(seuratObj)
}

#' @title NanoString_Housekeeping_Normalization
#'
#' @description a normalization method that uses the housekeeping genes to generate size factors for normalization.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param targetAssayName The assay name that the normalized data should be stored in. 
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the "Housekeeping" assay. 
#' @export

NanoString_Housekeeping_Normalization <- function(seuratObj, assay = "RNA", targetAssayName = "Housekeeping"){
  
  #Normalization Steps from Nanostring 
  #1.  Calculate the geometric mean of the selected housekeeping genes for each lane.
  #2. Calculate the arithmetic mean of these geometric means for all sample lanes.
  #3. Divide this arithmetic mean by the geometric mean of each lane to generate a lane-specific normalization factor.
  #4. Multiply the counts for every gene by its lane-specific normalization factor.
  
  counts <- as.matrix(seuratObj[[assay]]@counts) 
  
  housekeeping_genes <- c("ABCF1", "ACTB", "ALAS1", "B2M", "CLTC", "G6PD", "GAPDH", 
                          "GUSB", "HPRT1", "LDHA", "PGK1", "POLR1B", "POLR2A",
                          "RPL19", "RPLP0", "SDHA", "TBP", "TUBB")
  #Check to ensure housekeeping genes present in the count matrix.
  if (sum(rownames(counts) %in% housekeeping_genes) == 0){
    stop("No housekeeping genes detected in the count matrix.")
  }
  #Report how many housekeeping genes are present in the count matrix.
  print(paste0(sum(rownames(counts) %in% housekeeping_genes), " housekeeping genes detected in the count matrix."))
  
  housekeeping_counts <- counts[rownames(counts) %in% housekeeping_genes, ]
  
  #Ensure no downstream division by zero errors
  if(any(compositions::geometricmeanRow(as.matrix(t(housekeeping_counts)) == 0))){
    stop("Error: Geometric mean of housekeeping counts for one or more samples is zero. One of the housekeeping genes likely has a value of zero. Either add a pseudocount (if appropriate) or choose a different normalization")
  }
  
  #steps 1-3 of the normalization are computed here
  Sample_Normalization_Factors <- mean(compositions::geometricmeanRow(as.matrix(t(housekeeping_counts)))) / compositions::geometricmeanRow(as.matrix(t(housekeeping_counts)))
  
  #step 4 of the normalization is computed here
  hk_normalized_counts <- housekeeping_counts %*% diag(Sample_Normalization_Factors)
  colnames(hk_normalized_counts) <- colnames(counts)
  seuratObj[[targetAssayName]] <- Seurat::CreateAssayObject(counts = methods::as(hk_normalized_counts, "dgCMatrix"))
  return(seuratObj)
}

#' @title SpatialNormalizeByGroup
#'
#' @description a wrapper function to iteratively apply normalization methods that may be sensitive to batch variation.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param normalizationMethod The normalization method to be applied (One of: Q3, Housekeeping, or RUVg) across the batches/groups in groupField. 
#' @param groupField the metadata column that delineates how the samples within Seurat object should be grouped. 
#' @param targetAssayName The assay name that the normalized data should be stored in. 
#' @param inferAssayName Boolean to determine whether or not the function should infer the name of the assay that the normalized counts will be stored in based on the name of normalizationMethod.
#' @param k a passthrough variable for the RUVg normalization. Please see RUVg_Housekeeping_Normalization for details about this parameter.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the "Housekeeping" assay. 
#' @export

SpatialNormalizeByGroup <- function(seuratObj, assay = "RNA", normalizationMethod = NULL, inferAssayName = T, targetAssayName = NULL, groupField = "SlideName", k = NA){
  #Test that a normalizationMethod is set.
  if (is.null(normalizationMethod)){
    stop("Please specify a normalization method. Currently supported methods are: Q3, Housekeeping, or RUVg.")
  }
  #Ensure a resolvable targetAssayName, and infer targetAssayName if desired.
  if (inferAssayName != TRUE & is.null(targetAssayName)){
    stop("Please either set a targetAssayName or set inferAssayName = TRUE.")
  } else if (inferAssayName == TRUE & !is.null(targetAssayName)){
    print(paste0("inferAssayName is TRUE. Setting targetAssayName = ", normalizationMethod, "."))
    targetAssayName <- normalizationMethod
  } else if (inferAssayName != TRUE & !is.null(targetAssayName)){
    print(paste0("Normalizations will be stored in the ", targetAssayName, " assay."))
  }
    
  
  
  #split the seurat object according to groupField
  list_of_seurat_objects <- CellMembrane::SplitSeurat(seuratObj, splitField = groupField, minCellsToKeep = 1)
  print(paste("Normalization method:", normalizationMethod))
  #Apply one of the Normalization methods across the values of groupField
  for (group in names(list_of_seurat_objects)){
    if (normalizationMethod == "Q3"){
      list_of_seurat_objects[[group]] <- Q3_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay, targetAssayName = targetAssayName)
    } else if(normalizationMethod == "Housekeeping") {
      list_of_seurat_objects[[group]] <- NanoString_Housekeeping_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay, targetAssayName = targetAssayName)
    } else if(normalizationMethod == "RUVg"){
      list_of_seurat_objects[[group]] <- RUVg_Housekeeping_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay, targetAssayName = targetAssayName)
    } else{
      #Test for a supported normalizationMethod.
      stop("Error: no supported normalization method selected. Please select one of: Q3, Housekeeping, or RUVg.")
    }
  }
  
  print(paste("Merging", length(names(list_of_seurat_objects)), "Seurat Objects."))
  #instantiate the merged seuratObj
  mergedSeuratObj <- list_of_seurat_objects[[names(list_of_seurat_objects)[[1]]]]
  #iteratively merge the seurat objects 
  for (group in names(list_of_seurat_objects)[-1]){
    mergedSeuratObj <- merge(mergedSeuratObj, list_of_seurat_objects[[group]])
  }
  return(mergedSeuratObj)
}

#' @title RLE_Plot
#'
#' @description a function to compute relative log expression plots to qualitatively test whether or not a normalization successfully normalized the data.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param sampleIdentifier The metadata column that defines individual samples. 
#' @param colorVariable The metadata column that determines the color of the boxplots of RLE. Optimally, this should be some type of phenotypic description.
#' @return Returns a relative log expression plot
#' @export

RLE_Plot <- function(seuratObj, assay = "RNA", sampleIdentifier = "SegmentDisplayName", colorVariable){
  metadata <- seuratObj@meta.data
  
  #Calculate a dataframe of relative log expression. 
  deviations <- as.matrix(log2(seuratObj[[assay]]@counts+1)) - matrixStats::rowMedians(as.matrix(log2(seuratObj[[assay]]@counts+1)))
  deviations <- dplyr::as_tibble(as.data.frame(deviations))
  deviations <- tidyr::gather(deviations)
  colnames(deviations) <- c("Sample", "RLE")
  
  #merge the deviations with the metadata to color the final plot by colorVariable
  deviations_with_metadata <- merge(deviations, metadata, by.x = "Sample", by.y = sampleIdentifier )
  
  plot <- ggplot(deviations_with_metadata, aes(x = Sample, y = RLE, fill = !!sym(colorVariable))) +
    geom_boxplot(outlier.shape = NA, lwd = 0.5) + 
    theme_bw()+  
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) + 
    ggtitle(paste0(assay," RLE")) + 
    scale_color_manual(values = "gray90")
  
  #if the data is numeric, color it with a 4 scale "Low-MediumLow-MediumHigh-High" color scheme. 
  if(is.numeric(metadata[,colorVariable])){
    plot <- plot + scale_fill_gradientn(colors = c("navy", "dodgerblue", "gold", "red"))
  }
  return(plot)
}


