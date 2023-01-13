#' @import Seurat ggplot2
#' @importFrom compositions geometricmean geometricmeanRow 
#' @importFrom RUVSeq RUVg
#' @importFrom matrixStats rowMedians

#' @title Q3_Normalization
#'
#' @description a normalization method that uses the 75th percentile to find size factors for data normalization.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the Q3 assay. 
#' @export
Q3_Normalization <- function(seuratObj, assay = "RNA"){
  #from here: https://www.dlongwood.com/wp-content/uploads/2021/05/WP_MK2833_Introduction_to_GeoMx_Normalization_CTA_R5-1.pdf page 7
  # first, divide all the genes per AOI by their respective Q3 count, 
  #ii) second, multiply all the genes in all AOIs by a constant, defined as the geometric mean of Q3 counts for all AOIs.
  
  #calculate column (ROI) Q3s
  counts <- seuratObj[[assay]]@counts
  Q3_Sample_Counts <- apply(counts,MARGIN = 2,FUN = quantile,prob=c(.75))
  #divide count matrix by vector of Q3s
  Q3_Divided_Counts <- counts/Q3_Sample_Counts
  #multiply Q3 divided counts by geometric mean of Q3s (scalar)
  Q3_Final_Counts <- Q3_Divided_Counts * compositions::geometricmean(Q3_Sample_Counts)
  
  #stash into seurat object under new assay
  seuratObj[['Q3']] <- Seurat::CreateAssayObject(Q3_Final_Counts)
  
  return(seuratObj)
}

#' @title RUVg_Housekeeping_Normalization
#'
#' @description a normalization method that uses residuals from a fitted model to Remove Unwanted Variation by leveraging housekeeping genes (https://www.nature.com/articles/nbt.2931).
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param k the number of factors of unwanted variation. See section "Choice of number of factors of unwanted variation." in paper cited in the description.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the RUVg assay. 
#' @export
RUVg_Housekeeping_Normalization <- function(seuratObj, assay = "RNA", k = 1){
  housekeeping_genes <- c("ABCF1", "ACTB", "ALAS1", "B2M", "CLTC", "G6PD", "GAPDH", 
                          "GUSB", "HPRT1", "LDHA", "PGK1", "POLR1B", "POLR2A",
                          "RPL19", "RPLP0", "SDHA", "TBP", "TUBB")
  counts <- seuratObj[[assay]]@counts
  rounded_counts <- round(as.matrix(counts))
  RUVg <- RUVSeq::RUVg(rounded_counts,
                       housekeeping_genes[housekeeping_genes %in% rownames(seuratObj)], 
                       k = k)
  seuratObj[['RUVg']] <- Seurat::CreateAssayObject(RUVg$normalizedCounts)
  return(seuratObj)
}



#' @title NanoString_Housekeeping_Normalization
#'
#' @description a normalization method that uses the housekeeping genes to generate size factors for normalization.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the "Housekeeping" assay. 
#' @export

NanoString_Housekeeping_Normalization <- function(seuratObj, assay = "RNA"){
  
  #Normalization Steps from Nanostring 
  #1.  Calculate the geometric mean of the selected housekeeping genes for each lane.
  #2. Calculate the arithmetic mean of these geometric means for all sample lanes.
  #3. Divide this arithmetic mean by the geometric mean of each lane to generate a lane-specific normalization factor.
  #4. Multiply the counts for every gene by its lane-specific normalization factor.
  
  counts <- seuratObj[[assay]]@counts 
  
  housekeeping_genes <- c("ABCF1", "ACTB", "ALAS1", "B2M", "CLTC", "G6PD", "GAPDH", 
                          "GUSB", "HPRT1", "LDHA", "PGK1", "POLR1B", "POLR2A",
                          "RPL19", "RPLP0", "SDHA", "TBP", "TUBB")
  housekeeping_counts <- counts[rownames(counts) %in% housekeeping_genes, ]
  
  #steps 1-3 of the normalization are computed here
  Sample_Normalization_Factors <- mean(compositions::geometricmeanRow(as.matrix(t(housekeeping_counts)))) / compositions::geometricmeanRow(as.matrix(t(housekeeping_counts)))
  #step 4 of the normalization is computed here
  hk_normalized_counts <- housekeeping_counts %*% diag(Sample_Normalization_Factors)
  colnames(hk_normalized_counts) <- colnames(counts)
  seuratObj[["Housekeeping"]] <- Seurat::CreateAssayObject(counts = hk_normalized_counts)
  return(seuratObj)
}

#' @title SpatialNormalizeByGroup
#'
#' @description a wrapper function to iteratively apply normalization methods that may be sensitive to batch variation.
#' @param seuratObj The Seurat object that holds the data to be normalized.
#' @param assay The assay holding the raw data - the counts slot will be assumed to hold raw data.
#' @param normalizationMethod The normalization method to be applied (One of: Q3, Housekeeping, or RUVg) across the batches/groups in groupField. 
#' @param groupField the metadata column that delineates how the samples within Seurat object should be grouped. 
#' @param k a passthrough variable for the RUVg normalization. Please see RUVg_Housekeeping_Normalization for details about this parameter.
#' @return Returns a Seurat object containing normalized counts stored in the counts slot of the "Housekeeping" assay. 
#' @export

SpatialNormalizeByGroup <- function(seuratObj, assay = "RNA", normalizationMethod= NULL, groupField = "SlideName", k = NA){
  #Ensure that groupField is the currently set Idents value.
  Seurat::Idents(seuratObj) <- groupField
  #Test that a normalizationMethod is set.
  if (is.null(normalizationMethod)){
    stop("Please specify a normalization method. Currently supported methods are: Q3, Housekeeping, or RUVg.")
  }
  
  #split the seurat object according to groupField
  list_of_seurat_objects <- CellMembrane::SplitSeurat(seuratObj, splitField = groupField, minCellsToKeep = 1)
  print(paste("Normalization method:", normalizationMethod))
  #Apply one of the Normalization methods across the values of groupField
  for (group in names(list_of_seurat_objects)){
    if (normalizationMethod == "Q3"){
      list_of_seurat_objects[[group]] <- Q3_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay)
    } else if(normalizationMethod == "Housekeeping") {
      list_of_seurat_objects[[group]] <- NanoString_Housekeeping_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay)
    } else if(normalizationMethod == "RUVg"){
      list_of_seurat_objects[[group]] <- RUVg_Housekeeping_Normalization(seuratObj = list_of_seurat_objects[[group]], assay = assay)
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


