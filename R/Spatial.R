#' @import Seurat ggplot2 

utils::globalVariables(
  names = c('Sample', 'RLE'),
  package = 'CellMembrane',
  add = TRUE
)

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
  counts <- as.matrix(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts'))
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
  counts <- as.matrix(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts'))
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
  
  counts <- as.matrix(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts'))
  
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
  if(any(compositions::geometricmeanRow(as.matrix(t(housekeeping_counts))) == 0)){
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
  deviations <- as.matrix(log2(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts')+1)) - matrixStats::rowMedians(as.matrix(log2(Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts')+1)))
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

#CosMx Functions 

#' @title DetectCellStructuresBasedOnCellType
#'
#' @description A DBSCAN-based spatial cluster detection method to find dense cellular structures within spatial data. This function is currently implemented with FOV-based, iterative, classification in mind, but is technically extendable to global analysis of a single frame. 
#' @param seuratObjectMetadata Metadata dataframe storing, at minimum, a cell type or "transcriptomic clustering" type field, a field denoting FOV, and two spatial coordinates. 
#' @param cellTypeField A character field storing cell type annotations, but could be any discrete cell clustering assignment. 
#' @param minimumClusterSizeCoefficient The percentage of the input data (after cell type splitting) to be considered as a feasible minimum cluster size. Determines how many "noise data points" will be added by DBSCAN.  
#' @param fovField the metadata column that stores the Field of View information.
#' @param fovWhitelist An optional whitelist of FOVs. By default, the function will loop over all FOVs, which could be time consuming. 
#' @param cellTypeWhiteList A vector of genes that constitute your substructure (e.g. c("Bcell", "BCell") for B cell follicles).
#' @param xCoordinateField The metadata column that stores the x coordinate information within the Field of View. 
#' @param yCoordinateField The metadata column that stores the y coordinate information within the Field of View. 
#' @param substructureMetaDataFieldName An annotation that will be concatenated during the results. "Local" FOV information will be concatenated using "fov + substructureMetaDataFieldName + a substructure index" within the columns of the metadata. 
#' @param summarizeLocalResults An optional boolean that will wrap up the various substructureMetaDataFieldName columns into two single columns. One, which determines if a cell is within ANY of the defined substructures, stored in the output column "Within_Local + substructureMetaDataFieldName". The second is a metadata column that displays which of the local substructures the cell belongs in, concatenated as "Local + substructureMetaDataFieldName". "Local + substructureMetaDataFieldName + 0" is always the noise designation. 
#' @return Returns a dataframe containing columns related to the substructures found within the images at varying scopes. With summarizeLocalResults = FALSE, (number of FOVs) x (number of subtructures + 1) columns will be added. summarizeLocalResults rolls these high resolution results into two additional columns relative to the fovField. 
#' @examples
#' \dontrun{
#' #Perform Cell Structure detection for B cell follicles. 
#' 
#' metadata <- DetectCellStructuresBasedOnCellType(seuratObjectMetadata, 
#' cellTypeField = "cell_type", 
#' minimumClusterSizeCoefficient = 0.05,
#' fovField = "fov",
#' fovWhitelist = 1,
#' cellTypeWhiteList = c("Bcell", "B_cell", "B.cell"),
#' xCoordinateField = "x_FOV_px", 
#' yCoordinateField = "y_FOV_px", 
#' substructureMetaDataFieldName = "BCF",
#' summarizeLocalResults = TRUE
#' )
#' 
#' #load/install packages for plotting
#' library(pacman)
#' p_load(ggplot2, dplyr, egg, patchwork)
#' 
#' #define plotting layout
#' layout <- "
#' #AAAA#
#' BBBCCC
#' "
#' 
#' #Plot results
#' ggplot(metadata %>% filter(fov == 1), 
#' aes(x = x_FOV_px, y = y_FOV_px, color = simple_cellType)) + 
#' geom_point() + 
#' egg::theme_article() + 
#' ggtitle('Cell Type assignment')
#' ggplot(metadata %>% filter(fov == 1), aes(x = x_FOV_px, y = y_FOV_px, color = factor(Local_BCF))) + 
#' geom_point() + 
#' egg::theme_article() +  
#' ggtitle('Specific substructure cell assignment') + 
#' ggplot(metadata %>% filter(fov == 1), 
#' aes(x = x_FOV_px, y = y_FOV_px, color = Within_Local_BCF)) + 
#' geom_point() + 
#' egg::theme_article() + 
#' ggtitle('Non-specific substructure cell assignment') + 
#' plot_layout(design = layout, guides = "collect")
#' }
#' @export
DetectCellStructuresBasedOnCellType <- function(seuratObjectMetadata, 
                                                cellTypeField = "cell_type", 
                                                minimumClusterSizeCoefficient = 0.05,
                                                fovField = "fov",
                                                fovWhitelist = NULL,
                                                cellTypeWhiteList = c("Bcell", "B_cell", "B.cell"),
                                                xCoordinateField = "x_FOV_px", 
                                                yCoordinateField = "y_FOV_px", 
                                                substructureMetaDataFieldName = "BCF",
                                                summarizeLocalResults = TRUE
){
  #check that inputs exist and summarizeLocalResults is a boolean. 
  if (!(all(c(cellTypeField, fovField, xCoordinateField, yCoordinateField) %in% colnames(seuratObjectMetadata)))) {
    missingColumn <- c(cellTypeField, fovField, xCoordinateField, yCoordinateField)[which(!(c(cellTypeField, fovField, xCoordinateField, yCoordinateField) %in% colnames(seuratObjectMetadata)))]
    stop(paste0("Error: ", missingColumn , " was not found in the supplied metadata's column names. Please add ", missingColumn, " to the metadata."))
  } else if (!rlang::is_bool(summarizeLocalResults)) {
    stop("Please set summarizeLocalResults to either TRUE or FALSE.")
  }
  
  #determine FOVs to loop over
  if (is.null(fovWhitelist)) {
    fovs <- unique(seuratObjectMetadata[,fovField])
  } else if (!fovWhitelist %in% unique(seuratObjectMetadata[,fovField])){
    stop(paste0("The FOVs listed in whitelist are not detected in fovField: ", fovField, ". Please ensure the fovField is formatted as you expect (specifically: ensure integer vs character typing)."))
  } else {
    fovs <- fovWhitelist
  }
  
  #convert cellTypeWhitelist to regex
  #escape special characters in cell types
  escapedCellTypes <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", cellTypeWhiteList)
  cellTypeConstituentRegex <- paste0(escapedCellTypes, collapse = "|")
  
  #Detect substructures within each FOV
  for (fov in fovs) {
    fov_seuratObjectMetadata <- seuratObjectMetadata[seuratObjectMetadata[fovField] == fov,]
    #use regex to subset to the whitelisted cell type for substructure clustering
    coordinateDataFrameOfCellsOfInterest <- fov_seuratObjectMetadata[grepl(cellTypeConstituentRegex, fov_seuratObjectMetadata[,cellTypeField]), c(xCoordinateField, yCoordinateField)]
    #test to ensure regex matched cells within the current FOV. 
    if (nrow(coordinateDataFrameOfCellsOfInterest) > 0) {
      #cluster cells
      clusters <- dbscan::hdbscan(coordinateDataFrameOfCellsOfInterest, 
                                  minPts = nrow(coordinateDataFrameOfCellsOfInterest) * minimumClusterSizeCoefficient)
      #store cluster assignments for cells within the cell type of interest
      fov_seuratObjectMetadata[grepl(cellTypeConstituentRegex, fov_seuratObjectMetadata[,cellTypeField]), "dbscan_cluster"] <- clusters$cluster
      #iterate over clusters (e.g. dense structures of cell types)
      for (cluster_index in unique(fov_seuratObjectMetadata$dbscan_cluster)) {
        #cluster 0 is the "noise cluster" from DBSCAN
        if (cluster_index != 0 & !is.na(cluster_index)) {
          print(paste0("Identifying substructure ", cluster_index, " in FOV:", fov))
          #isolate cells from the current cluster/substructure
          coordinateDataFrameOfSubStructure <- fov_seuratObjectMetadata[
            grepl(cellTypeConstituentRegex, fov_seuratObjectMetadata[,cellTypeField]) & 
              fov_seuratObjectMetadata$dbscan_cluster == cluster_index, 
            c(xCoordinateField, yCoordinateField)]
          #form convex hull around substructure
          convexHull <- geometry::convhulln(coordinateDataFrameOfSubStructure, output.options = TRUE)
          #compare each cell and classify them as internal or external to the convex hull (in or out of the structure)
          mat <- as.matrix(fov_seuratObjectMetadata[, c(xCoordinateField, yCoordinateField)])
          #inhulln is picky about the formatting of the input matrix, so we need to null the col/rownames.
          colnames(mat) <- NULL
          rownames(mat) <- NULL
          #inhull_yes_no is a boolean vector determining if the point (cell centroid) is internal to the substructure or not. 
          inhull_yes_no <- geometry::inhulln(ch = convexHull, 
                                             p = cbind(as.double(mat[,1]), 
                                                       as.double(mat[,2])))
          
          #instantiate field as NULL
          seuratObjectMetadata[,paste0(fov, "_", substructureMetaDataFieldName, "_", cluster_index)] <- NULL
          #store results of hull detection 
          seuratObjectMetadata[seuratObjectMetadata[, fovField] == fov,
                               paste0(fov, "_", substructureMetaDataFieldName, "_", cluster_index)] <- inhull_yes_no
        }
      }
      
      if (summarizeLocalResults) {
        seuratObjectMetadata[,paste0("Local_", substructureMetaDataFieldName)] <- FALSE
        #if any of the values in the "fov_substructureMetaDataFieldName_cluster" columns are true (i.e. the cell is within any substructure) set to TRUE. 
        cells_within_structure <- apply(seuratObjectMetadata[,grepl(paste0(fov, "_", substructureMetaDataFieldName), colnames(seuratObjectMetadata))], FUN = any, MARGIN = 1)
        
        seuratObjectMetadata[seuratObjectMetadata[,fovField] == fov & 
                               cells_within_structure,
                             paste0("Within_Local_", substructureMetaDataFieldName)] <- TRUE
        #Similarly, figure out which of the substructures the cell is in. (H)DBSCAN yields unique cluster assignment, so which.max() across the boolean cluster columns is sufficient to determine identity. 
        seuratObjectMetadata[seuratObjectMetadata[,fovField] == fov & 
                               cells_within_structure,
                             paste0("Local_", substructureMetaDataFieldName)] <- apply(seuratObjectMetadata[seuratObjectMetadata[,fovField] == fov & cells_within_structure, grepl(paste0(fov, "_", substructureMetaDataFieldName), colnames(seuratObjectMetadata))], FUN = which.max, MARGIN = 1 )
      }
    } else {
      warning(paste0("No cells of cell type(s): ", cellTypeWhiteList, " detected in FOV ", fov, "."))
    }
  }
  return(seuratObjectMetadata) 
}
