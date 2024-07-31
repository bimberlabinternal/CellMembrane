#' @include Utils.R
#' @import Seurat
#' @import ggplot2



#' @title Triage ADTs and Classify Cells
#' 
#' @description This function will classify cells within a Seurat object based on characteristics relating to the bimodality of an ADT.
#' 
#' @param seuratObj A Seurat object.
#' @param adtToInspect The name of the ADT to inspect.
#' @param assay The assay to use.
#' @param layer The layer to use.
#' @param libraryIdentificationColumn The column in the metadata that identifies the atomic unit of sequencing data (e.g. cDNA library). 
#' @param minimumCounts A lower bound on the number of reads in the ADT assay to attempt classification.
#' @param minimumCells A lower bound on the number of cells recovered from a library to attempt classification.
#' @param plots If TRUE, plots will be generated.
#' @param whitelist A vector of library IDs to run the classification on. If NULL, all libraries will be run.
#' @export

triageADTsAndClassifyCells <- function(seuratObj, 
                                       adtToInspect = 'CD4', 
                                       assay = "ADT", 
                                       layer = 'data', 
                                       libraryIdentificationColumn = 'cDNA_ID', 
                                       minimumCounts = 200, 
                                       minimumCells = 20, 
                                       plots = T, 
                                       whitelist = NULL
){
  if (!is.logical(plots)) {
    stop("The 'plots' argument must be TRUE or FALSE.")
  }
  
  #TODO: If this is going to return a seurat object with calls, it should store the default ident and reset it after running. 
  
  #harvest data from seurat object, and fix Idents to be equal to the libraryIdentificationColumn
  adtMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)
  Seurat::Idents(seuratObj) <- libraryIdentificationColumn
  
  ids <- unique(seuratObj@meta.data[, libraryIdentificationColumn])
  
  #optionally, only select a few libraries to run on
  if (!is.null(whitelist)) {
    if (!all(whitelist %in% ids)) {
      warning(paste("Whitelist values:", paste0(unique(whitelist[!(whitelist %in% ids)]), collapse = ", "), "were not found in the", libraryIdentificationColumn, "column of the Seurat object's metadata."))
    }
    #subset library IDs
    ids <- ids[ids %in% whitelist]
  }
  #initalize vectors to hold statistics to be used for classification
  decay_vector <- c()
  diff_vector <- c()
  peakheight_ratio_vector <- c()
  scaled_mode_antimode_vector <- c()
  peak1_location_vector <- c()
  peak2_location_vector <- c()
  antimode_location_vector <- c()
  peak1_density_vector <- c()
  peak2_density_vector <- c()
  antimode_density_vector <- c()
  
  for (cid in ids) {
    cells_in_library <- Seurat::WhichCells(seuratObj, idents = cid)
    # check library size and skip if too small
    if (length(cells_in_library) <= minimumCells) {
      calls <- c(calls, "Low Counts")
      decay_vector <- c(decay_vector, NA)
      diff_vector <- c(diff_vector, NA)
      peakheight_ratio_vector <- c(peakheight_ratio_vector, NA)
      scaled_mode_antimode_vector <- c(scaled_mode_antimode_vector, NA)
      peak1_location_vector <- c(peak1_location_vector, NA)
      peak2_location_vector <- c(peak2_location_vector, NA)
      antimode_location_vector <- c(antimode_location_vector, NA)
      peak1_density_vector <- c(peak1_density_vector, NA)
      peak2_density_vector <- c(peak2_density_vector, NA)
      antimode_density_vector <- c(antimode_density_vector, NA)
      next
    }
    
    library_adt_submatrix <- adtMatrix[,cells_in_library]
    #debugging
    if (plots) {
      hist(library_adt_submatrix[adtToInspect,], breaks = 100, main = cid)
    }
    
    #try to find bimodality using locmodes
    multimodeResult <- multimode::locmodes(library_adt_submatrix[adtToInspect,], mod0 = 2, lowsup=-5, uppsup = 5, n = 10000, display = plots)
    #if mutlimode::locmodes succeeds (i.e. two modes and an antimode are not NA values), calculate heuristics
    if (all(!is.na(multimodeResult$fvalue))) {
      #TODO: sorting isn't really appropriate here, but maybe when Greg can implement his previous code that loops through the mod values, it will be necessary again. 
      sorted <- multimodeResult$fvalue |> sort(decreasing = T)
      peak1 <- sorted[1]
      peak2 <- sorted[2]
      antimode <- sorted[3]
      peakheight_ratio <- peak1/peak2
      scaled_mode_antimode <- abs((antimode - peak2)/(antimode - minval))
    } else {
      #if locmodes failed, set heuristics to NA
      peak1 <- NA
      peak2 <- NA
      antimode <- NA
      peakheight_ratio <- NA
      scaled_mode_antimode <- NA
      call = "Indeterminate"
    }
    
    #apply bimodality heuristics
    if (!is.na(peakheight_ratio) & !is.na(scaled_mode_antimode)) {
      df <- data.frame(locations = multimodeResult$locations, density = multimodeResult$fvalue) |> 
        dplyr::arrange(-density)
      
      #peak-distance divded by the standard deviation of the ADT
      diff <- ((df$locations[1] - df$locations[2])^2)/sd(library_adt_submatrix[adtToInspect,])
      if (diff < peakdist_sd_ratio) {
        call <- "Fail: Diff Threshold"
      } else if (scaled_mode_antimode < scaled_mode_antimode_threshold) {
        call <- "Fail: Scaled Mode Antimode Threshold"
      } else if (peakheight_ratio > peakheight_ratio_threshold) {
        call <- "Fail: Peak Height Ratio Threshold"
      } else {
        call <- "Pass"
      }
      
      #independent of locmodes, calculate autocorrelation and fit an exponential to model the decay rate of the autocorrelation wrt lag
      p_acf <- pacf(quantile(library_adt_submatrix, probs = seq(0.01,1,0.001)), plot = plots)
      data <- data.frame(p_acf = p_acf$acf, n = p_acf$lag)
      
      #this will warn that we're taking logs of negative numbers, but our primary interest is in the first few lags. After that, the autocorrelation will have smaller magnitudes and not matter nearly as much. 
      #TODO: consider dropping data after the first lag where p_acf = 0
      model <- suppressWarnings(lm(log(p_acf) ~ 0 + n, data = data))
      decay = as.double(coef(model))
    }
    #store bimodality statistics
    decay_vector <- c(decay_vector, decay)
    diff_vector <- c(diff_vector, diff)
    peakheight_ratio_vector <- c(peakheight_ratio_vector, peakheight_ratio)
    scaled_mode_antimode_vector <- c(scaled_mode_antimode_vector, scaled_mode_antimode)
    peak1_density_vector <- c(peak1_density_vector, multimodeResult$fvalue[1])
    peak2_density_vector <- c(peak2_density_vector, multimodeResult$fvalue[3])
    antimode_density_vector <- c(antimode_density_vector, multimodeResult$fvalue[2])
    peak1_location_vector <- c(peak1_location_vector, multimodeResult$locations[1])
    peak2_location_vector <- c(peak2_location_vector, multimodeResult$locations[3])
    antimode_location_vector <- c(antimode_location_vector, multimodeResult$locations[2])
    calls <- c(calls, call)
    
  }
  
  return(data.frame(cDNA_ID = ids, 
                    call = calls, 
                    decay = decay_vector, 
                    diff = diff_vector, 
                    peakheight_ratio = peakheight_ratio_vector, 
                    scaled_mode_antimode = scaled_mode_antimode_vector, 
                    peak1_density = peak1_density_vector, 
                    peak2_density = peak2_density_vector, 
                    antimode_density = antimode_density_vector, 
                    peak1_location = peak1_location_vector, 
                    peak2_location = peak2_location_vector, 
                    antimode_location = antimode_location_vector))
  #TODO: use these features to classify cells, potentially returning a seurat object with classified cells. 
}