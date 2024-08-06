#' @include Utils.R
#' @import Seurat
#' @import ggplot2



#' @title Triage ADTs and Classify Cells
#' 
#' @description This function will classify cells within a Seurat object based on characteristics relating to the bimodality of an ADT.
#' 
#' @param seuratObj A Seurat object.
#' @param whitelist A vector of library IDs to run the classification on. If NULL, all libraries will be run.
#' @param adtwhitelist An optional vector of ADTs to be triaged and classified.
#' @param assay The assay to use.
#' @param layer The layer to use.
#' @param libraryIdentificationColumn The column in the metadata that identifies the atomic unit of sequencing data (e.g. cDNA library). 
#' @param minimumCounts A lower bound on the number of reads in the ADT assay to attempt classification.
#' @param minimumCells A lower bound on the number of cells recovered from a library to attempt classification.
#' @param plots If TRUE, plots will be generated.
#' @param peakdist_sd_ratio An upper bound of acceptable peak distance divided by the standard deviation of the ADT.
#' @param peakheight_ratio_threshold The threshold for the peak height ratio.
#' @export

triageADTsAndClassifyCells <- function(seuratObj,
                                       adtwhitelist = NULL, 
                                       assay = "ADT", 
                                       layer = 'data', 
                                       libraryIdentificationColumn = 'cDNA_ID', 
                                       minimumCounts = 200, 
                                       minimumCells = 20, 
                                       plots = T, 
                                       whitelist = NULL, 
                                       peakdist_sd_ratio = 0.1,
                                       peakheight_ratio_threshold = 50
){
  if (!is.logical(plots)) {
    stop("The 'plots' argument must be TRUE or FALSE.")
  }
  
  #TODO: If this is going to return a seurat object with calls, it should store the default ident and reset it after running. 
  
  #harvest data from seurat object, and fix Idents to be equal to the libraryIdentificationColumn
  adtMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)
  
  ids <- unique(seuratObj@meta.data[, libraryIdentificationColumn])
  
  #optionally, only select a few libraries to run on
  if (!is.null(whitelist)) {
    if (!all(whitelist %in% ids)) {
      warning(paste("cDNA whitelist values:", paste0(unique(whitelist[!(whitelist %in% ids)]), collapse = ", "), "were not found in the", libraryIdentificationColumn, "column of the Seurat object's metadata."))
    }
    #subset library IDs
    ids <- ids[ids %in% whitelist]
  }
  #initalize vectors to hold statistics to be used for classification
  calls <- c()
  cells_vector <- c()
  decay_vector <- c()
  diff_vector <- c()
  peakheight_ratio_vector <- c()
  peak1_location_vector <- c()
  peak2_location_vector <- c()
  antimode_location_vector <- c()
  peak1_density_vector <- c()
  peak2_density_vector <- c()
  antimode_density_vector <- c()
  reason_vector <- c()
  cid_vector <- c()
  adt_vector <- c()
  
  for (cid in ids) {
    cells_in_library <- colnames(seuratObj)[seuratObj@meta.data[[libraryIdentificationColumn]] == cid]
    # check library size and skip if too small
    if (!is.null(adtwhitelist)){
      if (!all(adtwhitelist %in% rownames(adtMatrix))) {
        warning(paste("ADT whitelist values:", paste0(unique(adtwhitelist[!(adtwhitelist %in% rownames(adtMatrix))]), collapse = ", "), "were not found in the available ADTs of the Seurat object. These will be ignored"))
      }
      adts <- adtwhitelist[(adtwhitelist %in% rownames(adtMatrix))]
    }
    else {
      adts <- rownames(adtMatrix)
    }
    if (length(cells_in_library) <= minimumCells) {
      cells_vector <- c(cells_vector, length(cells_in_library))
      reason_vector <- c(reason_vector, rep("Low Counts", length(adts)))
      calls <- c(calls, rep("Fail", length(adts)))
      decay_vector <- c(decay_vector, rep(NA, length(adts)))
      diff_vector <- c(diff_vector, rep(NA, length(adts)))
      peakheight_ratio_vector <- c(peakheight_ratio_vector, rep(NA, length(adts)))
      peak1_location_vector <- c(peak1_location_vector, rep(NA, length(adts)))
      peak2_location_vector <- c(peak2_location_vector, rep(NA, length(adts)))
      antimode_location_vector <- c(antimode_location_vector, rep(NA, length(adts)))
      peak1_density_vector <- c(peak1_density_vector, rep(NA, length(adts)))
      peak2_density_vector <- c(peak2_density_vector, rep(NA, length(adts)))
      antimode_density_vector <- c(antimode_density_vector, rep(NA, length(adts)))
      adt_vector <- c(adt_vector, adts)
      cid_vector <- c(cid_vector, rep(cid, length(adts)))
      next
    }
    # Iterate over all ADTs in cDNA_Id
    for (adt in adts) {
      
      library_adt_vector <- adtMatrix[adt,cells_in_library]
      #debugging
      if (plots) {
        hist(library_adt_vector, breaks = 100, main = paste(cid, ", ", adt))
      }
      for (mod0val in seq(2,8,1)) {
        #try to find bimodality using locmodes
        multimodeResult <- multimode::locmodes(library_adt_vector, mod0 = mod0val, lowsup=0, uppsup = 5, n = 10000, display = plots)
        #if mutlimode::locmodes succeeds (i.e. two modes and an antimode are not NA values), calculate heuristics
        if (sum(!is.na(multimodeResult$fvalue))>=3) {
          #TODO: sorting isn't really appropriate here, but maybe when Greg can implement his previous code that loops through the mod values, it will be necessary again. 
          density_ordered_modes <- order(multimodeResult$fvalue, decreasing = T)
          sorted <- multimodeResult$fvalue[density_ordered_modes]
          peak1 <- sorted[1]
          peak2 <- sorted[2]
          peak1_location <- multimodeResult$locations[density_ordered_modes][1]
          peak2_location <- multimodeResult$locations[density_ordered_modes][2]
          midpoint <- (peak1_location + peak2_location)/2
          antimode_index <- which.min(abs(multimodeResult$locations - midpoint))
          antimode <- multimodeResult$fvalue[antimode_index]
          antimode_location <- multimodeResult$locations[antimode_index]
          peakheight_ratio <- (peak1-antimode)/(peak2-antimode)
          diff <- ((peak1_location - peak2_location)^2)/sd(library_adt_vector)
          # if a valid solution is found to antimode identification, stop iterating mod0val and accept modes
          if ((diff > peakdist_sd_ratio) & (peakheight_ratio < peakheight_ratio_threshold))  {
            break
          }
          
        } else if (sum(!is.na(multimodeResult$fvalue)) < 3 | mod0val < 10) {
          mod0val <- mod0val + 1
        }
        else {
          #if locmodes failed, set heuristics to NA
          peak1 <- NA
          peak2 <- NA
          antimode <- NA
          peakheight_ratio <- NA
          call <- "Fail"
          reason <- "Indeterminate"
        }
      }
      
      #apply bimodality heuristics
      if (!is.na(peakheight_ratio)) {
        df <- data.frame(locations = multimodeResult$locations, density = multimodeResult$fvalue) |> 
          dplyr::arrange(-density)
        
        #peak-distance divded by the standard deviation of the ADT
        diff <- ((df$locations[1] - df$locations[2])^2)/sd(library_adt_vector)
        if (diff < peakdist_sd_ratio) {
          call <- "Fail"
          reason <- "Diff Threshold"
        } else if (peakheight_ratio > peakheight_ratio_threshold) {
          call <- "Fail"
          reason <- "Peak Height Ratio Threshold"
        } else {
          call <- "Pass"
          reason <- "Pass"
        }
        
        if (length(unique(library_adt_vector)) > round(minimumCells/2, 0)) {
          #independent of locmodes, calculate autocorrelation and fit an exponential to model the decay rate of the autocorrelation wrt lag
          p_acf <- pacf(quantile(library_adt_vector, probs = seq(0.01,1,0.001)), plot = plots)
          data <- data.frame(p_acf = p_acf$acf, n = p_acf$lag)
          
          #this will warn that we're taking logs of negative numbers, but our primary interest is in the first few lags. After that, the autocorrelation will have smaller magnitudes and not matter nearly as much. 
          #TODO: consider dropping data after the first lag where p_acf = 0
          model <- suppressWarnings(lm(log(p_acf) ~ 0 + n, data = data))
          decay = as.double(coef(model))}
      } else {
        decay <- NA
        call <- "Fail"
        reason <- "Probable ADT failure"
      }
      #store bimodality statistics
      cells_vector <- c(cells_vector, length(cells_in_library))
      decay_vector <- c(decay_vector, decay)
      diff_vector <- c(diff_vector, diff)
      peakheight_ratio_vector <- c(peakheight_ratio_vector, peakheight_ratio)
      peak1_density_vector <- c(peak1_density_vector, peak1)
      peak2_density_vector <- c(peak2_density_vector, peak2)
      antimode_density_vector <- c(antimode_density_vector, antimode)
      peak1_location_vector <- c(peak1_location_vector, peak1_location)
      peak2_location_vector <- c(peak2_location_vector, peak2_location)
      antimode_location_vector <- c(antimode_location_vector, antimode_location)
      calls <- c(calls, call)
      reason_vector <- c(reason_vector, reason)
      cid_vector <- c(cid_vector, cid)
      adt_vector <- c(adt_vector, adt)
    }
  }
  return(data.frame(cDNA_ID = cid_vector,
                    ADT = adt_vector,
                    call = calls,
                    reason = reason_vector,
                    cells = cells_vector,
                    decay = decay_vector, 
                    diff = diff_vector, 
                    peakheight_ratio = peakheight_ratio_vector, 
                    peak1_density = peak1_density_vector, 
                    peak2_density = peak2_density_vector, 
                    antimode_density = antimode_density_vector, 
                    peak1_location = peak1_location_vector, 
                    peak2_location = peak2_location_vector, 
                    antimode_location = antimode_location_vector))
  #TODO: use these features to classify cells, potentially returning a seurat object with classified cells. 
}