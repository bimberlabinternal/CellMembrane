#' @include Utils.R
#' @import Seurat
#' @import ggplot2

utils::globalVariables(	
  names = c('cDNA_ID', 'ADT', 'reason', 'antimode_location'),	
  package = 'CellMembrane',	
  add = TRUE	
)

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
#' @param xdist_ratio_threshold The threshold for the peak1_antimode_distance/peak2_antimode_distance ratio.
#' @param allowMissingADTs If TRUE, the function will attempt to remove ADTs from the whitelist that are not present in the Seurat object. If FALSE, the function will throw an error if any ADTs in the whitelist are not present in the Seurat object.
#' @param troughheight_ratio_threshold The threshold for the peak2 - antimode to peak2 ratio.
#' @param verboseplots if TRUE, intermediate plots will be generated
#' @param antimode_loc_min The minimum antimode location.
#' @export

triageADTsAndClassifyCells <- function(seuratObj,
                                       whitelist = NULL,
                                       adtwhitelist = NULL, 
                                       assay = "ADT", 
                                       layer = 'data', 
                                       libraryIdentificationColumn = 'cDNA_ID', 
                                       minimumCounts = 200, 
                                       minimumCells = 200, 
                                       plots = T,
                                       verboseplots = F,
                                       peakdist_sd_ratio = 1,
                                       peakheight_ratio_threshold = 50,
                                       troughheight_ratio_threshold = 0.3,
                                       xdist_ratio_threshold = 6,
                                       allowMissingADTs = TRUE,
                                       antimode_loc_min = 1){
  
  #TODO: check the whitelist (rename to cellmask or something? TBD)	
  
  #check rules for allowing missing ADTs	
  if (!is.logical(allowMissingADTs)) {	
    stop("The 'allowMissingADTs' argument must be TRUE or FALSE.")	
  }	
  #check ADT whitelist, and optionally allow code to run if ADTs are missing	
  if (!is.null(adtwhitelist) && !all(adtwhitelist %in% rownames(Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)))) {	
    if (!any(adtwhitelist %in% rownames(Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)))) {	
      stop(paste0("None of the ADTs in the whitelist were found in the Seurat object."))	
    } else if (!allowMissingADTs) {	
      stop(paste0("Some of the ADTs (", paste0(adtwhitelist[!(adtwhitelist %in% rownames(Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)))], collapse = ', ') , ") in the whitelist were not found in the Seurat object."))	
    } else {	
      warning(paste("ADT whitelist values:", paste0(unique(adtwhitelist[!(adtwhitelist %in% rownames(Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)))]), collapse = ", "), "were not found in the available ADTs of the Seurat object. These will be removed."))	
      adtwhitelist <- adtwhitelist[(adtwhitelist %in% rownames(Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)))]	
    }	
  }	
  #check Seurat arguments (assays, layers) 	
  if (!assay %in% Seurat::Assays(seuratObj)) {	
    stop(paste0("The assay: ", assay, " is not present in the Seurat object."))	
  }	
  if (!layer %in% SeuratObject::Layers(seuratObj)) {	
    stop(paste0("The layer: ", layer, " is not present in the Seurat object."))	
  }	
  if (!libraryIdentificationColumn %in% colnames(seuratObj@meta.data)) { 	
    stop(paste0("The libraryIdentificationColumn: ", libraryIdentificationColumn, " is not present in the Seurat object's metadata."))	
  }	
  #check validity of numeric arguments	
  if (!all(is.numeric(c(minimumCounts, minimumCells, peakdist_sd_ratio, peakheight_ratio_threshold, xdist_ratio_threshold)))) {	
    failing_arguments <- c("minimumCounts", "minimumCells", "peakdist_sd_ratio", "peakheight_ratio_threshold", "xdist_ratio_threshold")[!sapply(list(minimumCounts, minimumCells, peakdist_sd_ratio, peakheight_ratio_threshold, xdist_ratio_threshold), FUN = is.numeric)]	
    failing_values <- c(minimumCounts, minimumCells, peakdist_sd_ratio, peakheight_ratio_threshold, xdist_ratio_threshold)[!sapply(list(minimumCounts, minimumCells, peakdist_sd_ratio, peakheight_ratio_threshold, xdist_ratio_threshold), FUN = is.numeric)]	
    stop(paste0("The ", paste0(failing_arguments, collapse = ', '), " argument(s) must be numeric. Current value(s) is/are : " , paste0(failing_arguments ,":", failing_values, collapse = ', ')))	
  }
  
  if (!is.logical(plots)) {
    stop("The 'plots' argument must be TRUE or FALSE.")
  }
  
  df <- .CalculateStatistics(seuratObj = seuratObj, assay = assay, layer=layer,
                             libraryIdentificationColumn=libraryIdentificationColumn,
                             minimumCounts=minimumCounts, minimumCells=minimumCells,
                             plots=plots, verboseplots = verboseplots, adtwhitelist=adtwhitelist,
                             whitelist=whitelist, peakdist_sd_ratio=peakdist_sd_ratio,
                             xdist_ratio_threshold = xdist_ratio_threshold,
                             troughheight_ratio_threshold = troughheight_ratio_threshold,
                             antimode_loc_min = antimode_loc_min)
  df <- .TriageADTs(df = df, minimumCells=minimumCells, peakdist_sd_ratio=peakdist_sd_ratio,
                    peakheight_ratio_threshold=peakheight_ratio_threshold,
                    xdist_ratio_threshold=xdist_ratio_threshold,
                    troughheight_ratio_threshold = troughheight_ratio_threshold,
                    plots=plots, antimode_loc_min = antimode_loc_min)
  seuratObj <- .ClassifyCells(seuratObj=seuratObj, df=df, libraryIdentificationColumn=libraryIdentificationColumn, assay=assay, layer=layer, adtwhitelist=adtwhitelist)
  
  return(seuratObj)
}

#' @title Calculate Statistics for ADT distributions
#' 
#' @description This function will calculate metrics and produce a dataframe that will be used for classification by .TriageADTs
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
#' @param verboseplots if TRUE, intermediate plots will be generated
#' @param peakdist_sd_ratio An upper bound of acceptable peak distance divided by the standard deviation of the ADT.
#' @param peakheight_ratio_threshold The threshold for the peak height ratio.
#' @param troughheight_ratio_threshold The threshold for the peak2 - antimode to peak2 ratio.
#' @param xdist_ratio_threshold The threshold for the peak1_antimode_distance/peak2_antimode_distance ratio.
#' @param antimode_loc_min The minimum antimode location.

.CalculateStatistics <- function(seuratObj,
                                 whitelist = NULL,
                                 adtwhitelist = NULL, 
                                 assay = "ADT", 
                                 layer = 'data', 
                                 libraryIdentificationColumn = 'cDNA_ID', 
                                 minimumCounts = 200, 
                                 minimumCells = 20, 
                                 plots = T,
                                 verboseplots = F,
                                 peakdist_sd_ratio = 1,
                                 peakheight_ratio_threshold = 50,
                                 troughheight_ratio_threshold = 0.3,
                                 xdist_ratio_threshold = 6,
                                 antimode_loc_min = 1){
  
  #harvest data from seurat object, and fix Idents to be equal to the libraryIdentificationColumn
  adtMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)
  countsMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, layer = "counts")
  
  ids <- unique(seuratObj@meta.data[, libraryIdentificationColumn])
  
  #optionally, only select a few libraries to run on
  if (!is.null(whitelist)) {
    if (!all(whitelist %in% ids)) {
      warning(paste("cDNA whitelist values:", paste0(unique(whitelist[!(whitelist %in% ids)]), collapse = ", "), "were not found in the", libraryIdentificationColumn, "column of the Seurat object's metadata."))
    }
    #subset library IDs
    ids <- ids[ids %in% whitelist]
  }
  
  if (!is.null(adtwhitelist)){
    if (!all(adtwhitelist %in% rownames(adtMatrix))) {
      warning(paste("ADT whitelist values:", paste0(unique(adtwhitelist[!(adtwhitelist %in% rownames(adtMatrix))]), collapse = ", "), "were not found in the available ADTs of the Seurat object. These will be ignored"))
    }
    adts <- adtwhitelist[(adtwhitelist %in% rownames(adtMatrix))]
  } else {
    adts <- rownames(adtMatrix)
  }
  
  #initalize vectors to hold statistics to be used for classification
  cells_vector <- c()
  diff_vector <- c()
  peakheight_ratio_vector <- c()
  troughheight_ratio_vector <- c()
  peak1_location_vector <- c()
  peak2_location_vector <- c()
  antimode_location_vector <- c()
  peak1_density_vector <- c()
  peak2_density_vector <- c()
  antimode_density_vector <- c()
  cid_vector <- c()
  adt_vector <- c()
  geom_mean_vector <- c()
  stdev_vector <- c()
  med_vector <- c()
  mod0vec <- c()
  bandwidth_vec <- c()
  
  #TODO: trim these if unused: 
  under_threshold_density_vector <- c()
  over_threshold_density_vector <- c()
  under_threshold_percentage_vector <- c()
  over_threshold_percentage_vector <- c()
  
  
  for (cid in ids) {
    cells_in_library <- colnames(seuratObj)[seuratObj@meta.data[[libraryIdentificationColumn]] == cid]
    # check library size and skip if too small
    if (length(cells_in_library) <= minimumCells) {
      cells_vector <- c(cells_vector, rep(length(cells_in_library), length(adts)))
      diff_vector <- c(diff_vector, rep(NA, length(adts)))
      peakheight_ratio_vector <- c(peakheight_ratio_vector, rep(NA, length(adts)))
      troughheight_ratio_vector <- c(troughheight_ratio_vector, rep(NA, length(adts)))
      peak1_location_vector <- c(peak1_location_vector, rep(NA, length(adts)))
      peak2_location_vector <- c(peak2_location_vector, rep(NA, length(adts)))
      antimode_location_vector <- c(antimode_location_vector, rep(NA, length(adts)))
      peak1_density_vector <- c(peak1_density_vector, rep(NA, length(adts)))
      peak2_density_vector <- c(peak2_density_vector, rep(NA, length(adts)))
      antimode_density_vector <- c(antimode_density_vector, rep(NA, length(adts)))
      adt_vector <- c(adt_vector, adts)
      cid_vector <- c(cid_vector, rep(cid, length(adts)))
      geom_mean_vector <- c(geom_mean_vector, rep(NA, length(adts)))
      stdev_vector <- c(stdev_vector, rep(NA, length(adts)))
      med_vector <- c(med_vector, rep(NA, length(adts)))
      mod0vec <- c(mod0vec, rep(NA, length(adts)))
      bandwidth_vec <- c(bandwidth_vec, rep(NA, length(adts)))
      
      #TODO: trim these if unused
      under_threshold_density_vector <- c(under_threshold_density_vector, rep(NA, length(adts)))
      over_threshold_density_vector <- c(over_threshold_density_vector, rep(NA, length(adts)))
      under_threshold_percentage_vector <- c(under_threshold_percentage_vector, rep(NA, length(adts)))
      over_threshold_percentage_vector <- c(over_threshold_percentage_vector, rep(NA, length(adts)))
      
      next
    }
    # Iterate over all ADTs in cDNA_Id
    for (adt in adts) {
      library_adt_vector <- adtMatrix[adt,cells_in_library]
      library_counts_vector <- countsMatrix[adt, cells_in_library]
      #debugging
      for (mod0val in seq(2,8,1)) {
        #try to find bimodality using locmodes
        multimodeResult <- multimode::locmodes(library_adt_vector, mod0 = mod0val, lowsup=0.2, uppsup = 10, n = 10000, display = verboseplots)
        #if mutlimode::locmodes succeeds (i.e. two modes and an antimode are not NA values), calculate heuristics
        if (sum(!is.na(multimodeResult$fvalue))>=3) {
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
          troughheight_ratio <- (peak2-antimode)/peak2
          diff <- ((peak1_location - peak2_location)^2)/(sd(library_adt_vector)^2)
          bandwidth <- multimodeResult$cbw$bw
          # if a valid solution is found to antimode identification, stop iterating mod0val and accept modes
          if ((diff > peakdist_sd_ratio) & (peakheight_ratio < peakheight_ratio_threshold) & 
              troughheight_ratio > troughheight_ratio_threshold & antimode_location > antimode_loc_min)  {
            break
          }
        } else {
          #if locmodes failed, set measurement values to NA
          peak1 <- NA
          peak2 <- NA
          peak1_location <- NA
          peak2_location <- NA
          antimode <- NA
          antimode_location <- NA
          peakheight_ratio <- NA
          troughheight_ratio <- NA
          diff <- NA
          bandwidth <- NA
        }
      }
      #compute library statistics (independent of bimodality results)
      geom_mean <- exp(mean(log(library_counts_vector + 1)))
      stdev <- sd(library_adt_vector)
      med <- stats::median(library_adt_vector)
      
      #compute bimodality-dependent statistics
      #TODO: remove whatever doesn't get used for thresholding from here
      if (!is.na(antimode_location)) {
        #number of cells above and below the accepted threshold
        under_threshold_density <- sum(library_adt_vector < antimode_location)
        over_threshold_density <- sum(library_adt_vector > antimode_location)
        #as percentages
        under_threshold_percentage <- under_threshold_density/length(library_adt_vector)
        over_threshold_percentage <- over_threshold_density/length(library_adt_vector)
      } else {
        under_threshold_density <- NA
        over_threshold_density <- NA
        under_threshold_percentage <- NA
        over_threshold_percentage <- NA
      }
      
      #store bimodality statistics
      cells_vector <- c(cells_vector, length(cells_in_library))
      diff_vector <- c(diff_vector, diff)
      peakheight_ratio_vector <- c(peakheight_ratio_vector, peakheight_ratio)
      troughheight_ratio_vector <- c(troughheight_ratio_vector, troughheight_ratio)
      peak1_density_vector <- c(peak1_density_vector, peak1)
      peak2_density_vector <- c(peak2_density_vector, peak2)
      antimode_density_vector <- c(antimode_density_vector, antimode)
      peak1_location_vector <- c(peak1_location_vector, peak1_location)
      peak2_location_vector <- c(peak2_location_vector, peak2_location)
      antimode_location_vector <- c(antimode_location_vector, antimode_location)
      cid_vector <- c(cid_vector, cid)
      adt_vector <- c(adt_vector, adt)
      geom_mean_vector <- c(geom_mean_vector, geom_mean)
      stdev_vector <- c(stdev_vector, stdev)
      med_vector <- c(med_vector, med)
      mod0vec <- c(mod0vec, mod0val)
      bandwidth_vec <- c(bandwidth_vec, bandwidth)
      
      #TODO: trim these if unused:
      under_threshold_density_vector <- c(under_threshold_density_vector, under_threshold_density)
      over_threshold_density_vector <- c(over_threshold_density_vector, over_threshold_density)
      under_threshold_percentage_vector <- c(under_threshold_percentage_vector, under_threshold_percentage)
      over_threshold_percentage_vector <- c(over_threshold_percentage_vector, over_threshold_percentage)
      
      if (plots) {
        # graphics::hist(library_adt_vector, breaks = 100, main = paste(cid, ", ", adt))
        dx <- density(library_adt_vector, bw =multimodeResult$cbw$bw)
        df_dx <- data.frame(x = dx$x, y = dx$y)
        max_y1 <- max(df_dx$y)
        max_y2 <- max(hist(library_adt_vector, breaks = 100, plot = F)$counts)
        xdist_ratio <- round((abs(peak1_location - antimode_location))/(abs(peak2_location - antimode_location)),2)
        
        xlabel <- ifelse(length(cells_in_library) > minimumCells,
                         "<span style = 'color:cornflowerblue;'>Number of Cells: </span>",
                         "<span style = 'color:red;'>Number of Cells: </span>")
        
        titlepeak <- ifelse(peakheight_ratio < peakheight_ratio_threshold,
                            "<span style = 'color:cornflowerblue;'>Height: </span>",
                            "<span style = 'color:red;'>Height: </span>")
        
        titledist <- ifelse(xdist_ratio < xdist_ratio_threshold,
                            "<span style = 'color:cornflowerblue;'>Distance: </span>",
                            "<span style = 'color:red;'>Distance: </span>")
        
        titlediff <- ifelse(diff > peakdist_sd_ratio,
                            "<span style = 'color:cornflowerblue;'>Variance: </span>",
                            "<span style = 'color:red;'>Variance: </span>")
        
        titletrough <- ifelse(troughheight_ratio > troughheight_ratio_threshold,
                              "<span style = 'color:cornflowerblue;'>Valley: </span>",
                              "<span style = 'color:red;'>Valley: </span>")
        
        title_minloc <- ifelse(antimode_location > antimode_loc_min,
                               "<span style = 'color:cornflowerblue;'>Threshold: </span>",
                               "<span style = 'color:red;'>Threshold: </span>"
                               )
        
        plt <- ggplot(as.data.frame(library_adt_vector), aes(x = library_adt_vector)) + 
          geom_histogram(bins = 100) + geom_line(data = df_dx, aes(x = dx$x, y=(max_y2/max_y1)*dx$y, color = "Density")) + 
          geom_vline(xintercept = c(peak1_location, peak2_location), color = "skyblue") +
          geom_vline(xintercept = c(antimode_location), linetype = "dashed", color = "dodgerblue") +
          # labs(x = paste0("ADT Signal (", length(cells_in_library)," cells)"), y = "Density",
          labs(x = paste0("ADT Signal (", xlabel, length(cells_in_library),")"), y = "Density",
               title = paste0(cid, "-", adt), subtitle = paste0(titlediff, round(diff, 2),", ", titlepeak, round(peakheight_ratio,2),
                                                                ", ", titledist, xdist_ratio, ", ", titletrough, round(troughheight_ratio, 3),
                                                                ", ", title_minloc, round(antimode_location, 2))) + 
          egg::theme_article() + 
          
          theme(
            axis.title.x = ggtext::element_markdown(),
            plot.subtitle = ggtext::element_markdown())
        print(plt)
      }
    }
  }
  
  return(data.frame(cDNA_ID = cid_vector,
                    ADT = adt_vector,
                    cells = cells_vector,
                    diff = diff_vector, 
                    peakheight_ratio = peakheight_ratio_vector,
                    troughheight_ratio = troughheight_ratio_vector,
                    peak1_density = peak1_density_vector, 
                    peak2_density = peak2_density_vector, 
                    antimode_density = antimode_density_vector, 
                    peak1_location = peak1_location_vector, 
                    peak2_location = peak2_location_vector, 
                    antimode_location = antimode_location_vector,
                    geom_mean = geom_mean_vector,
                    stdev = stdev_vector,
                    med = med_vector,
                    mod0val = mod0vec,
                    bandwidth = bandwidth_vec,
                    under_density = under_threshold_density_vector,
                    over_density = over_threshold_density_vector,
                    under_percentage = under_threshold_percentage_vector,
                    over_percentage = over_threshold_percentage_vector
                    ))
}

#' @title Triage ADTs for classification
#' 
#' @description This function will take in a dataframe and make calls on bimodality
#' 
#' @param df A dataframe resulting from running .CalculateStatistics on a Seurat object
#' @param minimumCells A lower bound on the number of cells recovered from a library to attempt classification.
#' @param peakdist_sd_ratio An upper bound of acceptable peak distance divided by the standard deviation of the ADT.
#' @param peakheight_ratio_threshold The threshold for the peak height ratio.
#' @param xdist_ratio_threshold The threshold for the peak1_antimode_distance/peak2_antimode_distance ratio.
#' @param troughheight_ratio_threshold The threshold for the peak2 - antimode to peak2 ratio.
#' @param plots If TRUE, plots will be generated.
#' @param antimode_loc_min The minimum antimode location.

.TriageADTs <- function(df, minimumCells, peakdist_sd_ratio, peakheight_ratio_threshold,
                        xdist_ratio_threshold, troughheight_ratio_threshold, plots, antimode_loc_min){
  print(paste0("minimumCells >  ", minimumCells))
  print(paste0("peakdist_sd_ratio > ", peakdist_sd_ratio))
  print(paste0("peakheight_ratio < ", peakheight_ratio_threshold))
  print(paste0("xdist_ratio < ", xdist_ratio_threshold))
  print(paste0("troughheight_ratio > ", troughheight_ratio_threshold))
  print(paste0("antimode_loc > ", antimode_loc_min))
  
  
  df$xdist_ratio <- (abs(df$peak1_location - df$antimode_location))/(abs(df$peak2_location - df$antimode_location))
  df$reason <- dplyr::case_when(
    df$cells < minimumCells ~ "Too few cells",
    is.na(df$antimode_location) ~ "Probable ADT failure",
    df$antimode_location < antimode_loc_min ~ "Low Antimode Location",
    df$peakheight_ratio > peakheight_ratio_threshold ~ "Peak Height Ratio",
    df$troughheight_ratio < troughheight_ratio_threshold ~ "Trough Height Ratio",
    df$diff < peakdist_sd_ratio ~ "Diff Threshold",
    df$xdist_ratio > xdist_ratio_threshold ~ "Distance Ratio",
    .default = "Pass"
  )
  
  df$call <- dplyr::case_when(
    df$reason == "Pass" ~ "Pass",
    .default = "Fail"
  )
  
  if (plots) {
    print(ggplot2::ggplot(df, aes(x = peak1_density, y = peak1_location, color = reason, shape = call)) +
            geom_point() + egg::theme_article() + facet_wrap(~ADT))
    
    print(ggplot2::ggplot(df, aes(x = peak2_density, y = peak2_location, color = reason, shape = call)) + 
            geom_point() + egg::theme_article() + facet_wrap(~ADT))
    
    print(ggplot2::ggplot(df, aes(x = antimode_density, y = antimode_location, color = reason, shape = call)) + 
            geom_point() + egg::theme_article() + facet_wrap(~ADT))
    
    print(ggplot2::ggplot(df, aes(y = peakheight_ratio, x = troughheight_ratio, color = reason, shape = call)) +  
            geom_point() + egg::theme_article() + scale_y_log10() + facet_wrap(~ADT))
    
    print(ggplot2::ggplot(df, aes(y = peakheight_ratio, x = diff, color = reason, shape = call))  + 
            geom_point() + egg::theme_article() + scale_y_log10()  + scale_x_log10() + facet_wrap(~ADT))
    
    print(ggplot2::ggplot(df, aes(x = xdist_ratio, y = peakheight_ratio, color = reason, shape = call))  + 
            geom_point() + egg::theme_article() + facet_wrap(~ADT) + scale_y_log10()  + scale_x_log10()) 
    
    print(ggplot2::ggplot(df, aes(x = troughheight_ratio, y = diff, color = reason, shape = call)) +
            geom_point() + egg::theme_article() + facet_wrap(~ADT)  + scale_y_log10())
    
    print(ggplot2::ggplot(df, aes(x = troughheight_ratio, y = xdist_ratio, color = reason, shape = call)) +
            geom_point() + egg::theme_article() + facet_wrap(~ADT)  + scale_y_log10())
    
    print(ggplot2::ggplot(df, aes(x = diff, y = xdist_ratio, color = reason, shape = call)) +
            geom_point() + egg::theme_article() + facet_wrap(~ADT)  + scale_y_log10()  + scale_x_log10())
  }
  return(df)
}

#' @title Classify Cells
#' 
#' @description This function will classify cells within a Seurat object based on output from .TriageADTs.
#' 
#' @param seuratObj A Seurat object.
#' @param df A dataframe resulting from running .TriageADTs
#' @param adtwhitelist An optional vector of ADTs to be triaged and classified.
#' @param assay The assay to use.
#' @param layer The layer to use.
#' @param libraryIdentificationColumn The column in the metadata that identifies the atomic unit of sequencing data (e.g. cDNA library). 

.ClassifyCells <- function(seuratObj, df, libraryIdentificationColumn, assay, layer, adtwhitelist) {
  adtMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, layer = layer)
  if (!is.null(adtwhitelist)){
    if (!all(adtwhitelist %in% rownames(adtMatrix))) {
      warning(paste("ADT whitelist values:", paste0(unique(adtwhitelist[!(adtwhitelist %in% rownames(adtMatrix))]), collapse = ", "), "were not found in the available ADTs of the Seurat object. These will be ignored"))
    }
    adts <- adtwhitelist[(adtwhitelist %in% rownames(adtMatrix))]
  }
  else {
    adts <- rownames(adtMatrix)
  }
  df_wider <- tidyr::pivot_wider(df, id_cols = cDNA_ID, names_from = ADT,
                                 values_from = c(reason, call, antimode_location))
  meta <- seuratObj@meta.data
  meta$barcode <- rownames(meta)
  meta <- merge(meta, df_wider, by.y = "cDNA_ID", by.x = libraryIdentificationColumn)
  rownames(meta) <- meta$barcode
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = meta, col.name = colnames(df_wider |> dplyr::select(-cDNA_ID)))
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = data.frame(t(as.matrix(adtMatrix))), col.name = adts)
  for (adt in adts) {
    seuratObj@meta.data[,paste0(adt, "_cellcall")] <- dplyr::case_when(
      seuratObj@meta.data[,adt] >= seuratObj@meta.data[,paste0("antimode_location_", adt)] & seuratObj@meta.data[,paste0("call_", adt)] == "Pass" ~ "Positive",
      seuratObj@meta.data[,adt] <  seuratObj@meta.data[,paste0("antimode_location_", adt)] & seuratObj@meta.data[,paste0("call_", adt)] == "Pass" ~ "Negative",
      seuratObj@meta.data[,paste0("call_", adt)] == "Fail" ~ "Failed"
    )
  }
  return(seuratObj)
}

