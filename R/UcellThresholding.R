
#' @title FindElbow
#' @description
#' A short description...
#' 

# Find Elbow Function 
findElbow <- function(x_sorted) {
  # x_sorted must be ascending
  n  <- length(x_sorted)
  x1 <- 1
  y1 <- x_sorted[1]
  x2 <- n
  y2 <- x_sorted[n]
  
  # Straight line from (x1, y1) to (x2, y2)
  m <- (y2 - y1) / (x2 - x1)
  b <- y1 - m * x1
  
  idx    <- seq_len(n)
  line_y <- m * idx + b
  dist   <- abs(x_sorted - line_y)
  
  elbow_idx <- which.max(dist)
  elbow_val <- x_sorted[elbow_idx]
  return(elbow_val)
}


# Ucell Thresholding function 
refine_activation_thresholds <- function(
    seuratObj,
    ucellCols,
    conditionCol  = "Antigen",
    noStimLev     = "NoStim",
    stimLev       = "SEB",
    plotElbow     = FALSE,
    scaleScores   = FALSE,
    noStimQuantile = NA  #0.75 for the top 25% of NoStim-Refined
) {
  
  # Data frame to store thresholds
  # We'll also store the NoStim quantile if used
  resultsDF <- data.frame(
    UCellCol         = ucellCols,
    ROC_Threshold    = NA_real_,
    Elbow_Threshold  = NA_real_,
    NoStim_Quantile  = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(ucellCols)) {
    scoreCol <- ucellCols[i]
    
    
    # Extract raw scores
    
    raw_scores <- seuratObj@meta.data[[scoreCol]]
    
    # Optionally scale
    if (scaleScores) {
      scaledColName <- paste0(scoreCol, "_Scaled")
      seuratObj@meta.data[[scaledColName]] <- as.numeric(scale(raw_scores))
      scores <- seuratObj@meta.data[[scaledColName]]
    } else {
      scores <- raw_scores
    }
    
    
    # Binary classification via ROC
    conditions <- factor(
      seuratObj@meta.data[[conditionCol]],
      levels = c(noStimLev, stimLev)  
    )
    
    roc_obj <- pROC::roc(
      response  = conditions,
      predictor = scores,
      levels    = c(noStimLev, stimLev)
    )
    best_threshold <- pROC::coords(
      roc     = roc_obj,
      x       = "best",
      best.method = "youden"
    )$threshold
    resultsDF$ROC_Threshold[i] <- best_threshold
    
    # Step 1 classification
    binaryColName <- paste0(scoreCol, "_Binary")
    seuratObj@meta.data[[binaryColName]] <- (scores >= best_threshold)
    
    ## 3) Elbow threshold on Stim+ Step1-positives
    idx_stim_act <- which(
      conditions == stimLev & 
        seuratObj@meta.data[[binaryColName]] == TRUE
    )
    
    if (length(idx_stim_act) < 5) {
      warning(sprintf(
        "Not enough %s-activated cells (%d) to find elbow in %s",
        stimLev, length(idx_stim_act), scoreCol
      ))
      next
    }
    
    sub_scores <- sort(scores[idx_stim_act], decreasing = FALSE)
    elbow_val  <- findElbow(sub_scores)
    resultsDF$Elbow_Threshold[i] <- elbow_val
    
    # Create a preliminary "Refined" vector
    idx_all_act <- which(seuratObj@meta.data[[binaryColName]] == TRUE)
    
    refinedColName <- paste0(scoreCol, "_Refined")
    refinedVec     <- rep(FALSE, nrow(seuratObj@meta.data))
    
    # Mark as TRUE those that were >= ROC threshold AND exceed the elbow
    refinedVec[idx_all_act] <- (scores[idx_all_act] >= elbow_val)
    
    ## 4) NoStim quantile threshold AFTER refining
    
    if (!is.na(noStimQuantile) && noStimQuantile > 0 && noStimQuantile < 1) {
      
      # Among the NoStim cells, look at those that are currently "Refined == TRUE"
      noStim_idx_refined <- which(conditions == noStimLev & refinedVec == TRUE)
      
      if (length(noStim_idx_refined) >= 5) {
        # Get their scores
        noStim_refined_scores <- scores[noStim_idx_refined]
        
        # e.g 0.75 => 75th percentile for NoStim-Refined
        noStim_thr <- as.numeric(quantile(noStim_refined_scores, probs = noStimQuantile))
        resultsDF$NoStim_Quantile[i] <- noStim_thr
        
      } else {
        warning(sprintf(
          "Not enough NoStim-refined cells (%d) to apply quantile in %s",
          length(noStim_idx_refined), scoreCol
        ))
      }
    }
    
    
    # Elbow Plot
    if (plotElbow && length(idx_stim_act) >= 5) {
      plot(
        x    = seq_along(sub_scores), 
        y    = sub_scores,
        type = "b",
        main = paste("Elbow:", scoreCol, "|", stimLev, "+ Step1 subset"),
        xlab = "Index (sorted subset)",
        ylab = paste0(scoreCol, if (scaleScores) " (scaled)" else " (raw)"),
        pch  = 16
      )
      # Add reference line (from first to last point)
      abline(
        a   = sub_scores[1], 
        b   = (tail(sub_scores,1) - sub_scores[1]) / (length(sub_scores) - 1),
        col = "red", lty = 2
      )
      # Mark the elbow
      elbow_idx <- which.min(abs(sub_scores - elbow_val))
      points(elbow_idx, elbow_val, col = "blue", pch = 19, cex = 1.5)
      text(elbow_idx, elbow_val,
           labels = formatC(elbow_val, digits = 3, format = "f"),
           pos    = 4, col = "blue")
    }
  }
  
  return(resultsDF)
  
}
