#' @import ggplot2 Seurat dplyr
#' @importFrom stats cor 

utils::globalVariables(
  names = c('ClusterProportion', 'Proportion', 'SizeFactor', 'XY_Key', 'Y_Key'),
  package = 'CellMembrane',
  add = TRUE
)

#' @title .GenerateSizeFactor 
#'
#' @description A helper function for ConstructEnrichmentDataFrameAndDoStatistics that calculates size factors for differently sized scRNA-Seq datasets.
#' @param seuratObj The seurat object that holds the data.
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 

.GenerateSizeFactor <- function(seuratObj, normalizationField = 'cDNA_ID', sizeFactorField = 'SizeFactor', maxSizeFactor = 100){
  if (!normalizationField %in% names(seuratObj@meta.data)) {
    stop(paste0('Field missing: ', normalizationField))
  }
 if (any(is.na(seuratObj@meta.data[,normalizationField]))){
   stop(paste0("Detected NAs in ", normalizationField,". Please remove them."))
 }
  dat <- as.character(seuratObj@meta.data[,normalizationField])

  #check max size factor for data
  dataMaxSizeFactor <- max(table(dat))/min(table(dat))
  if (!maxSizeFactor >= dataMaxSizeFactor){
    stop(paste0("Greater than max allowable size factor: ", round(dataMaxSizeFactor,2)))
  }

  #Calculate the ratio between maximum size
  print("Calculating size factors")
  for (dataset in unique(dat)){
    datasetSizeFactor <- max(table(dat)) / length(dat[dat == dataset])
    seuratObj@meta.data[[sizeFactorField]][dat == dataset] <- round(datasetSizeFactor, 4)
  }

  return(seuratObj)
}

#' @title ConstructEnrichmentDataFrameAndDoStatistics
#'
#' @description Helper function for makeDotPlot that does the dataframe construction and statistics. 
#' @param seuratObj The seurat object that holds the data.
#' @param xField The x axis for the dotplots.
#' @param yField The y axis for the dotplots.
#' @param extraGroupingFields An optional vector of additional fields to use in grouping
#' @param colorField The column of metadata that is used for the colorField
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 
#' @param independentVariableTestField This and dependentVariableTestField automatically define a small statistical test to see if your size factors are correlated. Ideally, they should not be if independentVariableTestField is not the same value as normalizationField. 
#' @param dependentVariableTestField This value should be equal to sizeFactorField initially, but can be changed to interactively see other correlations in the metadata.
#' @export

ConstructEnrichmentDataFrameAndDoStatistics <- function(seuratObj,
                                                        yField = 'ClusterNames_0.2',
                                                        xField = 'Timepoint',
                                                        extraGroupingFields = NULL,
                                                        colorField = 'Population',

                                                        normalizationField = 'cDNA_ID',
                                                        sizeFactorField = 'SizeFactor',
                                                        maxSizeFactor = 100,
                                                        independentVariableTestField = colorField,
                                                        dependentVariableTestField = sizeFactorField
){
  
  # Calculate Size Factor
  seuratObj <- .GenerateSizeFactor(seuratObj, normalizationField = normalizationField, sizeFactorField = sizeFactorField, maxSizeFactor = maxSizeFactor)

  if (!colorField %in% names(seuratObj@meta.data)) {
    stop(paste0('colorField not found: ', colorField))
  }

  colorData <- seuratObj@meta.data[[colorField]]
  if (!is.factor(colorData)) {
    stop('Expected colorField to be a factor')
  }

  if (length(unique(colorData)) != 2) {
    stop('Expected colorField to be a factor with two levels only')
  }

  baseValue <- levels(colorData)[2]

  if (is.null(extraGroupingFields)) {
  	extraGroupingFields <- c()
  }
  
  # Construct DataFrame:
  rawData <- as_tibble(seuratObj@meta.data[,c(colorField, yField, xField, normalizationField, sizeFactorField, extraGroupingFields)])
  names(rawData) <- c('colorField', 'yField', 'xField', normalizationField, sizeFactorField, extraGroupingFields)
  rawData$xField <- as.character(rawData$xField)
  rawData$yField <- as.character(rawData$yField)

  # Make concatenated columns for grouping:
  rawData <- rawData %>% tidyr::unite("XY_Key", all_of(c('xField', 'yField', extraGroupingFields)), remove = FALSE)
  rawData <- rawData %>% tidyr::unite("Y_Key", all_of(c('yField', extraGroupingFields)), remove = FALSE)
  
  # Calculate the weighted total of cells in each X/Y group
  xyTotals <- rawData %>% dplyr::count(XY_Key, wt = SizeFactor, name = 'TotalPerXY')
  
  # Now calculate the weighted total, including the color field
  colorProportions <- rawData %>% dplyr::count(XY_Key, colorField, wt = SizeFactor, name = 'TotalPerGroup')

  # Merge and calculate proportion of cells:
  colorProportions <- merge(colorProportions, xyTotals, all.x = T, by = 'XY_Key')
  colorProportions$Proportion <- colorProportions$TotalPerGroup / colorProportions$TotalPerXY
  colorProportions <- colorProportions[colorProportions$colorField == baseValue,] # Note: this could leave zeros for certain X/Y pairs where there is data for one colorValue but not another
  if (nrow(colorProportions) == 0) {
  	stop(paste0('baseValue not found: ', baseValue))
  }
  
  colorProportions <- colorProportions[c('XY_Key', 'Proportion')]

  #initial cluster enrichment tibble
  clusterProportions <- rawData %>% dplyr::count(XY_Key, Y_Key, wt = SizeFactor)
  clusterProportions <- clusterProportions %>% group_by(Y_Key) %>% mutate(ClusterProportion = prop.table(n))

  # Merge cluster enrichment and category enrichment tibbles.
  finalData <- merge(colorProportions, clusterProportions, by = "XY_Key")
  metadata <- unique(rawData[all_of(c('XY_Key', 'xField', 'yField', extraGroupingFields))])
  finalData <- merge(finalData, metadata, by = "XY_Key")
  finalData$xField <- naturalsort::naturalfactor(finalData$xField)
  finalData$yField <- naturalsort::naturalfactor(finalData$yField)

  #Do statistics
  
  #Initial Visualization
  ggplot(seuratObj@meta.data) + 
    geom_boxplot(aes_string(x = independentVariableTestField, y = dependentVariableTestField))
  
  #Fit linear model relating "independent variable" and SizeFactor
  #Interpret the first intercept as mean SizeFactor and subsequent intercepts as a difference in average size factor
  model.lm <- lm(paste(dependentVariableTestField, "~", independentVariableTestField), data = seuratObj@meta.data)
  print(model.lm)
  
  #Use this linear model to correlate observed data vs fitted values and correlate 
  summary(model.lm)
  
  #Create correlation plot
  pCorr <- ggplot(data.frame(x = model.lm$fitted.values, y = seuratObj@meta.data[,dependentVariableTestField]), aes(x = x, y = y)) +
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) + 
    ggtitle(paste(paste0("Correlation between Fitted ", independentVariableTestField, " and ", dependentVariableTestField,":") , round(cor(seuratObj@meta.data[,dependentVariableTestField], model.lm$fitted.values),2))) + 
    ylab(paste0("Observed ", dependentVariableTestField)) + 
    xlab(paste0("Fitted ", dependentVariableTestField))
  
  print(pCorr)

  finalData <- finalData[c('xField', 'yField', extraGroupingFields, 'Proportion', 'ClusterProportion')]
  
  return(finalData)
}

#' @title MakeEnrichmentDotPlot
#'
#' @description An extremely overloaded function that calculates statistics and enrichment in Seurat Objects. Please see an example dot plot before using this function. Note the aggregated data can be obtained from the ggplot object (i.e. P1$data)
#' @param seuratObj The seurat object that holds the data.
#' @param xField The x axis for the dotplots.
#' @param yField The y axis for the dotplots.
#' @param colorField The column of metadata that is used for the colorField. This should be a factor with two levels. The highest level, reported by levels(), will be treated as the 'highest' value.
#' @param colorLabels Vector of length 2-3 that defines the extremes and midpoint of the colorField axis. If null, the values of colorField will be used.
#' @param extraGroupingFields An optional vector of additional fields to use in grouping
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 
#' @param independentVariableTestField This and dependentVariableTestField automatically define a small statistical test to see if your size factors are correlated. Ideally, they should not be if independentVariableTestField is not the same value as normalizationField. 
#' @param dependentVariableTestField This value should be equal to sizeFactorField initially, but can be changed to interactively see other correlations in the metadata.
#' @export
MakeEnrichmentDotPlot <- function(seuratObj,
                                  yField = 'ClusterNames_0.2',
                                  xField = 'Timepoint',
                                  colorField = 'Tissue',
                                  colorLabels = NULL,
                                  extraGroupingFields = NULL,
                                  normalizationField = 'cDNA_ID',
                                  sizeFactorField = 'SizeFactor',
                                  maxSizeFactor = 100,
                                  independentVariableTestField = colorField,
                                  dependentVariableTestField = sizeFactorField

){
  
  metacounts <- ConstructEnrichmentDataFrameAndDoStatistics(seuratObj = seuratObj,
                                                            xField = xField,
                                                            yField = yField,
                                                            colorField = colorField,
                                                            extraGroupingFields = extraGroupingFields,
                                                            normalizationField = normalizationField,
                                                            sizeFactorField = sizeFactorField,
                                                            maxSizeFactor = maxSizeFactor,
                                                            independentVariableTestField = independentVariableTestField ,
                                                            dependentVariableTestField = dependentVariableTestField
  )

  if (is.null(colorLabels)) {
    colorLabels <- levels(seuratObj@meta.data[[colorField]])
    fillBreaks <- c(0,1)
  } else {
    if (length(colorLabels) > 3) {
      stop('Expected colorLabels to be either 2 or 3 values')
    }

    fillBreaks <- seq(0,1, 1/(length(colorLabels) - 1))
  }

  P1 <- ggplot(metacounts, aes(y = yField, x = xField)) +
    geom_point(aes(size = ClusterProportion, fill = Proportion), alpha = 1, shape = 21) +
    scale_size_continuous(
      limits = c(0, max(metacounts$ClusterProportion)),
      range = c(0,8),
      breaks = c(0.01,max(metacounts$ClusterProportion)/2,max(metacounts$ClusterProportion)),
      labels = c("1%", paste0(round(max(metacounts$ClusterProportion)/2 *100),"%"), paste0(round(max(metacounts$ClusterProportion) *100),"%"))
    ) +
    scale_fill_gradient2(low = "cadetblue2", high = "red", limits = c(0,1), midpoint = .5, breaks = fillBreaks, labels = colorLabels) +
    scale_y_discrete(limits = rev(levels(metacounts$yField))) +
    theme(
      axis.text.x = element_text(color = "grey20", size = 16, angle = 45, hjust = 1, vjust = 1, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 16, face = "plain"),
      axis.title.x = element_text(color = "grey20", size = 20, face = "plain"),
      axis.title.y = element_text(color = "grey20", size = 20, angle = 90, face = "plain"),
      panel.background = element_rect(fill = NA)
    ) +
    xlab(xField) +
    ylab(yField) +
    guides(fill = guide_colorbar(order = 1))
  
  return(P1)
}
