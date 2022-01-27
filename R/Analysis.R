#' @import ggplot2 Seurat dplyr
#' @importFrom stats cor 


#' @title .GenerateSizeFactor 
#'
#' @description A helper function for .ConstructDataFrameAndDoStatistics that calculates size factors for differently sized scRNA-Seq datasets.
#' @param seuratObj The seurat object that holds the data.
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 

.GenerateSizeFactor <- function(seuratObj, normalizationField = 'cDNA_ID', sizeFactorField = 'SizeFactor', maxSizeFactor = 100){
  #check max size factor for data
  dataMaxSizeFactor <- max(table(seuratObj@meta.data[,normalizationField]))/min(table(seuratObj@meta.data[,normalizationField]))
  if(!maxSizeFactor >= dataMaxSizeFactor){
    print(paste0("Max size factor: ", round(dataMaxSizeFactor,2)))
    stop()
  }
  #Calculate the ratio between maximum size 
  print("Calculating size factors")
  for (dataset in unique(seuratObj@meta.data[,normalizationField])){
    datasetSizeFactor <- max(table(seuratObj@meta.data[,normalizationField]))/ table(seuratObj@meta.data[seuratObj@meta.data[,normalizationField] == dataset ,normalizationField])
    seuratObj@meta.data[seuratObj@meta.data[,normalizationField] == dataset , sizeFactorField] <- round(datasetSizeFactor, 4)
  }
  return(seuratObj)
}

#' @title .ConstructDataFrameAndDoStatistics
#'
#' @description Helper function for makeDotPlot that does the dataframe construction and statistics. 
#' @param seuratObj The seurat object that holds the data.
#' @param colorValue The "value that should be ranked "High" on the color axis.
#' @param colorField The column of metadata that is used for the colorField
#' @param colorLabels Vector of length 3 that defines the extremes and midpoint of the colorField axis.
#' @param yField The y axis for the dotplots.
#' @param groupField The x axis for the dotplots.
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 
#' @param independentVariableTestField This and dependentVariableTestField automatically define a small statistical test to see if your size factors are correlated. Ideally, they should not be if independentVariableTestField is not the same value as normalizationField. 
#' @param dependentVariableTestField This value should be equal to sizeFactorField initially, but can be changed to interactively see other correlations in the metadata.
#' @export

.ConstructDataFrameAndDoStatistics <- function(seuratObj, colorField = 'Population', yField = 'ClusterNames_0.2', groupField = 'Timepoint', colorValue = "ld-LN-Right", colorLabels = c("Left", "Even", "Right"), normalizationField = 'cDNA_ID', sizeFactorField = 'SizeFactor', maxSizeFactor = 100, independentVariableTestField = "Tissue", dependentVariableTestField = "SizeFactor"){
  
  # Calculate Size Factor
  seuratObj <- .GenerateSizeFactor(seuratObj, normalizationField = normalizationField, sizeFactorField = sizeFactorField, maxSizeFactor = maxSizeFactor)
  
  # Construct Data Frame
  metatibble <- as_tibble(seuratObj@meta.data[c(colorField, yField, groupField, 'Vaccine', normalizationField, sizeFactorField)])
  names(metatibble) <- c('colorField', 'yField', 'groupField', 'Vaccine', normalizationField, sizeFactorField)
  
  metatibble$groupField <- as.character(metatibble$groupField)
  
  #initial cluster enrichment tibble
  metacounts <- metatibble %>% dplyr::count(yField, groupField, Vaccine, wt = SizeFactor)
  print(metacounts)
  #parallel count tibble to track lung/tissue enrichment
  metacounts.tissue <- metatibble %>% dplyr::count(colorField, groupField, yField, Vaccine, wt = SizeFactor) %>% group_by(colorField, groupField, yField,Vaccine,n) %>% group_by(groupField, yField, Vaccine) %>% mutate(TissueProportion = (n/sum(n)))
  
  #bisect tissue tibble to report only one tissue
  metacounts.tissue <- metacounts.tissue[metacounts.tissue$colorField == colorValue,]
  
  #calulate cluster enrichment
  metacounts <- metacounts %>% group_by(yField, Vaccine) %>% mutate(ClusterProportion = prop.table(n))
  
  #merge cluster enrichment and lung enrichment tibbles
  metacounts <- merge(metacounts, metacounts.tissue[ , c("groupField", "yField", "Vaccine", "TissueProportion")], by = c("groupField", "Vaccine", "yField"), )
  metacounts$groupField <- naturalsort::naturalfactor(metacounts$groupField)
  metacounts$yField <- naturalsort::naturalfactor(metacounts$yField)
  
  
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
  pCorr <- ggplot(data.frame(x = model.lm$fitted.values, y = seuratObj@meta.data[,dependentVariableTestField]), aes(x = x, y = y))+
    geom_point() + 
    geom_smooth(method='lm', formula= y~x) + 
    ggtitle(paste(paste0("Correlation between Fitted ", independentVariableTestField, " and ", dependentVariableTestField,":") , round(cor(seuratObj@meta.data[,dependentVariableTestField], model.lm$fitted.values),2))) + 
    ylab(paste0("Observed ", dependentVariableTestField)) + 
    xlab(paste0("Fitted ", dependentVariableTestField))
  
  print(pCorr)
  
  return(metacounts)
}

#' @title makeDotPlot
#'
#' @description An extremely overloaded function that calculates statistics and enrichment in Seurat Objects. Please see an example dot plot before using this function. 
#' @param seuratObj The seurat object that holds the data.
#' @param colorValue The "value that should be ranked "High" on the color axis.
#' @param colorField The column of metadata that is used for the colorField
#' @param colorLabels Vector of length 3 that defines the extremes and midpoint of the colorField axis.
#' @param yField The y axis for the dotplots.
#' @param groupField The x axis for the dotplots.
#' @param normalizationField The metadata column that is used for size factor calculation (equivalently, normalization between different datasets.) Recommended to be cDNA_ID.
#' @param sizeFactorField The column name of the seruat object metadata that size factors should be stored in.
#' @param maxSizeFactor The maximum allowable SizeFactor before an error is automatically thrown. For instance, a size factor of 100 means you have a group of 8000 cells that you're comparing to 80 cells. 
#' @param independentVariableTestField This and dependentVariableTestField automatically define a small statistical test to see if your size factors are correlated. Ideally, they should not be if independentVariableTestField is not the same value as normalizationField. 
#' @param dependentVariableTestField This value should be equal to sizeFactorField initially, but can be changed to interactively see other correlations in the metadata.
#' @export

makeDotPlot <- function(seuratObj, colorValue = "Lung-Right", colorField = 'Tissue', colorLabels = c("Left", "Even", "Right"), yField = 'ClusterNames_0.2', groupField = 'Timepoint', normalizationField = 'cDNA_ID', sizeFactorField = 'SizeFactor', maxSizeFactor = 100, independentVariableTestField = "Tissue", dependentVariableTestField = "SizeFactor"){
  
  metacounts <- .ConstructDataFrameAndDoStatistics(seuratObj = seuratObj, colorValue = colorValue, colorField = colorField, yField = yField, groupField = groupField, normalizationField = normalizationField, sizeFactorField = sizeFactorField, maxSizeFactor = maxSizeFactor, independentVariableTestField = independentVariableTestField , dependentVariableTestField = dependentVariableTestField)
  
  p1 <- ggplot(metacounts, aes(y = yField, x = groupField)) +
    geom_point(aes(size = ClusterProportion, fill = TissueProportion), alpha = 1, shape = 21) +
    scale_size_continuous(limits = c(0, max(metacounts$ClusterProportion)), range = c(0,8), breaks = c(0.01,max(metacounts$ClusterProportion)/2,max(metacounts$ClusterProportion)) ,labels = c("1%", paste0(round(max(metacounts$ClusterProportion)/2 *100),"%"), paste0(round(max(metacounts$ClusterProportion) *100),"%"))) +
    scale_fill_gradient2(low = "cadetblue2", high = "red", limits = c(0,1), midpoint = .5, breaks = seq(0,1,.5), labels = colorLabels) +
    scale_y_discrete(limits = rev(levels(metacounts$yField))) +
    theme(
      axis.text.x = element_text(color = "grey20", size = 16, angle = 45, hjust = 1, vjust = 1, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 16, face = "plain"),
      axis.title.x = element_text(color = "grey20", size = 20, face = "plain"),
      axis.title.y = element_text(color = "grey20", size = 20, angle = 90, face = "plain"),
      panel.background = element_rect(fill = NA)
    ) +
    ylab("Cluster") +
    xlab("")+
    labs(fill = "Tissue Bias", size = "Cluster Bias") +
    guides(fill = guide_colorbar(order = 1)) +
    facet_grid(. ~ Vaccine)
  
  
  return(p1)
}