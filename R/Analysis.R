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
  rawData <- rawData %>% tidyr::unite("XY_Key", tidyr::all_of(c('xField', 'yField', extraGroupingFields)), remove = FALSE)
  rawData <- rawData %>% tidyr::unite("Y_Key", tidyr::all_of(c('yField', extraGroupingFields)), remove = FALSE)
  
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
  metadata <- unique(rawData[c('XY_Key', 'xField', 'yField', extraGroupingFields)])
  finalData <- merge(finalData, metadata, by = "XY_Key")
  finalData$xField <- naturalsort::naturalfactor(finalData$xField)
  finalData$yField <- naturalsort::naturalfactor(finalData$yField)

  #Do statistics
  
  #Initial Visualization
  ggplot(seuratObj@meta.data) + 
    geom_boxplot(aes(x = !!rlang::sym(independentVariableTestField), y = !!rlang::sym(dependentVariableTestField)))
  
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

#' @title CalculateClusterEnrichment
#'
#' @description A function that calculates the enrichment of a cluster under a given treatment variable. 
#' @param seuratObj The Seurat object containing a subjectField, clusterField, and treatmentField. Please see the individual arguments for more information.
#' @param subjectField The column of the Seurat object's metadata that contains the subject field. This field should denote individual samples that are independently collected.
#' @param clusterField The column of the Seurat object's metadata that contains the clustering field. This field should denote cluster membership, generally given by louvain/leiden clustering, but any subject-independent clustering method is valid.  
#' @param treatmentField The column of the Seurat object's metadata that contains the treatment field. This field should denote the treatment of the subject, and should be the primary variable of interest within your study. 
#' @param alternative A passthrough variable to wilcox.test. If "greater", the alternative hypothesis is that the difference in medians is greater than the null hypothesis. If "less", the alternative hypothesis is that the difference in medians is less than the null hypothesis. If "two.sided", the alternative hypothesis is that the difference in medians is simply "different" from the null hypothesis. In the case of the wilcoxon rank sum (e.g. paired = FALSE), this will test the difference of the medians, rather than the medians themselves. 
#' @param pValueCutoff The p-value cutoff for significance.
#' @param showPlots A boolean that determines if the cluster significance should be shown in a DimPlot.
#' @param paired A passthrough variable to wilcox.test. If TRUE, the function will perform a paired Wilcoxon test. If FALSE, the function will perform an unpaired Wilcoxon test. If you're testing (for instance, timepoint) enrichment on repeated measures, this should be TRUE. If you're testing different treatments on different subjects, this should be FALSE. If set to "infer", the function will attempt to infer the correct value based on the name of the treatment field. Specifically, this will search for "time" in your treatment field. If it finds it, it will set paired = TRUE. If it doesn't, it will set paired = FALSE.
#' @param removePriorPvalues A boolean that determines if the prior p-values should be removed from the Seurat object metadata. It's likely that you'll want to iteratively compute significance on different metadata fields, so this is set to TRUE by default and will remove the Cluster_pValue and Cluster_p_adj fields from the Seurat object's metadata.
#' @param postHocTest A boolean that determines if a post-hoc test should be performed. If TRUE, the function will perform a Conover-Iman post-hoc test to determine which pairs of treatmentField groups are significantly different from each other.
#' @return A Seurat object with the p-values of the clusters in the metadata columns Cluster_pValue and Cluster_p_adj. If showPlots = TRUE, a DimPlot will be shown with significant clusters highlighted.
#' @examples 
#'  \dontrun{
#'  seuratObj <- CalculateClusterEnrichment(seuratObj,
#'                                        clusterField = "ClusterNames_0.4",
#'                                        treatmentField = "vaccine_cohort",
#'                                        subjectField = "SubjectId",
#'                                        paired = "infer", 
#'                                        showPlots = TRUE)
#'                                        }
#' @export

CalculateClusterEnrichment <- function(seuratObj,
                                       subjectField = 'SubjectId',
                                       clusterField = 'ClusterNames_0.2',
                                       treatmentField = NULL,
                                       alternative = 'two.sided',
                                       pValueCutoff = 0.05,
                                       showPlots = TRUE, 
                                       paired = "infer", 
                                       removePriorPvalues = TRUE, 
                                       postHocTest = TRUE
){
  # test for validity of metadata fields within the seurat object
  # treatmentField
  if (is.null(treatmentField)) {
    stop('treatmentField is set to NULL, please specify a valid treatmentField. treatmentField should should denote the treatment of the subject, such as drug or vaccine administration, and should be the primary variable of interest within your study. Timepoint-based metadata fields are also valid values of treatmentField.')
  } else if (!treatmentField %in% colnames(seuratObj@meta.data)) {
    stop(paste0('treatmentField: ', treatmentField, ' not found in the seuratObject metadata columns. Please check the spelling and case sensitivity.'))
  } else if (grepl(" - ", seuratObj@meta.data[[treatmentField]])) {
    stop(paste0('treatmentField: ', treatmentField , ' has entries that contain a " - " character. Please change the entries containing the  " - " characters from treatmentField to any other delimiting values.'))
  }
  # subjectField
  if (!subjectField %in% colnames(seuratObj@meta.data)) {
    stop(paste0('subjectField: ', subjectField, ' not found in the seuratObject metadata columns. For Prime-Seq Seurat objects, this should be "SubjectId". Please check the spelling and case sensitivity.'))
  }
  # clusterField 
  if (!clusterField %in% colnames(seuratObj@meta.data)) {
    stop(paste0('clusterField: ', clusterField, ' not found in the seuratObject metadata columns. For Prime-Seq Seurat objects, this should be "ClusterNames_X" where X is a resolution parameter between 0.2 and 1.2.'))
  }
  # alternative
  if (!alternative %in% c('greater', 'less', 'two.sided')) {
    stop(paste0('alternative: ', alternative, ' is not an valid value. Please use one of: "greater", "less", or "two.sided" based on your definition of "enrichment". See ?wilcox.test for more information.'))
  }
  # pvalueCutoff
  if (!is.numeric(pValueCutoff)) {
    stop(paste0('pValueCutoff: ', pValueCutoff, ' is not a numeric value. Please specify a numeric value for pValueCutoff. This field is only used for plotting a DimPlot if showPlots = TRUE.'))
  }
  # showPlots
  if (!is.logical(showPlots)) {
    stop(paste0('showPlots: ', showPlots, ' is not a boolean. Please specify showPlots = TRUE or showPlots = FALSE. If TRUE, a DimPlot will be shown with significant clusters highlighted.'))
  }
  # paired
  # If the paired variable is set to infer, attempt to infer the correct value based on the name of the treatment field. If it's a timepoint field, it's likely paired.
  if (paired == "infer") {
    if (grepl("time", tolower(treatmentField))) {
      message(paste0("Inferred paired = TRUE based on treatmentField = ", treatmentField))
      paired <- TRUE
    } else {
      message(paste0("Inferred paired = FALSE based on treatmentField = ", treatmentField))
      paired <- FALSE
    }
  }
  # removePriorPvalues
  if (!is.logical(removePriorPvalues)) {
    stop(paste0('removePriorPvalues: ', removePriorPvalues, ' is not a boolean. Please specify removePriorPvalues = TRUE or removePriorPvalues = FALSE. If TRUE, the p-values from prior runs of CalculateClusterEnrichment will be removed from the Seurat object metadata.'))
  }
  
  # Remove prior p-values added by CalculateClusterEnrichment if removePriorPvalues = TRUE.
  if (removePriorPvalues) {
    message(paste0("Removing Cluster_pValue and Cluster_p_adj columns from the metadata in preparation for computing new p-values based on treatmentField: ", treatmentField, " and clusterField: ", clusterField, "."))
    seuratObj$Cluster_pValue <- NULL
    seuratObj$Cluster_p_adj <- NULL
  }
  
  # Calculate the total number of cells in the cluster
  clusterProportionsDataFrame <- seuratObj@meta.data %>% 
    dplyr::group_by(!!rlang::sym(clusterField)) %>%
    dplyr::mutate(ClusterCount = n()) %>% 
    dplyr::group_by(!!rlang::sym(clusterField), !!rlang::sym(subjectField), !!rlang::sym(treatmentField)) %>% 
    dplyr::reframe(SubjectClusterProportion = n()/ClusterCount) %>% 
    unique.data.frame()
  #instantiate a dataframe to store the enrichment statistics
  enrichmentDataFrame <- data.frame()
  
  #check the number of treatment groups and infer a non-parametric test
  if (length(unique(seuratObj@meta.data[[treatmentField]])) == 1) {
    stop(paste0('Only one treatment group found in treatmentField: ', treatmentField, ". Comparisons between groups are impossible with only one group."))
  } else if (length(unique(seuratObj@meta.data[[treatmentField]])) > 2) { 
    message("More than two treatment groups found. Kruskal-Wallis rank sum test will be performed.")
    #loop through each cluster and perform kruskal-wallis test
    for (cluster in unique(clusterProportionsDataFrame[[clusterField]])) { 
      clusterProportionsSubset <- clusterProportionsDataFrame[clusterProportionsDataFrame[[clusterField]] == cluster,]
      clusterProportions <- clusterProportionsSubset$SubjectClusterProportion
      groupMembership <- clusterProportionsSubset[[treatmentField]]
      #calculate p values
      pValue <- kruskal.test(clusterProportions, groupMembership, alternative = alternative)$p.value
      enrichmentDataFrame <- rbind(enrichmentDataFrame, data.frame(clusterField = cluster, Cluster_pValue = pValue))
    }
    #adjust p values
    enrichmentDataFrame$Cluster_p_adj <- p.adjust(enrichmentDataFrame$Cluster_pValue, n = nrow(enrichmentDataFrame))
    #if specified, a post hoc test will be performed on each significant cluster from the Kruskal-Wallis test.
    if (postHocTest) { 
      for (cluster in unique(clusterProportionsDataFrame[[clusterField]])) { 
        if (enrichmentDataFrame[enrichmentDataFrame$clusterField == cluster,]$Cluster_p_adj < pValueCutoff) {
          # do testing 
          pairwise_test <- conover.test::conover.test(unlist(clusterProportionsDataFrame[clusterProportionsDataFrame[[clusterField]] == cluster, "SubjectClusterProportion"]),
                                             unlist(clusterProportionsDataFrame[clusterProportionsDataFrame[[clusterField]] == cluster, treatmentField]))
          pairwise_test_plotting_dataframe <- data.frame(comparisons = pairwise_test$comparisons, 
                                                T_statistic = pairwise_test$`T`, 
                                                P_val_adj = pairwise_test$`P.adjusted`) %>% 
            tidyr::separate(comparisons, into = c("Group1", "Group2"), sep = " - ") %>% 
            dplyr::mutate(stars = dplyr::case_when(P_val_adj < 0.001 ~ "***", 
                                     P_val_adj < 0.01 ~ "**", 
                                     P_val_adj < 0.05 ~ "*", 
                                     TRUE ~ ""))
          enrichmentPlot <- ggplot2::ggplot(pairwise_test_plotting_dataframe, ggplot2::aes(x = Group1, y = Group2, fill = T_statistic)) + 
            ggplot2::geom_tile() + 
            colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2) + 
            ggplot2::geom_text(aes(label=stars), color="black", size=5) +
            egg::theme_article() + 
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
            patchwork::plot_annotation(title = paste0("Cluster ", cluster, " Enrichment"),
                                       subtitle = "Contrast for T statistic is Group1 - Group2")
          if (showPlots){
            print(enrichmentPlot)
          }
        }
      }
    }
    
  } else if (length(unique(seuratObj@meta.data[[treatmentField]])) == 2) {
    if (paired) { 
      message("Two treatment groups found and paired data was either specified or inferred. Pairwise wilcoxon signed rank test will be performed.")
      #perform wilcoxon signed rank test for each cluster 
      for (cluster in unique(clusterProportionsDataFrame[[clusterField]])) { 
        clusterProportionsSubset <- clusterProportionsDataFrame[clusterProportionsDataFrame[[clusterField]] == cluster,]
        groupOneProportions <- clusterProportionsSubset[clusterProportionsSubset[[treatmentField]] == unique(clusterProportionsSubset[[treatmentField]])[1],]$SubjectClusterProportion
        groupTwoProportions <- clusterProportionsSubset[clusterProportionsSubset[[treatmentField]] == unique(clusterProportionsSubset[[treatmentField]])[2],]$SubjectClusterProportion
        #if a subject is missing from one of the groups, we need to impute the missing percentage with a zero. 
        #I can't figure out a way to pass a "TreatmentGroup by Subjects" table succinctly to dplyr::complete, but it will be really obvious at this stage, so we can impute. 
        if (length(groupOneProportions) != length(groupTwoProportions)) {
          if (length(groupOneProportions) > length(groupTwoProportions)) {
            groupTwoProportions <- c(groupTwoProportions, rep(0, length(groupOneProportions) - length(groupTwoProportions)))
          } else {
            groupOneProportions <- c(groupOneProportions, rep(0, length(groupTwoProportions) - length(groupOneProportions)))
          }
        }
        #calculate p values. 
        pValue <- wilcox.test(groupOneProportions, groupTwoProportions, paired = TRUE, alternative = alternative)$p.value
        enrichmentDataFrame <- rbind(enrichmentDataFrame, data.frame(clusterField = cluster, Cluster_pValue = pValue))
      }
      #adjust p values
      enrichmentDataFrame$Cluster_p_adj <- p.adjust(enrichmentDataFrame$Cluster_pValue, n = nrow(enrichmentDataFrame))
    } else {
      message("Two treatment groups found and unpaired data was either specified or inferred. Wilcoxon rank sum test will be performed.")
      #perform Wilcoxon rank sum test for each cluster 
      for (cluster in unique(clusterProportionsDataFrame[[clusterField]])) { 
        clusterProportionsSubset <- clusterProportionsDataFrame[clusterProportionsDataFrame[[clusterField]] == cluster,]
        groupOneProportions <- clusterProportionsSubset[clusterProportionsSubset[[treatmentField]] == unique(clusterProportionsSubset[[treatmentField]])[1],]$SubjectClusterProportion
        groupTwoProportions <- clusterProportionsSubset[clusterProportionsSubset[[treatmentField]] == unique(clusterProportionsSubset[[treatmentField]])[2],]$SubjectClusterProportion
        #if a subject is missing from one of the groups, we need to impute the missing percentage with a zero. 
        #I can't figure out a way to pass a "TreatmentGroup by Subjects" table succinctly to dplyr::complete, but it will be really obvious at this stage, so we can impute. 
        if (length(groupOneProportions) != length(groupTwoProportions)) {
          if (length(groupOneProportions) > length(groupTwoProportions)) {
            groupTwoProportions <- c(groupTwoProportions, rep(0, length(groupOneProportions) - length(groupTwoProportions)))
          } else {
            groupOneProportions <- c(groupOneProportions, rep(0, length(groupTwoProportions) - length(groupOneProportions)))
          }
        }
        #calculate p values. 
        pValue <- wilcox.test(groupOneProportions, groupTwoProportions, paired = FALSE, alternative = alternative)$p.value
        enrichmentDataFrame <- rbind(enrichmentDataFrame, data.frame(clusterField = cluster, Cluster_pValue = pValue))
      }
      #adjust p values
      enrichmentDataFrame$Cluster_p_adj <- p.adjust(enrichmentDataFrame$Cluster_pValue, n = nrow(enrichmentDataFrame))
    }
  }
  #add cell barcodes to metadata if they aren't present
  if (!"CellBarcode" %in% colnames(seuratObj@meta.data)) { 
    message("CellBarcode not found in metadata. Adding CellBarcode to metadata.")
    seuratObj@meta.data$CellBarcode <- rownames(seuratObj@meta.data)
  }
  #merge the enrichmentDataFrame with the Seurat object metadata to populate p values in the metadata
  metadata <- merge(seuratObj@meta.data, enrichmentDataFrame[, c("Cluster_pValue", "Cluster_p_adj", "clusterField")], by.x = clusterField, by.y =  'clusterField')
  rownames(metadata) <- metadata$CellBarcode
  #populate P values
  seuratObj <- AddMetaData(seuratObj, metadata = metadata)
  
  #show a DimPlot of significant clusters
  if (showPlots) { 
    print(Seurat::DimPlot(seuratObj, cells.highlight = metadata$CellBarcode[metadata$Cluster_p_adj < pValueCutoff]) + ggplot2::ggtitle("Significant Clusters"))
    }
  return(seuratObj)
}




