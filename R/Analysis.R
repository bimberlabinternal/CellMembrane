#' @import ggplot2 Seurat dplyr
#' @importFrom stats cor 
#' @importFrom NBZIMM glmm.nb glmm.zinb

utils::globalVariables(
  names = c('ClusterProportion', 'Proportion', 'SizeFactor', 'XY_Key', 'Y_Key', 'ClusterCount',  
            'comparisons', 'T_statistic', 'P_val_adj', 'Group1', 'Group2', 'stars', 'pct.exp', 'features.plot', 
            'lane_yield', 'Treatment', 'ClusterEnrichment_LFC', 'p.adj', 'Cluster', 'EnrichedPairs'),
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

#' @title CalculateClusterEnrichmentOmnibus
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
#'  seuratObj <- CalculateClusterEnrichmentOmnibus(seuratObj,
#'                                        clusterField = "ClusterNames_0.4",
#'                                        treatmentField = "vaccine_cohort",
#'                                        subjectField = "SubjectId",
#'                                        paired = "infer", 
#'                                        showPlots = TRUE)
#'                                        }
#' @export

CalculateClusterEnrichmentOmnibus <- function(seuratObj,
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
  } else if (any(grepl(" - ", seuratObj@meta.data[[treatmentField]]))) {
    stop(paste0('treatmentField: ', treatmentField , ' has entries (', paste0(unique(seuratObj@meta.data[grepl(" - ", seuratObj@meta.data[[treatmentField]]), treatmentField]), collapse = ', '), ') that contain a " - " character. Please change the entries containing the  " - " characters from treatmentField to any other delimiting values.'))
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
      pValue <- stats::kruskal.test(clusterProportions, groupMembership, alternative = alternative)$p.value
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

#' @title CalculateClusterEnrichmentPairwise
#' 
#' @description A function that calculates the enrichment of a cluster under a given treatment variable using a Generalized Linear Mixed Model (GLMM) or General Linear Model with quasi-likelihood. This allows for more complex treatment of enrichment where you want to account for subject-driven random effects and variable cellular recovery. 
#' @param seuratObj The Seurat object containing a subjectField, clusterField, and treatmentField. Please see the individual arguments for more information.
#' @param subjectField The column of the Seurat object's metadata that contains the subject field. This field should denote individual samples that are independently collected and represent biological replicates.
#' @param clusterField The column of the Seurat object's metadata that contains the clustering field. This field should denote cluster membership, generally given by louvain/leiden clustering, but any subject-independent clustering method is valid.
#' @param biologicalReplicateGroupingVariables A vector of metadata fields that should be used to group biological replicates. This can identical to subjectField, but in the case of repeated sample collections, the variables specified here will uniquely identify samples/collections.
#' @param treatmentField The column of the Seurat object's metadata that contains the treatment field. This field should denote the treatment of the subject, and should be the primary variable of interest within your study.
#' @param referenceValue The reference value for the treatment field. This is the value that will be used as the baseline for comparison. Since this method requires a baseline, a reference value must be provided. Consider using CalculateClusterEnrichment for an omnibus test style of enrichment if a change from an individual reference value is not particularly insightful. 
#' @param pValueCutoff The p-value cutoff for significance. This is used to determine which clusters are significantly enriched or depleted.
#' @param showPlots A boolean that determines if the function should print Volcano and DimPlots. 
#' @param returnSeuratObjectOrPlots A string that determines if the function should return a Seurat object with the enrichment results or a list of plots. If "SeuratObject", the function will return a Seurat object with the enrichment results in the metadata. If "Plots", the function will return a list of volcano plots and the dataframe used to plot them. 
#' @param includeDepletions A boolean that determines if the function should include depletions in the enrichment results. If TRUE, the function will include depletions in the enrichment results. If FALSE, the function will only include enrichments. Since we typically recover variable numbers of cells, depletions can be difficult to interpret/easily explained by low sampling volume compared to transcriptional complexity, so this is set to FALSE by default.
#' @param pvalueAdjustmentMethod The method to use for p-value adjustment. This is passed to the p.adjust function. Default is "holm", but can be set to "BH" or "bonferroni" if desired.
#' @param lowSampleSizeDetection A boolean that determines if the function should detect low sample sizes and warn the user. If TRUE, the function will warn the user if the sample size is below a the threshold specified by lowSampleSizeThreshold. This is useful for detecting low sample sizes that may prevent random effects from being estimated in a mixed model. 
#' @param lowSampleSizeThreshold The threshold for low sample size detection. If the number of samples in a group is below this threshold, the function will warn the user and use a typical GLM rather than a hierarchical model. Default is 3, but can be set to a different value if desired.
#' 
#' @examples
#'  \dontrun{
#'  # Define dummy variables
#'  seuratObj$cDNA_ID <- rep(1:4, ncol(seuratObj)) 
#'  seuratObj$Vaccine <- rep(c("Vaccine1", "Vaccine2"), each = ncol(seuratObj)/2) 
#'  seuratObj$SubjectId <- rep(1:4, each = ncol(seuratObj)/4) 
#'  
#'  seuratObj <- CalculateClusterEnrichmentPairwise(seuratObj,
#'                                                  subjectField = 'SubjectId',
#'                                                  clusterField = 'ClusterNames_0.2',
#'                                                  biologicalReplicateGroupingVariables = c("cDNA_ID"),
#'                                                  treatmentField = "Vaccine",
#'                                                  referenceValue = "Vaccine1",
#'                                                  pValueCutoff = 0.05,
#'                                                  showPlots = TRUE, 
#'                                                  returnSeuratObjectOrPlots = "SeuratObject", 
#'                                                  includeDepletions = FALSE, 
#'                                                  lowSampleSizeDetection = TRUE, 
#'                                                  lowSampleSizeThreshold = 3)
#'   DimPlot(seuratObj, group.by = "Cluster_Enrichment", label = TRUE) + NoLegend()
#'   }
#'  
#' @export
#' @return A SeuratObject or a list of plots + dataframe, depending on the value of returnSeuratObjectOrPlots.

CalculateClusterEnrichmentPairwise <- function(seuratObj,
                                               subjectField = 'SubjectId',
                                               clusterField = 'ClusterNames_0.2',
                                               biologicalReplicateGroupingVariables = c("cDNA_ID"),
                                               treatmentField = NULL,
                                               referenceValue = NULL,
                                               pValueCutoff = 0.05,
                                               showPlots = TRUE, 
                                               returnSeuratObjectOrPlots = "SeuratObject", 
                                               includeDepletions = FALSE, 
                                               pvalueAdjustmentMethod = "holm", 
                                               lowSampleSizeDetection = TRUE, 
                                               lowSampleSizeThreshold = 3) {
  #basic checks for valid inputs
  if (!inherits(seuratObj, "Seurat")) {
    stop("seuratObj must be a Seurat object.")
  } else if (typeof(subjectField) != "character") {
    stop("subjectField must be a character string representing the column name in the Seurat object's metadata that contains the subject IDs.")
  } else if (typeof(clusterField) != "character") {
    stop("clusterField must be a character string representing the column name in the Seurat object's metadata that contains the cluster names.")
  } else if (is.null(treatmentField)) {
    stop("treatmentField must be specified as a character string representing the column name in the Seurat object's metadata that contains the treatment information.")
  } else if (typeof(treatmentField) != "character") {
    stop("treatmentField must be a character string representing the column name in the Seurat object's metadata that contains the treatment information.")
  } else if (!is.null(referenceValue) && typeof(referenceValue) != "character") {
    stop("referenceValue must be a character string representing the reference value for the treatment field.")
  } else if (!is.numeric(pValueCutoff) || pValueCutoff <= 0 || pValueCutoff >= 1) {
    stop("pValueCutoff must be a numeric value between 0 and 1.")
  } else if (!is.logical(showPlots)) {
    stop("showPlots must be a boolean value (TRUE or FALSE).")
  } else if (!returnSeuratObjectOrPlots %in% c("SeuratObject", "Plots")) {
    stop("returnSeuratObjectOrPlots must be either 'SeuratObject' or 'Plots'.")
  } else if (!is.logical(includeDepletions)) {
    stop("includeDepletions must be a boolean value (TRUE or FALSE).")
  } else if (!is.character(pvalueAdjustmentMethod) || !pvalueAdjustmentMethod %in% stats:::p.adjust.methods) {
    stop("pvalueAdjustmentMethod must be one of 'holm', 'BH', or 'bonferroni'.")
  } else if (!is.logical(lowSampleSizeDetection)) {
    stop("lowSampleSizeDetection must be a boolean value (TRUE or FALSE).")
  } else if (!is.numeric(lowSampleSizeThreshold) || lowSampleSizeThreshold <= 0) {
    stop("lowSampleSizeThreshold must be a positive numeric value.")
  }
  
  #more specific checks to make sure modeling is possible
  #subjectField exists & is coerce-able to a non-degenerate factor
  if (!(subjectField %in% colnames(seuratObj@meta.data))) {
    stop(paste0("subjectField: ", subjectField, " not found in the Seurat object metadata. Please check the spelling and case sensitivity."))
  } else if (length(levels(factor(unique(seuratObj@meta.data[,subjectField])))) <= 1) {
    stop(paste0(subjectField, " must be a field with more than one unique value. Please check the values in the " , subjectField, " column."))
  }
  #treatment field exists and is coerce-able to a non-degenerate factor
  if (!(treatmentField %in% colnames(seuratObj@meta.data))) {
    stop(paste0("treatmentField: ", treatmentField, " not found in the Seurat object metadata. Please check the spelling and case sensitivity."))
  } else if (length(levels(factor(unique(seuratObj@meta.data[,treatmentField])))) <= 1) {
    stop(paste0(treatmentField, " must be a field with more than one unique value. Please check the values in the " , treatmentField, " column."))
  }
  
  #initialize metadata
  grouping_variables <- c(biologicalReplicateGroupingVariables, subjectField, treatmentField, clusterField)
  
  #coerce metadata
  metadata <- seuratObj@meta.data %>% 
    dplyr::group_by(dplyr::across(dplyr::all_of(biologicalReplicateGroupingVariables))) %>% 
    dplyr::mutate(lane_yield = n()) %>% 
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_variables))) %>% 
    dplyr::mutate(total = dplyr::n(), lane_yield = lane_yield) %>% 
    dplyr::select(dplyr::any_of(c(grouping_variables, "total", "lane_yield"))) %>% 
    dplyr::distinct() %>% 
    dplyr::ungroup() %>% 
    tidyr::complete(
      tidyr::nesting(!!sym(subjectField),!!!syms(biologicalReplicateGroupingVariables), !!sym(treatmentField), lane_yield),
      !!sym(clusterField),
      fill = list(
        total = 0
      )
    )
  
  #force the reference level
  if (!is.null(referenceValue)) {
    metadata[[treatmentField]] <- forcats::fct_relevel(metadata[[treatmentField]], referenceValue, after = 0) 
  } else {
    stop("referenceValue must be specified as a character string representing the reference value for the treatment field. This is used to set the 'background' or 'control' level for the treatment field factor. Intercept-less modeling is not supported in this wrapper.")
  }
  
  lowSampleSize <- FALSE
  if (lowSampleSizeDetection) {
    #check for low sample sizes
    if (length(unlist(unique(metadata[, subjectField]))) < lowSampleSizeThreshold) {
      lowSampleSize <- TRUE
      warning(paste0("There are fewer than ", lowSampleSizeThreshold, " subjects in this dataset. This may lead to unreliable results, especially if the number of cells per subject is low. To adjust for this, subject-specific cluster bias will not be estimated, in favor of a simpler model. If you would like to use a more complex model, please ensure that you have at enough replication within the experiment to support subject-specific effects. Alternatively, you can set lowSampleSizeDetection = FALSE or raise the lowSampleSizeThreshold argument, but singularities may occur that prevent model fitting."))
    }
  }
  
  
  #fit models
  results <- data.frame()
  for (cluster in unique(unlist(metadata[, clusterField]))) {
    metadata_subset <- metadata %>% filter(!!sym(clusterField) == cluster)
    
    if (lowSampleSize) {
      #if low sample size, use a simple quasipoisson model
      model_type <- "Low Sample Size QuasiPoisson"
      model <- stats::glm(
        stats::as.formula(
          paste0("total ~ ", treatmentField, " + offset(log(lane_yield))")
        ),
        data = metadata_subset,
        family = stats::quasipoisson()
      )
      
      #store results for GLMs & handle missing levels.
      if (!all(unique(unlist(metadata[,treatmentField])) %in% gsub(paste0("^",treatmentField), "", names(model$coefficients)))) {
        #the reference level will always be missing, as it's the intercept, so just take the non-first elements
        missing_coefs <- unique(unlist(metadata[,treatmentField]))[!unique(unlist(metadata[,treatmentField])) %in% gsub(paste0("^",treatmentField), "", names(model$coefficients))]
        coefs <- model$coefficients
        errors <- sqrt(diag(sandwich::vcovHC(model, type="HC3")))
        pvalues <- stats::coef(summary(model))[,4]
        for (missing_coef in missing_coefs) { 
          
          if (missing_coef == referenceValue) {
            next #skip the reference level
          }
          coefs[paste0(treatmentField, missing_coef)] <- NA
          errors[paste0(treatmentField, missing_coef)] <- NA
          pvalues[paste0(treatmentField, missing_coef)] <- 1
        }
        
      } else {
        coefs <- model$coefficients
        errors <- sqrt(diag(sandwich::vcovHC(model, type="HC3")))
        pvalues <- stats::coef(summary(model))[,4]
      }
      
    } else {
      model <- tryCatch({
        model_type <- "Zero-Inflated Negative Binomial"
        model <- NBZIMM::glmm.zinb(
          stats::as.formula(
            paste0("total ~ ", treatmentField, " + offset(log(lane_yield))")
          ),
          data = metadata_subset,
          random = stats::as.formula(paste0("~ 1|", subjectField)),
          zi_fixed = ~ 1,
          zi_random = NULL
        )
        
      }, error = function(e) {
        model_type <- "Negative Binomial"
        model <- NBZIMM::glmm.nb(
          stats::as.formula(
            paste0("total ~ ", treatmentField, " + offset(log(lane_yield))")
          ),
          data = metadata_subset,
          random = stats::as.formula(paste0("~ 1|", subjectField))
        )
        return(model)
      })
      
      #store results for GLMMs & handle missing levels. 
      if (!all(unique(unlist(metadata[,treatmentField])) %in% gsub(paste0("^",treatmentField), "", names(model$coefficients$fixed)))) {
        #the reference level will always be missing, as it's the intercept, so just take the non-first elements
        missing_coefs <- unique(unlist(metadata[,treatmentField]))[!unique(unlist(metadata[,treatmentField])) %in% gsub(paste0("^",treatmentField), "", names(model$coefficients$fixed))]
        coefs <- model$coefficients$fixed
        #try to compute robust standard errors, but fall back to standard errors if unable
        errors <- tryCatch(
          {
            sqrt(diag(clubSandwich::vcovCR(model, type = "CR3")))
          },
          error = function(e) {
            message("\nProblem in estimating robust standard errors in Cluster: ", cluster,  " using clubSandwich::vcovCR. \n Returning standard errors.")
            sqrt(diag(vcov(model)))
          }
        )
        pvalues <- summary(model)$tTable[,5]
        for (missing_coef in missing_coefs) { 
          
          if (missing_coef == referenceValue) {
            next #skip the reference level
          }
          coefs[paste0(treatmentField, missing_coef)] <- NA
          errors[paste0(treatmentField, missing_coef)] <- NA
          pvalues[paste0(treatmentField, missing_coef)] <- 1
        }
        
      } else {
        coefs <- model$coefficients$fixed
        #try to compute robust standard errors, but fall back to standard errors if unable
        errors <- tryCatch(
          {
            sqrt(diag(clubSandwich::vcovCR(model, type = "CR3")))
          },
          error = function(e) {
            message("\nProblem in estimating robust standard errors in Cluster: ", cluster,  " using clubSandwich::vcovCR. \n Returning standard errors.")
            sqrt(diag(vcov(model)))
          }
        )
        pvalues <- summary(model)$tTable[,5]
      }
    }
    
    #store results into dataframe
    results <- rbind(results, 
                     data.frame(Cluster = cluster, 
                                Treatment =  gsub("\\(Intercept\\)", referenceValue, gsub(paste0("^", treatmentField),"", names(coefs))),
                                ClusterEnrichment_LFC = coefs[grepl(paste0(treatmentField, "|Intercept"), names(coefs))], 
                                StdError = errors[grepl(paste0(treatmentField, "|Intercept"), names(errors))],
                                pValue = pvalues[grepl(paste0(treatmentField, "|Intercept"), names(pvalues))], 
                                model_type = model_type
                     )
    )
    
  }
  
  #adjust p values
  results$p.adj <- p.adjust(results$pValue, method = pvalueAdjustmentMethod)
  
  #construct volcano plot
  volcano_plot <- results %>% 
    filter(Treatment != referenceValue) %>% 
    ggplot2::ggplot(ggplot2::aes(x = ClusterEnrichment_LFC, y = -log10(p.adj), 
                                 color = Treatment, 
                                 shape = model_type)) + 
    ggplot2::geom_point() + 
    ggrepel::geom_label_repel(aes(label = Cluster), show.legend = F) + 
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    egg::theme_article() + 
    ggplot2::labs(x = "Log Fold Change", y = "-log10 Adjusted P-value", 
                  title = paste0("Differential Abundance of Clusters by ", treatmentField), 
                  subtitle = paste0("Estimation of ", treatmentField, " enrichment in ", clusterField, " clustering")) 
  #optionally edit volcano plot to only show positive values
  if (!includeDepletions) {
    
    only_enrichments <- results %>% 
      dplyr::filter(ClusterEnrichment_LFC > 0) %>% 
      dplyr::filter(Treatment != referenceValue)
    
    
    volcano_plot <- volcano_plot + 
      ggplot2::xlim(c(0, max(only_enrichments %>% pull(ClusterEnrichment_LFC), na.rm = TRUE)+1)) + 
      ggplot2::ylim(c(0, max(-log10(only_enrichments %>% dplyr::pull(p.adj)), na.rm = TRUE)+1))
  }
  if (showPlots) {
    print(volcano_plot)
  }
  #returns
  if (returnSeuratObjectOrPlots == "SeuratObject") { 
    #add the results to the Seurat object metadata
    metadata_to_add <- seuratObj@meta.data
    
    filtered_results <- results %>% 
      dplyr::filter(p.adj < pValueCutoff) %>%
      dplyr::filter(Treatment != referenceValue) %>%
      dplyr::mutate(EnrichedPairs = paste0(Treatment, ":", Cluster))
    
    #check if ClusterEnrichment_LFC already exists in the seurat object's metadata, and if so, remove it before the join. 
    if ("ClusterEnrichment_LFC" %in% colnames(metadata_to_add)) {
      metadata_to_add <- metadata_to_add %>% dplyr::select(-ClusterEnrichment_LFC)
    }
    
    metadata_to_add <- metadata_to_add %>% 
      dplyr::mutate(Enrichment_VariablePair = paste0(!!sym(treatmentField), ":", !!sym(clusterField))) %>%
      dplyr::left_join(filtered_results %>% 
                         dplyr::select(EnrichedPairs, ClusterEnrichment_LFC), 
                       by = c("Enrichment_VariablePair" = "EnrichedPairs")) %>%
      dplyr::mutate(Cluster_Enrichment = case_when(
        !is.na(ClusterEnrichment_LFC) & ClusterEnrichment_LFC > 0 ~ paste0("Enriched: ", Enrichment_VariablePair),
        !is.na(ClusterEnrichment_LFC) & ClusterEnrichment_LFC < 0 ~ paste0("Depleted: ", Enrichment_VariablePair),
        TRUE ~ "Not Enriched"
      ))
    
    if (!includeDepletions) {
      metadata_to_add$Cluster_Enrichment <- ifelse(grepl("Depleted", metadata_to_add$Cluster_Enrichment), 
                                                   "Not Enriched", 
                                                   metadata_to_add$Cluster_Enrichment)
    }
    #append metadata
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = metadata_to_add)
    
    if (showPlots) {
      print(Seurat::DimPlot(seuratObj, group.by = "Cluster_Enrichment") + 
              ggplot2::ggtitle("Cluster Enrichment"))
      
    }
    print("Returning Seurat object with enrichment results added to metadata.")
    return(seuratObj)
  } else if (returnSeuratObjectOrPlots == "Plots") {
    
    #guard for interactive use to prevent overwriting a Seurat object if following the example.
    if (Sys.info()[['sysname']] != "Linux") {
      if (utils::menu(c("Yes", "No"), title="Running this function using returnSeuratObjectOrPlots = 'Plots' will return a list of plots and a data frame, and NOT a Seurat Object.\nThis may overwrite your Seurat Object if you followed the example verbatim.\nDo you want this?") == 2) {
        stop("Please set returnSeuratObjectOrPlots to 'SeuratObject' to return a Seurat object with the GLMM enrichment results added to metadata.")
      }
    }
    
    #return a list of volcano plots and the model coefficients dataframe 
    return(list(volcano_plot = volcano_plot, 
                model_coefficients = results))
  }
}



#' @title ClusteredDotPlot
#' 
#' @description A function that generates a clustered dot plot with a heatmap of scaled expression.
#' @param seuratObj The Seurat object that holds the data.
#' @param features The features to plot.
#' @param groupFields The metadata column that is used for grouping.
#' @param assay The assay to plot.
#' @param scaling The scaling method for the heatmap. Options are "row", "column", or none.
#' @param layer The layer of the Seurat object that holds the relevant expression data.
#' @param dotSizeScaleFactor An integer used to scale the dot size. Larger values result in smaler max dot size
#' @param forceRescaling A boolean that determines if the Seurat object should be rescaled to include entries in the features vector if any are missing from the scale.data layer. This might be costly to perform locally.
#' @param inferDefaultArguments If TRUE, the function will infer the default arguments for the ComplexHeatmap::Heatmap function.
#' @param printInferredArguments Boolean to control optional printing of the arguments inferred by inferDefaultArguments.
#' @param numberColumns Boolean controlling the behavior of column titling by inferDefaultArguments. If TRUE, this will label each column's K means clusters with numeric titles.
#' @param numberRows Boolean controlling the behavior of row titling by inferDefaultArguments. If TRUE, this will label each row's K means clusters with numeric titles.
#' @param ... Additional arguments to pass to ComplexHeatmap::Heatmap
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' #set the seurat Idents for FindAllMarkers
#' Seurat::Idents(seuratObj) <- "ClusterNames_0.2"
#' 
#' markers <- Seurat::FindAllMarkers(seuratObj)
#' 
#' #filter markers to display the largest cluster identity markers according to average log fold change and differences in pct expression.
#' strong_markers <- markers[abs(markers$avg_log2FC) > 3 & abs(markers$pct.1 - markers$pct.2) > 0.25, "gene"]
#' 
#' dotPlot <- ClusteredDotPlot(seuratObj, features = strong_markers, groupFields = "ClusterNames_0.2", scaling = 'column')
#' print(dotPlot)
#' 
#' #Additionally, you can establish logical groupings based on a mixture of metadata fields
#'  #create individual classifications based on lineage markers
#'  CD3E_positive <- Seurat::WhichCells(seuratObj, expression = CD3E > 0)
#'  CD20_positive <- Seurat::WhichCells(seuratObj, expression = MS4A1 > 0)
#'  CD14_positive <- Seurat::WhichCells(seuratObj, expression = CD14 > 0)
#'  seuratObj$CD3E_positive <- ifelse(rownames(seuratObj@meta.data) %in% CD3E_positive, "CD3E_Pos", "CD3E_Neg")
#'  seuratObj$CD20_positive <- ifelse(rownames(seuratObj@meta.data) %in% CD20_positive, "CD20_Pos", "CD20_Neg")
#'  seuratObj$CD14_positive <- ifelse(rownames(seuratObj@meta.data) %in% CD14_positive, "CD14_Pos", "CD14_Neg")
#'  #roll lineage markers up into a single classification
#'  seuratObj$CellClassification <- paste0(seuratObj$CD3E_positive, "_", seuratObj$CD20_positive, "_" , seuratObj$CD14_positive)
#'  suppressWarnings(ClusteredDotPlot(seuratObj, 
#'                                  features = c("CD3E", "MS4A1", "CD14"),
#'                                  groupFields = c("CellClassification"),
#'                                  layer = 'data',
#'                                  scaling = "none", 
#'                                  height = grid::unit(10, "cm"), 
#'                                  width = grid::unit(4.5, "cm"), 
#'                                  show_row_dend = FALSE, 
#'                                  inferDefaultArguments = F
#'                                ))
#' 
#' }

ClusteredDotPlot <- function(seuratObj, 
                             features, 
                             groupFields = "ClusterNames_0.2", 
                             assay = "RNA", 
                             scaling = 'column', 
                             layer = 'data',
                             dotSizeScaleFactor = 50,
                             forceRescaling = FALSE, 
                             inferDefaultArguments = TRUE, 
                             printInferredArguments = FALSE,
                             numberColumns = T,
                             numberRows = T,
                             ...) {
  ## BEGIN ARGUMENT CHECKING
  #If you do some filtering upstream that removes all of the genes in your features vector, this doesn't error in an obvious way, so throw a specific error if you feed an empty vector into the features argument.
  if (length(features) == 0) {
    stop("The features argument is empty. Please specify a non-empty vector of features.")
  }
  #check that features are both valid and force them to be unique. 
  if (!all(features %in% rownames(Seurat::GetAssayData(seuratObj, layer = 'data')))) {
    warning(paste0('Features not found in Seurat object: ', paste0(features[!features %in% rownames(Seurat::GetAssayData(seuratObj, layer = layer))], collapse = ', ')))
  }
  if (any(duplicated(features))) {
    warning(paste0('Features supplied are not unique. Duplicates will be removed.'))
    features <- unique(features)
  }
  #check that groupFields are in the metadata
  if (!all(groupFields %in% colnames(seuratObj@meta.data))) { 
    stop(paste0('The following groupFields were not found in Seurat object metadata: ', paste0(groupFields[!groupFields %in% colnames(seuratObj@meta.data)], collapse = ', ')))
  }
  #check that scaling is supported
  if (!scaling %in% c('row', 'column', 'none')) {
    stop(paste0('Scaling method not supported: ', scaling, '. Please use one of: "row", "column", or "none".'))
  }
  #check assay
  if (!assay %in% Seurat::Assays(seuratObj)) {
    stop(paste0('Assay not found in Seurat object: ', assay))
  }
  #check if forceRescaling is a boolean
  if (!is.logical(forceRescaling)) {
    stop(paste0('forceRescaling: ', forceRescaling, ' is not a boolean. Please specify forceRescaling = TRUE or forceRescaling = FALSE. If TRUE, the Seurat object will be rescaled to include the features in the scale.data layer if any are missing.'))
  }
  # we could support this, but I think if someone wants to do this, they should be visualizing the covariance matrix instead. 
  if (layer == 'scale.data' & scaling != 'none') {
    stop("Further scaling of the scale.data layer is not supported. Please set scaling = 'none' or use the 'counts' or 'data' layers.")
  }
  #if we need to interact with the scale.data layer, we need to perform a bunch of checks to ensure features are or are not in the layer. Store this for now, since Seurat::GetAssayData() might be slow 
  if (layer == 'scale.data') {
    scaleDataFeatures <- rownames(GetAssayData(seuratObj, layer = 'scale.data'))
  }
  # ensure features (already sanitized above) are in the scale.data layer for automatic scaling in the Seurat::AverageExpression() function below if using the scale.data layer.
  # if not, offer to scale the data to include the features via the forceRescaling argument.
  # INFO: the feature selection for AverageSeurat reports the features already present in the given layer, so we need to recompute ScaleData to add them if they're missing.
  # see: https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/utilities.R#L1511C32-L1511C38
  # and: https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/utilities.R#L1478C5-L1478C17
  if (layer == 'scale.data' & !forceRescaling && !all(features %in% scaleDataFeatures)) {
    warning("Some features (reported below) are missing from the Seurat object's scale.data layer. These will be omitted from the final dotplot. To include them, please set the forceRescaling = TRUE and re-run the ClusteredDotPlot function.")
  } else if (layer == 'scale.data' & forceRescaling && !all(features %in% scaleDataFeatures)) {
    message(paste0("Features ", paste0(features[!features %in% scaleDataFeatures], collapse = ', '), " were not found in the scale.data layer. Rescaling the Seurat object to include these features."))
    seuratObj <- Seurat::ScaleData(seuratObj, features = features)
    #reset scaleDataFeatures if the data was rescaled
    scaleDataFeatures <- rownames(GetAssayData(seuratObj, layer = 'scale.data'))
  } 
  #AverageSeurat will fail with a cryptic matrix dimensionality error caused by accessing 0 or 1 features in the scale.data layer. Throw a more explicit error if that would fail.
  if (layer == 'scale.data' && length(intersect(features, scaleDataFeatures)) <= 2){
    stop("Less than two features would be present in the dot plot. Please set forceRescaling = TRUE to proceed with the scale.data layer, or use the 'data' or 'counts' and set the scaling argument to one of: 'column', 'row', or 'none'.")
  }
  #check booleans around argument inference and printing
  if (!is.logical(inferDefaultArguments)) {
    stop(paste0('inferDefaultArguments: ', inferDefaultArguments, ' is not a boolean. Please specify inferDefaultArguments = TRUE or inferDefaultArguments = FALSE. If TRUE, the function will infer the default arguments for the ComplexHeatmap::Heatmap function.'))
  }
  if (!is.logical(printInferredArguments)) {
    stop(paste0('printInferredArguments: ', printInferredArguments, ' is not a boolean. Please specify printInferredArguments = TRUE or printInferredArguments = FALSE. If TRUE, the function will print the arguments inferred by inferDefaultArguments.', collapse = ' '))
  }
  #Run a bunch of external checks on the ComplexHeatmap::Heatmap function to ensure that the arguments are valid.
  .CheckComplexHeatmapArguments(features = features, groupFields = groupFields, ...)
  
  ## END ARGUMENT CHECKING
  ## INFER HEATMAP DEFAULTS
  if (inferDefaultArguments) {
    inferred_heatmap_args <- .setComplexHeatmapDefaults(features = features, groupFields = groupFields, numberColumns = numberColumns, numberRows = numberRows, printInferredArguments = printInferredArguments, ...)
  } else {
    inferred_heatmap_args <- list(...)
  }
  ## END DEFAULT ARGUMENT INFERENCE
  ## BEGIN HEATMAP CONSTRUCTION
  
  #create averaged Seurat object for mean expression and subset features
  avgSeurat <- Seurat::AverageExpression(seuratObj, 
                                         features = features,
                                         group.by = c(groupFields),
                                         layer = layer, 
                                         return.seurat = T,
                                         assays = assay)
  
  #initialize the matrix for the heatmap
  mat <- Seurat::GetAssayData(avgSeurat, layer = layer)
  
  #scale, if requested (densifying in the process, if the matrix happens to be sparse)
  if (scaling %in% c('row', 'column')) {
    mat <- as.matrix(mat)  %>%
      Matrix::t() %>%
      pheatmap:::scale_mat(scale = scaling) 
  } else if (scaling == "none") {
    mat <- as.matrix(mat) %>%
      Matrix::t()
  }
  #establish the ordering of the expression heatmap, since one may want manual control over the ordering of the features using column_split.
  if (!is.null(inferred_heatmap_args[['column_split']])) {
    featureorder <- features
  } else {
    featureorder <- colnames(mat)
  }
  
  #harvest the percentage of cells expressing genes within the features vector from the Seurat::DotPlot output.
  #TODO: this works fine, but we have a version of this in the pseudobulking code. We could replace it if Seurat changes their dotplot. 
  plt <- Seurat::DotPlot(seuratObj, features = features, group.by = groupFields, assay = assay)
  pct <- plt$data %>%
    dplyr::select(pct.exp, id, features.plot) |>
    tidyr::pivot_wider(id_cols = features.plot, names_from = id, values_from = pct.exp) |>
    as.data.frame()
  row.names(pct) <- pct$features.plot
  pct <- pct |>
    dplyr::select(-features.plot) |>
    as.matrix() |>
    Matrix::t()
  #Establish symmetric color scaling based on the extremes in the heatmap
  col_RNA <- circlize::colorRamp2(quantile(c(-max(abs(mat)), 0, max(abs(mat)))),
                                  c("#0000FFFF", #blue
                                    "#7F53FDFF", #purple (pale)
                                    "gray90", #very light gray
                                    "#FF6948FF", #orangered/brick red
                                    "#FF0000FF"), #red
                                  space = "sRGB")
  #check if "g" characters were added to numeric rows in Seurat. 
  #This only happens if the original rownames are purely numeric. 
  
  if (all(paste0("g", rownames(pct)) %in% rownames(mat))) {
    rownames(pct) <- paste0("g", rownames(pct))
  } else if ( all(gsub("-", "_", rownames(mat)) %in% rownames(pct))) {
    rownames(mat) <- gsub("-", "_", rownames(mat))
  } else {
    message("There's a parsing error between Seurat and ComplexHeatmap. Try to avoid special characters in the name of your groupField.")
    message("Printing rownames to demonstrate the parsing issues:")
    print(unique(c(rownames(pct)[!rownames(pct) %in% rownames(mat)], rownames(mat)[!rownames(mat) %in% rownames(pct)])))
    stop("Parsing error, please see above.")
  }
  #enforce the order of the dotplot to match the heatmap (or the features vector argument if the heatmap uses column_split)
  pct <- pct[rownames(mat),featureorder]
  
  #Set the heatmap name according to scaling
  if (scaling == 'row') {
    legendName <- 'Scaled\nExpr. (Row)'
  } else if (scaling == 'column') {
    legendName <- 'Scaled\nExpr. (Column)'
  } else {
    #Average Seurat correctly averages the counts depending on the layer argument, but the transformation afterwards is also layer dependent and should be reported
    legendName <- paste0(c('Unscaled\nExpr.\n(layer = ', layer,')'), collapse = '')
  }
  
  #create the heatmap
  staticHeatmapArguments <- list(
    matrix = mat,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.circle(x = x, y = y, r = sqrt(pct[i,j])/dotSizeScaleFactor, default.units = "cm",
                        gp = grid::gpar(fill = col_RNA(mat[i, j])))
    },
    rect_gp = grid::gpar(type="none"),
    border_gp = grid::gpar(col = "black", lty = 1),
    name = legendName,
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    row_names_side = "left", 
    column_names_rot = 45
  )
  #merge the heatmap arguments - preferring the inferred arguments (which always includes args supplied in the ellipses)
  heatmapArgs <- .mergeComplexHeatmapArguments(inferred_args = inferred_heatmap_args,
                                               static_args = staticHeatmapArguments)
  
  suppressMessages(comp_heatmap <- do.call(ComplexHeatmap::Heatmap, heatmapArgs))
  print(comp_heatmap)
  return(comp_heatmap)
}

.CheckComplexHeatmapArguments <- function(features, groupFields, ...) { 
  passthrough_args <- list(...)
  #check that the arguments exist in ComplexHeatmap::Heatmap
  if (!all(names(passthrough_args) %in% names(formals(ComplexHeatmap::Heatmap)))) {
    stop(paste0('The following arguments are not supported by ComplexHeatmap::Heatmap: ', paste0(names(passthrough_args)[!names(passthrough_args) %in% names(formals(ComplexHeatmap::Heatmap))], collapse = ', ')))
  }
  #check column_km parameter 
  if (!is.null(passthrough_args[['column_km']]) && !(passthrough_args[['column_km']] %% 1 == 0)) {
    stop(paste0('K means column clustering parameter (column_km): ', passthrough_args[['column_km']], ' is not an integer. Please specify an integer value for column_km.'))
  } else if (!is.null(passthrough_args[['column_km']]) && passthrough_args[['column_km']] < 1) {
    stop(paste0('K means column clustering parameter (column_km): ', passthrough_args[['column_km']], ' is less than 1. Please specify an integer value greater than 1 for column_km.'))
  } 
  #check row_km parameter
  if (!is.null(passthrough_args[['row_km']]) && !(passthrough_args[['row_km']] %% 1 == 0) ) {
    stop(paste0('K means row clustering parameter (row_km): ', passthrough_args[['row_km']], ' is not an integer. Please specify an integer value for row_km.'))
  } else if (!is.null(passthrough_args[['row_km']]) && passthrough_args[['row_km']] < 1) {
    stop(paste0('K means row clustering parameter (row_km): ', passthrough_args[['row_km']], ' is less than 1. Please specify an integer value greater than 1 for row_km.'))
  } 
  #check that the length of the column_title matches the number of k-means clusters if specified
  if (!is.null(passthrough_args[['column_title']])) {
    if ((length(passthrough_args[['column_title']]) != 1) && !is.null(passthrough_args[['column_km']]) && !is.null(passthrough_args[['column_title']]) && length(passthrough_args[['column_title']]) != passthrough_args[['column_km']]) {
      stop(paste0('The length of column_title: ', length(passthrough_args[['column_title']]), ' does not match the number of k-means clusters (column_km): ', passthrough_args[['column_km']], '. Please specify a single title, NULL, or a vector of column titles that has elements equal to the value of column_km.'))
    } else if (length(passthrough_args[['column_title']]) == length(passthrough_args[['column_km']])) {
      warning('Please manually ensure that the order of the column_title matches the order of intended the k-means clusters in the heatmap.')
    }
  }
  #check that the length of the row_title matches the number of k-means clusters if specified
  if (!is.null(passthrough_args[['row_title']])) {
    if ((length(passthrough_args[['row_title']]) != 1) && !is.null(passthrough_args[['row_km']]) && !is.null(passthrough_args[['row_title']]) && length(passthrough_args[['row_title']]) != passthrough_args[['row_km']]) {
      stop(paste0('The length of row_title: ', length(passthrough_args[['row_title']]), ' does not match the number of k-means clusters (row_km): ', passthrough_args[['row_km']], '. Please specify a single title, NULL, or a vector of row titles that has elements equal to the value of row_km.'))
    } else if (length(passthrough_args[['row_title']]) == length(passthrough_args[['km']])) {
      warning('Please manually ensure that the order of the row_title matches the order of intended the k-means clusters in the heatmap.', immediate. = T)
    }
  }
  #check that the numbering defaults are boolean
  if (!is.null(passthrough_args[['numberColumns']]) && !is.logical(passthrough_args[['numberColumns']])) {
    stop(paste0('numberColumns: ', passthrough_args[['numberColumns']], ' is not a boolean. Please specify numberColumns = TRUE or numberColumns = FALSE. If TRUE and the column_title vector is not supplied, the column titles will be numbered.'))
  }
  if (!is.null(passthrough_args[['numberRows']]) && !is.logical(passthrough_args[['numberRows']])) {
    stop(paste0('numberRows: ', passthrough_args[['numberRows']], ' is not a boolean. Please specify numberRows = TRUE or numberRows = FALSE. If TRUE and the row_title vector is not supplied, the row titles will be numbered.'))
  }
  #check that height and width are valid unit variables
  if (!is.null(passthrough_args[['height']]) && !all(class(passthrough_args[['height']]) %in% c("simpleUnit", "unit", "unit_v2"))) { 
    stop(paste0('height: ', passthrough_args[['height']], ' is not a valid unit. Please specify a valid unit for the height argument, such as unit(7, "mm"). The default is the number of features multiplied by the function unit(25, "mm").'))
  }
  if (!is.null(passthrough_args[['width']]) && !all(class(passthrough_args[['width']]) %in% c("simpleUnit", "unit", "unit_v2"))) { 
    stop(paste0('width: ', passthrough_args[['width']], ' is not a valid unit. Please specify a valid unit for the width argument, such as unit(7, "mm"). The default is the number of groups (groupFields) multiplied by the function unit(7, "mm").'))
  }
  #check that title rotations are valid
  if (!is.null(passthrough_args[['row_title_rot']]) && !is.numeric(passthrough_args[['row_title_rot']])) {
    stop(paste0('row_title_rot: ', passthrough_args[['row_title_rot']], ' is not a numeric value. Please specify a numeric value for the angle to rotate the row titles.'))
  }
  if (!is.null(passthrough_args[['column_title_rot']]) && !is.numeric(passthrough_args[['column_title_rot']])) {
    stop(paste0('column_title_rot: ', passthrough_args[['column_title_rot']], ' is not a numeric value. Please specify a numeric value for the angle to rotate the column titles.'))
  }
  #check the booleans for showing the dendrograms
  if (!is.null(passthrough_args[['show_row_dend']]) && !is.logical(passthrough_args[['show_row_dend']])) {
    stop(paste0('show_row_dend: ', passthrough_args[['show_row_dend']], ' is not a boolean. Please specify show_row_dend = TRUE or show_row_dend = FALSE. If TRUE, the row dendrogram will be shown.'))
  }
  if (!is.null(passthrough_args[['show_column_dend']]) && !is.logical(passthrough_args[['show_column_dend']])) {
    stop(paste0('show_column_dend: ', passthrough_args[['show_column_dend']], ' is not a boolean. Please specify show_column_dend = TRUE or show_column_dend = FALSE. If TRUE, the column dendrogram will be shown.'))
  }
  #check that row_split and column_split are valid
  if (!is.null(passthrough_args[['row_split']]) && is.vector(passthrough_args[['row_split']])) {
    if (length(passthrough_args[['row_split']]) != length(groupFields)) {
      stop(paste0('row_split: ', paste0(passthrough_args[['row_split']], collapse = ', '), ' is not the same length as the groupFields vector. Please specify an integer value for row_split or a vector of length equal to the groupFields vector.'))
    } else if (!(length(passthrough_args[['row_split']]) %% 1 == 0) && !(passthrough_args[['row_split']] > 0)) {
      stop(paste0('row_split: ', paste0(passthrough_args[['row_split']], collapse = ", "), ' is not a positive integer or vector. Please specify either a positive integer, or a vector of length equal to the features vector for the row_split argument'))
    }
  } 
  if (!is.null(passthrough_args[['column_split']]) && is.vector(passthrough_args[['column_split']])) {
    if (length(passthrough_args[['column_split']]) != length(features)) {
      stop(paste0('column_split: ', paste0(passthrough_args[['column_split']], collapse= ', '), ' is not the same length as the features vector. Please specify either a positive integer value or a vector of length equal to the features vector for the column_split argument.'))
      
    } else if (!(length(passthrough_args[['column_split']]) %% 1 == 0) && !(passthrough_args[['column_split']] > 0)) {
      stop(paste0('column_split: ', paste0(passthrough_args[['column_split']], collapse= ', '), ' is not a positive integer or vector. Please specify an integer value for column_split.'))
    }
  }
} 

.setComplexHeatmapDefaults <- function(features, groupFields, numberColumns = T, numberRows = T, printInferredArguments, ...) {
  passthrough_args <- list(...)
  
  #if the user doesn't specify row_km or column_km, infer that no clustering is desired
  if (is.null(passthrough_args[['row_km']])) {
    #set a flag for downstream warnings to indicate that we inferred this row_km
    inferred_row_km <- TRUE
    passthrough_args[['row_km']] <- 1
  }
  if (is.null(passthrough_args[['column_km']])) {
    #set a flag for downstream warnings to indicate that we inferred this column_km
    inferred_column_km <- TRUE
    passthrough_args[['column_km']] <- 1
  }
  #number columns if column_km isn't null, but the user didn't specify column_title.
  #however, if the user specifies column_split, we assume that they want to manually define the feature grouping, and we shouldn't infer the column titles.
  if (is.null(passthrough_args[['column_title']]) && !is.null(passthrough_args[['column_km']]) && numberColumns && is.null(passthrough_args[['column_split']])) {
    passthrough_args[['column_title']] <- 1:passthrough_args[['column_km']]
  }
  #number rows if row_km isn't null, but the user didn't specify row_title.
  #however, if the user specifies row_split, we assume that they want to manually define the sample-level grouping, and we shouldn't infer the row titles.
  if (is.null(passthrough_args[['row_title']]) && !is.null(passthrough_args[['row_km']]) && numberRows && is.null(passthrough_args[['row_split']])) {
    passthrough_args[['row_title']] <- 1:passthrough_args[['row_km']]
  }
  #set defaults for height and width based on the number of features and groups if left unspecified.
  if (is.null(passthrough_args[['height']])) {
    passthrough_args[['height']] <- length(groupFields)*unit(25, "mm")
  }
  if (is.null(passthrough_args[['width']])) {
    passthrough_args[['width']] <- length(features)*unit(7, "mm")
  }
  #set show_row_dend to TRUE if row_km is greater than 1 and show_row_dend is NULL
  if (is.null(passthrough_args[['show_row_dend']]) && passthrough_args[['row_km']] > 1) {
    passthrough_args[['show_row_dend']] <- TRUE
  }
  #set show_column_dend to TRUE if column_km is greater than 1 and show_column_dend is NULL
  if (is.null(passthrough_args[['show_column_dend']]) && passthrough_args[['column_km']] > 1) {
    passthrough_args[['show_column_dend']] <- TRUE
  }
  #if column/row split is defined and also column/row_km is defined, keep the split but warn the user. 
  if (!is.null(passthrough_args[['row_split']]) && !is.null(passthrough_args[['row_km']])) {
    #it's possible that this warning can trip because we inferred the row_km argument - if that's the case, suppress it.
    if (!inferred_row_km) {
      warning("Both row_split and row_km are defined. The row_km argument will be ignored, and the row_split argument will be used to group features into subpanels regardless of clustering.")
    }
    passthrough_args[['row_km']] <- NULL
  }
  if (!is.null(passthrough_args[['column_split']]) && !is.null(passthrough_args[['column_km']])) {
    #it's possible that this warning can trip because we inferred the column_km argument - if that's the case, suppress it.
    if (!inferred_column_km) {
      warning("Both column_split and column_km are defined. The column_km argument will be ignored, and the column_split argument will be used to group groupFields into subpanels regardless of clustering.")
    }
    passthrough_args[['column_km']] <- NULL
  }
  if (printInferredArguments) {
    print(passthrough_args)
  }
  return(passthrough_args)
}

.mergeComplexHeatmapArguments <- function(inferred_args, static_args) {
  intersecting_args <- intersect(names(inferred_args), names(static_args))
  #if there are intersecting arguments, prefer the inferred arguments
  if (length(intersecting_args) > 0) {
    heatmapArgs <- c(inferred_args, static_args[!names(static_args) %in% intersecting_args])
  } else {
    heatmapArgs <- c(inferred_args, static_args)
  }
  return(heatmapArgs)
}
