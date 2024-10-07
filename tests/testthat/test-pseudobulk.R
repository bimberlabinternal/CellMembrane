context("scRNAseq")

test_that("Pseudobulk works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  pseudo <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.2'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.2)), nrow(pseudo@meta.data))
  expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo@assays$RNA))
  expect_equal(64488, max(pseudo@assays$RNA$counts))
  expect_equal(21.9, mean(as.matrix(pseudo@assays$RNA$counts), na.rm = TRUE), tolerance = 0.0001)
  
  ld <- SeuratObject::LayerData(pseudo, layer = 'pct.expression', assay = 'RNA')
  expect_equal(1, max(ld))
  expect_equal(0.0261, mean(as.matrix(ld), na.rm = TRUE), tolerance = 0.001)
  
  pseudo2 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), assayToAggregate = c('RNA'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.4)), nrow(pseudo2@meta.data))
  expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo2@assays$RNA))

  pseudo3 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), assayToAggregate = c('RNA'), additionalFieldsToAggregate = c('G2M.Score', 'p.mito'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.4)), nrow(pseudo3@meta.data))
  expect_equal(max(pseudo3$G2M.Score_mean, na.rm = T), -0.007676878)
  expect_equal(min(pseudo3$G2M.Score_mean, na.rm = T), -0.02633076)
  
  expect_equal(max(pseudo3$p.mito_mean, na.rm = T), 0.06488353)
  expect_equal(min(pseudo3$p.mito_mean, na.rm = T), 0.04050757)
  
  pseudo4 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), nCountRnaStratification = T)
  testthat::expect_false(all(pseudo4$nCount_RNA_Stratification == "NormalRNACounts"))
  testthat::expect_true(pseudo4@meta.data[pseudo4$ClusterNames_0.4 == 5, "nCount_RNA_Stratification"] == "AbnormalRNACounts")
  testthat::expect_warning(PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), nCountRnaStratification = T))

})

test_that("Pseudobulk-based differential expression works", {
  testthat::expect_equal(length(colnames(pbmc_small)), expected = 80) #weak insurance that pbmc_small doesn't change.
  pbmc_small@meta.data[,"random_cohort"] <- base::rep(c(1,2), length.out = length(colnames(pbmc_small)))
  
  pbmc_pbulk <- PseudobulkSeurat(pbmc_small, groupFields = c("RNA_snn_res.0.8", "letter.idents", "groups", "random_cohort"))
  design <- DesignModelMatrix(pbmc_pbulk, contrast_columns = c("RNA_snn_res.0.8", "letter.idents","random_cohort"), sampleIdCol = "groups")
  fit <- PerformGlmFit(pbmc_pbulk, design = design, test.use = "QLF")
  testthat::expect_true("coefficients" %in% names(fit))
})

test_that("Logic gate study design works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  testthat::expect_equal(ncol(seuratObj), expected = 1557) #check that test seuratObj doesn't change
  #add fabricated study metadata
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day_4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  #set up pseudobulking
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint","subject"))
  design <- DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject")
  #create logical_dataframe/study design
  logic_list <- list( list("vaccine_cohort", "xor", "control"), list("timepoint", "any", "necropsy"))
  filtered_contrasts <- FilterPseudobulkContrasts(logicList = logic_list, 
                                                  design = design, 
                                                  useRequireIdenticalLogic = F, 
                                                  requireIdenticalFields = NULL, 
                                                  filteredContrastsOutputFile = tempfile())
  #test filtering worked as expected
  testthat::expect_equal(nrow(filtered_contrasts), expected =  15)
  #test that xor gate on vaccine_cohort field worked
  testthat::expect_true(all(apply(filtered_contrasts, 
                                  MARGIN = 1,
                                  FUN = function(x){xor(grepl("control",x["positive_contrast_vaccine_cohort"]), 
                                                        grepl("control",x["negative_contrast_vaccine_cohort"]))})))
  #test that OR gate worked (i.e. necropsy appears in at least one of the two sides of the contrast)
  testthat::expect_false(any(apply(filtered_contrasts, 
                                   MARGIN = 1,
                                   FUN = function(x){nor(grepl("necropsy",x["positive_contrast_timepoint"]), 
                                                         grepl("necropsy",x["negative_contrast_timepoint"]))})))
  #Test require_identical_field logic
  logic_list <- list( list("vaccine_cohort", "xor", "control"))
  require_identical_fields <- c("timepoint")
  filtered_contrasts_require_identical <- FilterPseudobulkContrasts(logicList = logic_list, 
                                                            design = design, 
                                                            useRequireIdenticalLogic = T, 
                                                            requireIdenticalFields = c("timepoint"), 
                                                            filteredContrastsOutputFile = tempfile())
  #test filtering worked as expected
  testthat::expect_equal(nrow(filtered_contrasts_require_identical), expected =  9)
  #test that timepoints are always equal
  testthat::expect_true(all(filtered_contrasts_require_identical[,"positive_contrast_timepoint"] == filtered_contrasts_require_identical[,"negative_contrast_timepoint"]))
  
  #Test that running differential expression on filtered constrasts works. 
  DE_results <- RunFilteredContrasts(seuratObj = pbulk, 
                       filteredContrastsDataframe = filtered_contrasts, 
                       design = design,
                       test.use = "QLF", 
                       logFC_threshold = 0,
                       FDR_threshold = 0.5,
                       minCountsPerGene = 1, 
                       assayName = "RNA")
  #15 total contrasts
  testthat::expect_equal(length(DE_results), expected = 15)
  #9008 "DEGs" in the first contrast (note overly permissive DEG thresholds)
  testthat::expect_equal(nrow(DE_results$`1`), expected = 9008)
  #test that pct.1 and pct.2 are present in the DE results
  testthat::expect_true(all(c("pct.1", "pct.2") %in%  colnames(DE_results$`1`)))
  
  barPlot <- PseudobulkingBarPlot(filteredContrasts = DE_results, 
                       metadataFilterList = NULL
                       )
  #test that PseudobulkingBarPlot yields a list with a barPlot element
  testthat::expect_true("barPlot" %in% names(barPlot))
  #test that the barPlot element is a ggplot object (list)
  testthat::expect_true("list" == typeof(barPlot$barPlot))
  
  genes <- rownames(pbulk)[1:10]
  
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "control", 
                                         sampleIdCol = 'subject'
  )
  
  #test that PseudobulkingDEHeatmap yields a list with a heatmap and matrix element
  testthat::expect_true(all(c("heatmap", "matrix") %in% names(heatmap_list)))
  #test that the barPlot element is a ComplexHeatmap object (S4)
  testthat::expect_true(typeof(heatmap_list$heatmap) == "S4")
  #test that the heatmap matrix has 3 columns
  testthat::expect_true(ncol(heatmap_list$matrix) == 3)
  #test that pseudobulk heatmap subsetting works. 
  testthat::expect_no_error(PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes[genes!="PRF1"], 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "control", 
                                         sampleIdCol = 'subject', 
                                         subsetExpression = "vaccine_cohort %in% c('unvax', 'vaccineOne')"
  ))
})

test_that("Non-alphanumeric characters do not break pseudobulking pipeline", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  testthat::expect_equal(ncol(seuratObj), expected = 1557) #check that test seuratObj doesn't change
  #add fabricated study metadata with odd characters in the metadata fields
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("cont-rol", "vaccine One", "vaccine_Two", "un$vax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day_4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint","subject"))
  genes <- rownames(pbulk)[1:10]
  
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "cont-rol", 
                                         sampleIdCol = 'subject')
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "vaccine One", 
                                         sampleIdCol = 'subject')
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "un$vax", 
                                         sampleIdCol = 'subject')
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "vaccine_Two", 
                                         sampleIdCol = 'subject')
  })

test_that("Feature Selection by GLM works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  testthat::expect_equal(ncol(seuratObj), expected = 1557) #check that test seuratObj doesn't change
  #add fabricated study metadata
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- c(base::rep(1, length.out = 1000), base::rep(2, length.out = 557))
  #set up pseudobulking
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint","subject"))
  print('DEBUG!!')
  table(colnames(pbulk))
  classification_results <- suppressWarnings(FitRegularizedClassificationGlm(pbulk, 
                                                          metadataVariableForClassification = "vaccine_cohort", 
                                                          returnModelAndSplits = T, 
                                                          rescale = F
                                                          ))
  #ensure a glmnet specific parameter exists in the fitted model, showing that it trained successfully
  expect_true("lambda.min" %in% names(classification_results$model))
})
