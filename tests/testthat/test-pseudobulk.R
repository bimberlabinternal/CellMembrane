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

test_that("nCount RNA stratification skips when data is too small for KL divergence", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  tiny <- seuratObj[, 1:2]
  tiny$kl_test_cluster <- 0L
  pseudo_tiny <- testthat::expect_no_error(
    PseudobulkSeurat(
      tiny,
      groupFields = c('ClusterNames_0.2'),
      nCountRnaStratification = TRUE,
      stratificationGroupingFields = 'kl_test_cluster'
    )
  )
  testthat::expect_false(is.null(pseudo_tiny))
  testthat::expect_false('nCount_RNA_Stratification' %in% colnames(pseudo_tiny@meta.data))
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
  CellMembrane::SetSeed(CellMembrane::GetSeed())
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
                       filterGenes = TRUE, 
                       assayName = "RNA")
  #15 total contrasts
  testthat::expect_equal(length(DE_results), expected = 15)
  #9008 "DEGs" in the first contrast (note overly permissive DEG thresholds)
  testthat::expect_equal(nrow(DE_results$`1`), expected = 900)
  #test that pct.1 and pct.2 are present in the DE results
  testthat::expect_true(all(c("pct.1", "pct.2") %in%  colnames(DE_results$`1`)))
  
  barPlot <- PseudobulkingBarPlot(filteredContrasts = DE_results, 
                       metadataFilterList = NULL, 
                       logFC_threshold = 0,
                       FDR_threshold = 0.5
                       )
  #test that PseudobulkingBarPlot yields a list with a barPlot element
  testthat::expect_true("barPlot" %in% names(barPlot))
  print('Class of barPlot$barPlot')
  print(class(barPlot$barPlot))
  testthat::expect_contains(class(barPlot$barPlot), c("gg", "ggplot"))
  
  genes <- rownames(pbulk)[1:10]
  
  heatmap_list <- PseudobulkingDEHeatmap(seuratObj = pbulk, 
                                         geneSpace = genes, 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "control", 
                                         sampleIdCol = 'subject', 
                                         filterGenes = FALSE
  )
  
  #test that PseudobulkingDEHeatmap yields a list with a heatmap and matrix element
  testthat::expect_true(all(c("heatmap", "matrix") %in% names(heatmap_list)))
  #test that the barPlot element is a ComplexHeatmap object (S4)
  testthat::expect_true(typeof(heatmap_list$heatmap) == "S4")
  #test that the heatmap matrix has 3 columns
  testthat::expect_true(ncol(heatmap_list$matrix) == 3)
  #test that pseudobulk heatmap subsetting works. 
  testthat::expect_no_error(suppressWarnings(PseudobulkingDEHeatmap(seuratObj = pbulk,
                                         geneSpace = genes[genes!="PRF1"], 
                                         contrastField = "vaccine_cohort", 
                                         negativeContrastValue = "control", 
                                         sampleIdCol = 'subject', 
                                         subsetExpression = "vaccine_cohort %in% c('unvax', 'vaccineOne')"
  )))
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

test_that("Technical covariates work in the DE pipeline", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  #add fabricated study metadata including a technical covariate
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  #derive PlateId from subject so it varies within vaccine groups (avoids confounding with vaccine_cohort)
  seuratObj@meta.data[,"PlateId"] <- ifelse(seuratObj@meta.data[,"subject"] %in% c(1, 2), "PlateA", "PlateB")

  #pseudobulk including the covariate in groupFields so it appears in sample-level metadata
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint", "subject", "PlateId"))

  #test DesignModelMatrix with technicalCovariates
  design_with_cov <- DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject", technicalCovariates = c("PlateId"))

  #the design matrix should contain group columns plus covariate columns
  group_cols <- levels(factor(interaction(pbulk$vaccine_cohort, pbulk$timepoint, sep = "_")))
  testthat::expect_true(any(grepl("PlateId", colnames(design_with_cov))))
  testthat::expect_equal(attr(design_with_cov, "technicalCovariates"), c("PlateId"))

  #test that the baseline design (no covariates) has fewer columns
  design_no_cov <- DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject")
  testthat::expect_true(ncol(design_with_cov) > ncol(design_no_cov))
  testthat::expect_null(attr(design_no_cov, "technicalCovariates"))

  #test that GLM fitting works with the covariate design
  fit <- PerformGlmFit(pbulk, design = design_with_cov, test.use = "QLF")
  testthat::expect_true("coefficients" %in% names(fit))

  #test erroring when the covariate is not in metadata
  testthat::expect_error(DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort"), sampleIdCol = "subject", technicalCovariates = c("NonExistentField")))
  #test validation: covariate overlaps with contrast_columns
  testthat::expect_error(DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort"), sampleIdCol = "subject", technicalCovariates = c("vaccine_cohort")))

  #test RunFilteredContrasts propagates covariates
  logic_list <- list(list("vaccine_cohort", "xor", "control"))
  filtered_contrasts <- FilterPseudobulkContrasts(logicList = logic_list,
                                                  design = design_with_cov,
                                                  useRequireIdenticalLogic = T,
                                                  requireIdenticalFields = c("timepoint"),
                                                  filteredContrastsOutputFile = tempfile())
  DE_results_cov <- RunFilteredContrasts(seuratObj = pbulk,
                       filteredContrastsDataframe = filtered_contrasts,
                       design = design_with_cov,
                       test.use = "QLF",
                       logFC_threshold = 0,
                       FDR_threshold = 0.5,
                       filterGenes = TRUE,
                       assayName = "RNA")
  testthat::expect_true(length(DE_results_cov) > 0)
  testthat::expect_true(all(c("logFC", "FDR", "gene") %in% colnames(DE_results_cov[[1]])))
})

test_that("Rank deficiency is detected and reported when covariates are confounded", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  #deliberately confound PlateId with vaccine_cohort: period-2 rep perfectly nests inside period-4 rep
  seuratObj@meta.data[,"PlateId"] <- base::rep(c("PlateA", "PlateB"), length.out = length(colnames(seuratObj)))

  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint", "subject", "PlateId"))

  #DesignModelMatrix should raise an error about rank deficiency
  testthat::expect_error(
    DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject", technicalCovariates = c("PlateId")),
    regexp = "rank-deficient"
  )

  #no error when covariates are not confounded (PlateId derived from randomly assigned subject)
  seuratObj@meta.data[,"PlateId"] <- ifelse(seuratObj@meta.data[,"subject"] %in% c(1, 2), "PlateA", "PlateB")
  pbulk_ok <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint", "subject", "PlateId"))
  testthat::expect_no_error(
    DesignModelMatrix(pbulk_ok, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject", technicalCovariates = c("PlateId"))
  )
})

test_that("Feature Selection by GLM works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  testthat::expect_equal(ncol(seuratObj), expected = 1557) #check that test seuratObj doesn't change
  #add fabricated study metadata
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::rep(c(1:10), length.out = length(colnames(seuratObj)))
  #set up pseudobulking
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint","subject"))
  classification_results <- suppressWarnings(FitRegularizedClassificationGlm(pbulk,
                                                          metadataVariableForClassification = "vaccine_cohort", 
                                                          returnModelAndSplits = T, 
                                                          rescale = F
                                                          ))
  #ensure a glmnet specific parameter exists in the fitted model, showing that it trained successfully
  expect_true("lambda.min" %in% names(classification_results$model))
})
