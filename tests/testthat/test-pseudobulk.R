context("scRNAseq")

test_that("Pseudobulk works", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
  
  pseudo <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.2'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.2)), nrow(pseudo@meta.data))
  expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo@assays$RNA))
  
  pseudo2 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), assays = c('RNA'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.4)), nrow(pseudo2@meta.data))
  expect_equal(nrow(seuratObj@assays$RNA), nrow(pseudo2@assays$RNA))
  
  pseudo3 <- PseudobulkSeurat(seuratObj, groupFields = c('ClusterNames_0.4'), assays = c('RNA'), additionalFieldsToAggregate = c('G2M.Score', 'p.mito'))
  expect_equal(length(unique(seuratObj$ClusterNames_0.4)), nrow(pseudo3@meta.data))
  expect_equal(max(pseudo3$G2M.Score_mean, na.rm = T), -0.007676878)
  expect_equal(min(pseudo3$G2M.Score_mean, na.rm = T), -0.02633076)
  
  expect_equal(max(pseudo3$p.mito_mean, na.rm = T), 0.06488353)
  expect_equal(min(pseudo3$p.mito_mean, na.rm = T), 0.04050757)
})

test_that("Pseudobulk-based differential expression works", {
  testthat::expect_equal(length(colnames(pbmc_small)), expected = 80) #weak insurance that pbmc_small doesn't change.
  pbmc_small@meta.data[,"random_cohort"] <- base::rep(c(1,2), length.out = length(colnames(pbmc_small)))
  
  pbmc_pbulk <- PseudobulkSeurat(pbmc_small, groupFields = c("RNA_snn_res.0.8", "letter.idents", "groups", "random_cohort"))
  design <- DesignModelMatrix(pbmc_pbulk, contrast_columns = c("RNA_snn_res.0.8", "letter.idents","random_cohort"), sampleIdCol = "groups")
  fit <- PerformGlmFit(pbmc_pbulk, design = design, test.use = "QLF")
  pairwise_de_results <- RunPairwiseContrasts(fit, test.use = "QLF", logFC_threshold = 1)
  testthat::expect_equal(sum(pairwise_de_results$`X0_A_1-X0_A_2`$differential_expression$table$FDR < 0.05), expected = 0) #this contrast should have no significantly differentially expressed genes
  testthat::expect_equal(sum(pairwise_de_results$`X0_A_1-X1_B_1`$differential_expression$table$FDR < 0.05), expected = 44) #this contrast should have 44 genes that pass the FDR threshold.
  
  bar_plot <- CreateStudyWideBarPlot(pairwise_de_results = pairwise_de_results, logFC_threshold = 1)
  testthat::expect_equal(typeof(bar_plot), expected = "list") #ensure that the ggplot was created.
})

test_that("Logic gate study design works", {
  seuratObj <- readRDS('../testdata/seuratOutput.rds')
  testthat::expect_equal(ncol(seuratObj), expected = 1557) #check that test seuratObj doesn't change
  #add fabricated study metadata
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy", "day4"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- c(base::rep(1, length.out = 1000), base::rep(2, length.out = 557))
  #set up pseudobulking
  pbulk <- PseudobulkSeurat(seuratObj, groupFields = c("vaccine_cohort", "timepoint","subject"))
  design <- DesignModelMatrix(pbulk, contrast_columns = c("vaccine_cohort", "timepoint"), sampleIdCol = "subject")
  #create logical_dataframe/study design
  fields_to_check <- c("vaccine_cohort", "timepoint")
  field_logic <- c("xor", "any") #check for equality. any can be used for an OR gate, all can be used as an AND gate. 
  criteria <- c("control", "necropsy")
  logical_dataframe <- data.frame(fields_to_check = fields_to_check, field_logic = field_logic, criteria = criteria)
  filtered_contrasts <- FilterPseudobulkContrasts(logical_dataframe = logical_dataframe, 
                                                  design = design, 
                                                  contrast_columns = c("vaccine_cohort", "timepoint"), 
                                                  use_require_identical_logic = F, 
                                                  require_identical_fields = NULL, 
                                                  filtered_contrasts_output_file = tempfile())
  #test filtering worked as expected
  testthat::expect_equal(nrow(filtered_contrasts), expected =  15)
  #test that xor gate worked
  testthat::expect_true(all(apply(filtered_contrasts, 
                                  MARGIN = 1,
                                  FUN = function(x){xor(grepl("control",x[1]), grepl("control",x[2]))})))
  #test that OR gate worked (i.e. necropsy appears in at least one of the two sides of the contrast)
  testthat::expect_false(any(apply(filtered_contrasts, 
                                   MARGIN = 1,
                                   FUN = function(x){nor(grepl("necropsy",x[1]), grepl("necropsy",x[2]))})))
  #Test require_identical_field logic
  fields_to_check <- c("vaccine_cohort")
  field_logic <- c("xor") #check for equality. any can be used for an OR gate, all can be used as an AND gate. 
  criteria <- c("control")
  logical_dataframe <- data.frame(fields_to_check = fields_to_check, field_logic = field_logic, criteria = criteria)
  require_identical_fields <- c("timepoint")
  filtered_contrasts_require_identical <- FilterPseudobulkContrasts(logical_dataframe = logical_dataframe, 
                                                            design = design, 
                                                            contrast_columns = c("vaccine_cohort", "timepoint"), 
                                                            use_require_identical_logic = T, 
                                                            require_identical_fields = c("timepoint"), 
                                                            filtered_contrasts_output_file = tempfile())
  #test filtering worked as expected
  testthat::expect_equal(nrow(filtered_contrasts_require_identical), expected =  9)
  #test invariance
  positive_side_timepoints <- unlist(strsplit(filtered_contrasts_require_identical[,1], split = '_'))
  negative_side_timepoints <- unlist(strsplit(filtered_contrasts_require_identical[,2], split = '_'))
  #test that every second value in the contrasts (timepoint) are always equal
  testthat::expect_true(all(positive_side_timepoints[1:length(positive_side_timepoints)%% 2 == 0] == negative_side_timepoints[1:length(negative_side_timepoints)%% 2 == 0]))
})
