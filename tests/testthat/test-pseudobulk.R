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
  pbmc_small@meta.data[,"random_cohort"] <- rep(c(1,2), size = length(colnames(pbmc_small)))
  
  pbmc_pbulk <- PseudobulkSeurat(pbmc_small, groupFields = c("RNA_snn_res.0.8", "letter.idents", "groups", "random_cohort"))
  design <- DesignModelMatrix(pbmc_pbulk, contrast_columns = c("RNA_snn_res.0.8", "letter.idents","random_cohort"), sampleIdCol = "groups")
  fit <- PerformGlmFit(pbmc_pbulk, design = design, test.use = "QLF")
  pairwise_de_results <- RunPairwiseContrasts(fit, test.use = "QLF", logFC_threshold = 1)
  testthat::expect_equal(sum(pairwise_de_results$`X0_A_1-X0_A_2`$differential_expression$table$FDR < 0.05), expected = 0) #this contrast should have no significantly differentially expressed genes
  testthat::expect_equal(sum(pairwise_de_results$`X0_A_1-X1_B_1`$differential_expression$table$FDR < 0.05), expected = 44) #this contrast should have 44 genes that pass the FDR threshold.
  
  bar_plot <- CreateStudyWideBarPlot(pairwise_de_results = pairwise_de_results, logFC_threshold = 1)
  testthat::expect_equal(typeof(bar_plot), expected = "list") #ensure that the ggplot was created.
})
