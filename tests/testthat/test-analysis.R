context("scRNAseq")

test_that("Cluster enrichment works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  #test that the Kruskal-Wallis rank sum test runs
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  seuratObj <- CalculateClusterEnrichment(seuratObj, 
                            clusterField = "ClusterNames_0.4", 
                            treatmentField = "vaccine_cohort",
                            subjectField = "subject",
                            paired = "infer")
  #test that the Kruskal-Wallis rank sum test runs
  testthat::expect_true("Cluster_p_adj" %in% colnames(seuratObj@meta.data))  
  testthat::expect_true(all(seuratObj$Cluster_p_adj <= 1))
  
  #erase results
  seuratObj$Cluster_p_adj <- NULL
  seuratObj$Cluster_pValue <- NULL
  
  #test that the Wilcoxon rank sum test runs 
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  seuratObj <- CalculateClusterEnrichment(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "vaccine_cohort",
                                          subjectField = "subject",
                                          paired = "infer")
  testthat::expect_true("Cluster_p_adj" %in% colnames(seuratObj@meta.data))  
  testthat::expect_true(all(seuratObj$Cluster_p_adj <= 1))
  
  #erase results
  seuratObj$Cluster_p_adj <- NULL
  seuratObj$Cluster_pValue <- NULL
  
  #test that the Wilcoxon signed rank test runs 
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  seuratObj <- CalculateClusterEnrichment(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "timepoint",
                                          subjectField = "subject",
                                          paired = "infer")
  testthat::expect_true("Cluster_p_adj" %in% colnames(seuratObj@meta.data))  
  testthat::expect_true(all(seuratObj$Cluster_p_adj <= 1))
})
