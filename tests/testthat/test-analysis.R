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
  #test that " - " breaks the function
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("base - line", "necropsy"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  testthat::expect_error(CalculateClusterEnrichment(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "timepoint",
                                          subjectField = "subject",
                                          paired = "infer"))
})

test_that("ClusteredDotPlot works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  #general testing of ClusteredDotPlot
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne"), length.out = length(colnames(seuratObj)))
  
  dotplot <- ClusteredDotPlot(seuratObj, 
                              features = c("CD3E", "CD8A"),
                              groupFields = c("vaccine_cohort"),
                              layer = 'data',
                              scaling = "column")
  #test that the plot was created and the name is correct
  testthat::expect_true(dotplot@name == "Scaled\nExpr. (Column)")
  #mixed scaling should fail
  testthat::expect_error(ClusteredDotPlot(seuratObj, 
                                features = c("CD3E", "CD8A"),
                                groupFields = c("vaccine_cohort"),
                                layer = 'scale.data',
                                scaling = "row"))
  #this should fail due to CD8A being the only variable present in the scale.data layer by default.
  testthat::expect_error(suppressWarnings(ClusteredDotPlot(seuratObj, 
                                                 features = c("CD3E", "CD8A", "MKI67"),
                                                 groupFields = c("vaccine_cohort"),
                                                 layer = 'scale.data',
                                                 scaling = "none")))
  #however, rescaling the matrix should fix the issues
  testthat::expect_no_error(suppressWarnings(ClusteredDotPlot(seuratObj, 
                                                    features = c("CD3E", "CD8A", "MKI67"),
                                                    groupFields = c("vaccine_cohort"),
                                                    layer = 'scale.data',
                                                    scaling = "none", 
                                                    forceRescaling = T)))
})
