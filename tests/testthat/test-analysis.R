context("scRNAseq")

test_that("Cluster enrichment works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  #test that the Kruskal-Wallis rank sum test runs
  seuratObj@meta.data[,"vaccine_cohort"] <- base::rep(c("control", "vaccineOne", "vaccineTwo", "unvax"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  seuratObj <- CalculateClusterEnrichmentOmnibus(seuratObj, 
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
  
  seuratObj <- CalculateClusterEnrichmentOmnibus(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "vaccine_cohort",
                                          subjectField = "subject",
                                          paired = "infer")
  testthat::expect_true("Cluster_p_adj" %in% colnames(seuratObj@meta.data))  
  testthat::expect_true(all(seuratObj$Cluster_p_adj <= 1))

  #test that the Wilcoxon signed rank test runs 
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("baseline", "necropsy"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  seuratObj <- CalculateClusterEnrichmentOmnibus(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "timepoint",
                                          subjectField = "subject",
                                          paired = "infer")
  testthat::expect_true("Cluster_p_adj" %in% colnames(seuratObj@meta.data))  
  testthat::expect_true(all(seuratObj$Cluster_p_adj <= 1))
  #test that " - " breaks the function
  seuratObj@meta.data[,"timepoint"] <- base::rep(c("base - line", "necropsy"), length.out = length(colnames(seuratObj)))
  seuratObj@meta.data[,"subject"] <- base::sample(c(1,2,3,4), size = 1557, replace = T)
  
  testthat::expect_error(CalculateClusterEnrichmentOmnibus(seuratObj, 
                                          clusterField = "ClusterNames_0.4", 
                                          treatmentField = "timepoint",
                                          subjectField = "subject",
                                          paired = "infer"))
  
  #test GLMM-based cluster enrichment
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  seuratObj$cDNA_ID <- rep(1:4, ncol(seuratObj)) 
  seuratObj$Vaccine <- rep(c("Vaccine1", "Vaccine2"), each = ncol(seuratObj)/2) 
  seuratObj$SubjectId <- rep(1:4, each = ncol(seuratObj)/4) 
  
  #induce a bias in cluster 1
  seuratObj$Vaccine <- ifelse(seuratObj$ClusterNames_0.2 == 1, 
                              sample(c("Vaccine1", "Vaccine2"), 
                                     size = sum(seuratObj$ClusterNames_0.2 == 1),
                                     replace = TRUE, 
                                     prob = c(0.9, 0.1)), 
                              seuratObj$Vaccine)
  
  
  testthat::expect_no_error(seuratObj <- CalculateClusterEnrichmentPairwise(seuratObj,
                                              subjectField = 'SubjectId',
                                              clusterField = 'ClusterNames_0.2',
                                              biologicalReplicateGroupingVariables = c("cDNA_ID"),
                                              treatmentField = "Vaccine",
                                              referenceValue = "Vaccine1",
                                              pValueCutoff = 0.05,
                                              showPlots = FALSE, 
                                              returnSeuratObjectOrPlots = "SeuratObject", 
                                              includeDepletions = TRUE))
  #test that the GLMM enrichment ran
  testthat::expect_true("Depleted: Vaccine2:4" %in% seuratObj$Cluster_Enrichment)
  testthat::expect_true("Depleted: Vaccine2:1" %in% seuratObj$Cluster_Enrichment)
  testthat::expect_false("Enriched: Vaccine2:3" %in% seuratObj$Cluster_Enrichment)
  testthat::expect_true(sum(is.na(seuratObj$Estimate)) == 1471)
  testthat::expect_equal(mean(seuratObj$Estimate, na.rm = TRUE), expected = -1, tolerance = 1)
  
  
  testthat::expect_no_error(plots <- CalculateClusterEnrichmentPairwise(seuratObj,
                                              subjectField = 'SubjectId',
                                              clusterField = 'ClusterNames_0.2',
                                              biologicalReplicateGroupingVariables = c("cDNA_ID"),
                                              treatmentField = "Vaccine",
                                              referenceValue = "Vaccine1",
                                              pValueCutoff = 0.05,
                                              showPlots = FALSE, 
                                              returnSeuratObjectOrPlots = "Plots", 
                                              includeDepletions = FALSE))
  testthat::expect_true(length(plots) == 2)
  testthat::expect_true(typeof(plots$model_coefficients) == "list")
  
  #test quasipoisson-based cluster enrichment
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  seuratObj$cDNA_ID <- rep(1:4, ncol(seuratObj)) 
  seuratObj$Vaccine <- rep(c("Vaccine1", "Vaccine2"), each = ncol(seuratObj)/2) 
  seuratObj$SubjectId <- rep(1:2, each = ncol(seuratObj)/2) 
  
  seuratObj$Vaccine <- ifelse(seuratObj$ClusterNames_0.2 == 1, 
                              sample(c("Vaccine1", "Vaccine2"), 
                                     size = sum(seuratObj$ClusterNames_0.2 == 1),
                                     replace = TRUE, 
                                     prob = c(0.9, 0.1)), 
                              seuratObj$Vaccine)
  
  
  testthat::expect_no_error(seuratObj <- CalculateClusterEnrichmentPairwise(seuratObj,
                                              subjectField = 'SubjectId',
                                              clusterField = 'ClusterNames_0.2',
                                              biologicalReplicateGroupingVariables = c("cDNA_ID"),
                                              treatmentField = "Vaccine",
                                              referenceValue = "Vaccine1",
                                              pValueCutoff = 0.05,
                                              showPlots = FALSE, 
                                              returnSeuratObjectOrPlots = "SeuratObject", 
                                              includeDepletions = TRUE))
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
