context("scRNAseq")

test_that("Pathway scoring works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  #pull Hallmark genesets from Msigdb, but only keep one for scoring. 
  geneSets <- escape::getGeneSets(library = "H")
  #this should fail if the geneSets provided are not named lists of vectors
  expect_error(AlternativeSsgseaSeurat(seuratObj, geneSets[[1]]))
  
  geneSet <- geneSets[1]
  seuratObj <- AlternativeSsgseaSeurat(seuratObj, geneSet)
  
  #check that the assay was added to the seurat object
  expect_true("ssGSEA.alternative" %in% Seurat::Assays(seuratObj))
  results <- GetAssayData(seuratObj, assay = "ssGSEA.alternative")
  #check that the results are approximately correct. This calculation should be deterministic. 
  expect_equal(colMeans(matrix(results)), 616297.8, tolerance = 1)
})