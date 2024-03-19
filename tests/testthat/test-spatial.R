context("spatialRNAseq")

test_that("Q3 normalization works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  
  testthat::expect_no_error(seuratObj[["RNA"]])
  seuratObj<- RUVg_Housekeeping_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["RUVg"]])
  expect_equal(max(Seurat::GetAssayData(seuratObj, assay = 'RUVg', slot = 'counts')), 2335)
  testthat::expect_error(object = NanoString_Housekeeping_Normalization(seuratObj),"Geometric mean of housekeeping counts for one or more samples is zero")
  testthat::expect_error(object = Q3_Normalization(seuratObj),"Error: 75th percentile is zero for some samples.")

  seuratObj <- Seurat::SetAssayData(seuratObj, assay = 'RNA', slot = 'counts', new.data = SeuratObject::as.sparse(as.matrix(Seurat::GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')) + 1))
  
  seuratObj<- Q3_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["Q3"]])
  testthat::expect_equal(max(Seurat::GetAssayData(seuratObj, assay = 'Q3', slot = 'counts')), 1663)

  seuratObj<- NanoString_Housekeeping_Normalization(seuratObj)
  testthat::expect_no_error(seuratObj[["Housekeeping"]])
  testthat::expect_equal(round(max(Seurat::GetAssayData(seuratObj, assay = 'Housekeeping', slot = 'counts'))), 76)

})

test_that("Spatial cell substructure detection works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  set.seed(GetSeed())
  seuratObj@meta.data$x_FOV_px <- c(rnorm(n = 1000, mean = 1000, sd = 100), 
                                    rnorm(n = 500, mean = 2000, sd = 100),
                                    rnorm(n = 57, mean = 1500, sd = 2000))
  seuratObj@meta.data$y_FOV_px <- c(rnorm(n = 1000, mean = 1000, sd = 100), 
                                    rnorm(n = 500, mean = 2000, sd = 100),
                                    rnorm(n = 57, mean = 1500, sd = 2000))
  seuratObj@meta.data$cell_type <- c(rep("Bcell", 1000), 
                                     rep("Bcell", 500), 
                                     rep("Myeloid", 57))
  seuratObj@meta.data$fov <- 1
  metadata <- DetectCellStructuresBasedOnCellType(seuratObjectMetadata = seuratObj@meta.data, 
                                      cellTypeField = "cell_type",
                                      fovField =  "fov", 
                                      cellTypeWhiteList = "Bcell", 
                                      substructureMetaDataFieldName = "BCF", 
                                      summarizeLocalResults = TRUE
                                      )
  #ensure metadata field detection errors
  testthat::expect_error(DetectCellStructuresBasedOnCellType(seuratObjectMetadata = seuratObj@meta.data, 
                                                             cellTypeField = "cell_type",
                                                             fovField =  "aFieldThatDoesNotExist", 
                                                             cellTypeWhiteList = "Bcell", 
                                                             substructureMetaDataFieldName = "BCF", 
                                                             summarizeLocalResults = TRUE)
                         )
  #test that cell types exist
  testthat::expect_error(DetectCellStructuresBasedOnCellType(seuratObjectMetadata = seuratObj@meta.data, 
                                                               cellTypeField = "cell_type",
                                                               fovField =  "fov", 
                                                               cellTypeWhiteList = "aCellTypeThatDoesn'tExist", 
                                                               substructureMetaDataFieldName = "BCF", 
                                                               summarizeLocalResults = TRUE)
                           )
  
  #lax test to make sure the code ran at all
  testthat::expect_true(c("Within_Local_BCF" %in% colnames(metadata)))
  #more stringent test to make sure the code yields two identified clusters of cells
  #This should work 100% of the time, but I can't ensure the seed passes to (h)dbscan correctly, so we can ditch this if it causes problems.
  testthat::expect_true(all(c(0,1,2) %in% unique(metadata[,"Local_BCF"])))
})