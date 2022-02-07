test_that("Seurat-merge works as expected", {
  dataDir <- '../testdata/CellRanger3/raw_feature_bc_matrix'
  rawFeatureMatrix <- Seurat::Read10X(data.dir = dataDir, strip.suffix = T)
  seuratObj <- RunCellBender(rawFeatureMatrix, epochs = 2)
  print(seuratObj)
})