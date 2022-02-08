# Disabled since this is very slow w/o GPUs
# test_that("cellbender works as expected", {
#   dataDir <- '../testdata/CellRanger3/raw_feature_bc_matrix'
#   rawFeatureMatrix <- Seurat::Read10X(data.dir = dataDir, strip.suffix = T)
#   seuratObj <- RunCellBender(rawFeatureMatrix, epochs = 1)
#   print(seuratObj)
# })