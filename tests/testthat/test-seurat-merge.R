context("scRNAseq")

test_that("Seurat-merge works as expected", {
  set.seed(CellMembrane::GetSeed())

  data <- list(
    'Set1' = '../testdata/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
    'Set2' = '../testdata/CellRanger3/raw_feature_bc_matrix'
  )

  seuratObjs <- list()
  for (datasetId in names(data)) {
    seuratObj <- ReadAndFilter10xData(testthat::test_path(data[[datasetId]]), datasetId = datasetId, emptyDropNIters=1000)
    expect_true('BarcodePrefix' %in% colnames(seuratObj@meta.data))
    expect_true('DatasetId' %in% colnames(seuratObj@meta.data))
    expect_false('DatasetName' %in% colnames(seuratObj@meta.data))
    seuratObjs[[datasetId]] <- seuratObj
  }

  expect_equal(length(seuratObjs), 2)

  seuratObj <- MergeSeuratObjs(seuratObjs, projectName = 'project')
  expect_true('BarcodePrefix' %in% colnames(seuratObj@meta.data))
  expect_true('DatasetId' %in% colnames(seuratObj@meta.data))

  expect_equal("RNA", Seurat::DefaultAssay(seuratObj))

  #Gene IDs preserved:
  expect_equal(nrow(seuratObj), length(seuratObj@assays$RNA@meta.features$GeneId))
  geneIds <- GetGeneIds(seuratObj, c('HES4', 'CALML6'))
  names(geneIds) <- NULL
  expect_equal(geneIds, c('ENSMMUG00000001817', 'ENSMMUG00000012392'))

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)

  # NOTE: this might not be deterministic.  expect about 7987
  print(paste0('cells: ', ncol(seuratObj)))
  expect_equal(7987, ncol(seuratObj), tolerance = 10)
  
  #Invalid method
  expect_error(MergeSeuratObjs(seuratObjs))

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)
})

