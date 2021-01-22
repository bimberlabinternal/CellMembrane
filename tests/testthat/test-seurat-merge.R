context("scRNAseq")

test_that("Seurat-merge works as expected", {
  data <- list(
  'Set1' = '../testdata/CellRanger2/raw_gene_bc_matrices/cellRanger-3204293',
  'Set2' = '../testdata/CellRanger3/raw_feature_bc_matrix'
  )

  seuratObjs <- list()
  for (datasetName in names(data)) {
    seuratObjs[[datasetName]] <- ReadAndFilter10xData(testthat::test_path(data[[datasetName]]), datasetName, emptyDropNIters=1000)
  }

  expect_equal(length(seuratObjs), 2)

  #Simple method
  seuratObj <- MergeSeuratObjs(seuratObjs, data, method = NULL)
  expect_equal("RNA", Seurat::DefaultAssay(seuratObj))
  expect_equal("simple", seuratObj@misc$MergeMethod)

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
  expect_error(MergeSeuratObjs(seuratObjs, data, method = 'bad'))

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)

  #Not all genes will be present:
  expect_equal(length(intersect(rownames(seuratObj), expectedSpike)), 46)

  #barcodes should have prefix:
  expect_equal(sum(!grepl(colnames(seuratObj), pattern = '^Set')), 0)
  
  #All spike genes should be present
  expect_equal(sort(intersect(rownames(seuratObj), expectedSpike)), sort(expectedSpike))
  
})

