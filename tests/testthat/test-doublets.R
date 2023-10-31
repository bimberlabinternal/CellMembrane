context("scRNAseq")

test_that("Doublet detection works as expected", {
  set.seed(CellMembrane::GetSeed())

  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

  fn <- 'doublets.txt'
  seuratObj <- FindDoublets(seuratObj, rawResultFile = fn)
  
  expect_equal(17, sum(seuratObj$scDblFinder.class == 'doublet'), tolerance = 3)
  expect_equal(1540, sum(seuratObj$scDblFinder.class == 'singlet'), tolerance = 3)
  
  df <- read.table(fn, header = T, sep = '\t')
  expect_equal(1557, nrow(df))
  unlink(fn)

  seuratObj <- FindDoublets(seuratObj, rawResultFile = fn, dropDoublets = TRUE)
  expect_equal(1540, ncol(seuratObj), tolerance = 3)
})

