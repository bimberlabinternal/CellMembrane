context("scRNAseq")

test_that("RunDecoupleR works as expected", {
  # This failed on github runners with a non-informative 'job terminated' error.
  if (version$minor > 4.3) {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    seuratObj <- RunDecoupleR(seuratObj)

    expect_true('tfsulm' %in% names(seuratObj@assays))
    expect_equal(nrow(seuratObj@assays$tfsulm), 726)
    expect_equal(ncol(seuratObj@assays$tfsulm), ncol(seuratObj))

    CellMembrane::PlotTfData(seuratObj, groupField = 'Phase')
  } else {
    print('Skipping decoupleR tests')
  }
})


