test_that("ADT triage works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  seuratObj <- CellMembrane::AppendCiteSeq(seuratObj, unfilteredMatrixDir = '../testdata/raw_feature_bc_matrix')
  ad <- seuratObj@assays[["ADT"]]
  ad <- Seurat::NormalizeData(ad, normalization.method = 'CLR', margin = 1, verbose = FALSE)
  seuratObj[["ADT"]] <- ad
  seuratObj <- triageADTsAndClassifyCells(seuratObj, adtwhitelist = c("TotalSeq-C-159", "TotalSeq-C-146", "TotalSeq-C-033"),
                                          libraryIdentificationColumn = "orig.ident", max_frac_under_threshold = 0.999)
  expect_equal(seuratObj$`reason_TotalSeq-C-033` |> unique(), "Antimode too high")
  expect_equal(seuratObj@meta.data[seuratObj$`TotalSeq-C-159_cellcall`=="Positive",] |> rownames(),
               c("AACTCAGTCGGCATCG", "AAGGTTCTCTTGACGA", "ACAGCTATCTCCCTGA","CCTACCACAATCCAAC", "TAGACCACAAGCGAGT", "TCGCGTTGTGAGCGAT"))
  expect_equal(seuratObj@meta.data[seuratObj$`TotalSeq-C-146_cellcall`=="Positive",] |> rownames(),
               c("AAGGTTCTCTTGACGA", "CCTACCACAATCCAAC", "GGCGACTAGGATTCGG", "TCGCGTTGTGAGCGAT"))
  expect_equal(seuratObj$`antimode_location_TotalSeq-C-033` |> unique(), 4.232, tolerance = 0.001)
  expect_equal(seuratObj$`antimode_location_TotalSeq-C-159` |> unique(), 4.969, tolerance = 0.001)
  expect_equal(seuratObj$`antimode_location_TotalSeq-C-146` |> unique(), 4.190, tolerance = 0.001)
})
