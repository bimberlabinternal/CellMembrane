context("scRNAseq")

test_that("CSCore works as expected", {
    seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
    csResults <- RunCsCore(seuratObj, genes_selected = Seurat::VariableFeatures(seuratObj))

    expect_equal(length(csResults), 51)
    expect_equal(paste0(sort(csResults[[40]]), collapse = ','), "ANKRD27,C10orf90,CDK5RAP2,CENPT,COA1,CRTC1,ELMSAN1,ENSMMUG00000005997,LACTB,MAP3K14,NTN4,PRICKLE3,PRPF40B,RNFT1,S100PBP,SLC25A42,STK38L,SYNE2,TOP3B,USP9Y,WDR11")
})
