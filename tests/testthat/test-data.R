library(ggplot2)

context("scRNAseq")

test_that("Cell cycle gene data OK", {
	expect_equal(length(cc.genes), 97)

	expect_equal(length(g2m.genes.orig), 200)
})

test_that("CD genes data OK", {
	expect_equal(nrow(cdGenes), 394)
})

test_that("CD genes translation works", {
	expect_equal(RenameUsingCD(c('ITGB1', 'EGLN3', 'Other', 'CD52')), c('ITGB1 (CD29)', 'EGLN3', 'Other', 'CD52 (CDW52)'))
})

# This primarily serves to demonstrate the code runs without overt errors
test_that("DotPlot works", {
	set.seed(CellMembrane::GetSeed())
	seuratObj <- readRDS('../testdata/seuratOutput.rds')

	seuratObj$XField <- sample(LETTERS[1:3], size = nrow(seuratObj@meta.data), replace = TRUE)
	seuratObj$YField <- sample(LETTERS[4:8], size = nrow(seuratObj@meta.data), replace = TRUE)
	seuratObj$Color <- as.factor(sample(LETTERS[9:10], size = nrow(seuratObj@meta.data), replace = TRUE))
	seuratObj$Cluster <- sample(LETTERS[16:21], size = nrow(seuratObj@meta.data), replace = TRUE)
	seuratObj$Group1 <- sample(LETTERS[25:26], size = nrow(seuratObj@meta.data), replace = TRUE)
	seuratObj$NormField <- sample(LETTERS[16:21], size = nrow(seuratObj@meta.data), replace = TRUE)

	metacounts <- ConstructEnrichmentDataFrameAndDoStatistics(seuratObj = seuratObj,
															  xField = 'XField',
															  yField = 'YField',
															  colorField = 'Color',
															  extraGroupingFields = c('Group1'),
															  normalizationField = 'NormField')

	expect_equal(nrow(metacounts), 30)
	
	P1 <- MakeEnrichmentDotPlot(seuratObj = seuratObj,
															xField = 'XField',
															yField = 'YField',
															colorField = 'Color',
															extraGroupingFields = c('Group1'),
															normalizationField = 'NormField')

	P1 + facet_grid(. ~ Group1)

})