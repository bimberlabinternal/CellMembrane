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