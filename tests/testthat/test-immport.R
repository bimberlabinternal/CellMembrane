context("ImmPort")

test_that("GeneToImmPortPathway rejects invalid input", {
  expect_error(GeneToImmPortPathway(123))
  expect_error(GeneToImmPortPathway(c("AHNAK", "CD3E")))
})

test_that("GeneToImmPortPathway returns NA for unknown genes", {
  result <- GeneToImmPortPathway("THISISNOTAREALGENE_XYZ")
  expect_true(is.na(result))
})

test_that("GeneToImmPortPathway is case-insensitive and returns pathways", {
  result_upper <- GeneToImmPortPathway("CD3E")
  result_lower <- GeneToImmPortPathway("cd3e")
  expect_identical(result_upper, result_lower)
  expect_true(length(result_upper) >= 1)
  expect_true(is.character(result_upper))
})

test_that("GetImmPortGenes returns a non-empty character vector of uppercase symbols", {
  immport_genes <- GetImmPortGenes()
  expect_true(is.character(immport_genes))
  expect_true(length(immport_genes) > 0)
  expect_true(all(immport_genes == toupper(immport_genes)))
})

test_that("SubsetSeuratToImmPortGenes reduces feature space", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  original_n_features <- nrow(seuratObj)

  subset_obj <- SubsetSeuratToImmPortGenes(seuratObj)
  retained_features <- rownames(subset_obj)

  expect_true(length(retained_features) > 0)
  expect_true(length(retained_features) <= original_n_features)
  expect_true(all(toupper(retained_features) %in% GetImmPortGenes()))
})

test_that("SubsetSeuratToImmPortGenes errors on missing assay", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
  expect_error(SubsetSeuratToImmPortGenes(seuratObj, assayName = "NONEXISTENT_ASSAY"))
})

test_that("PermutationTestForImmPortPathwayOccurrences rejects invalid input", {
  expect_error(PermutationTestForImmPortPathwayOccurrences(123))
  expect_error(PermutationTestForImmPortPathwayOccurrences(character(0)))
})

test_that("PermutationTestForImmPortPathwayOccurrences errors when no genes overlap ImmPort universe", {
  expect_error(
    PermutationTestForImmPortPathwayOccurrences(c("THISISNOTAREALGENE_XYZ", "ANOTHERFAKEGENE_ABC")),
    regexp = "None of the genes"
  )
})

test_that("PermutationTestForImmPortPathwayOccurrences returns expected data.frame structure", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  # use head() for a deterministic gene set that doesn't depend on an outer seed
  gene_sample <- head(GetImmPortGenes(), 50)
  result <- PermutationTestForImmPortPathwayOccurrences(gene_sample, nPermutations = 100, seed = CellMembrane::GetSeed())

  expect_s3_class(result, "data.frame")
  expect_true(all(c("Pathway", "ObservedCount", "MeanNullCount", "PValue", "PValueAdj") %in% colnames(result)))
  expect_true(nrow(result) > 0)
  expect_true(all(result$PValue    >= 0 & result$PValue    <= 1))
  expect_true(all(result$PValueAdj >= 0 & result$PValueAdj <= 1))
  expect_true(all(result$ObservedCount > 0))
  # result must be sorted by ascending PValue
  expect_true(all(diff(result$PValue) >= 0))
})

test_that("PermutationTestForImmPortPathwayOccurrences is reproducible with a fixed seed", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  gene_sample <- head(GetImmPortGenes(), 50)
  result1 <- PermutationTestForImmPortPathwayOccurrences(gene_sample, nPermutations = 100, seed = CellMembrane::GetSeed())
  result2 <- PermutationTestForImmPortPathwayOccurrences(gene_sample, nPermutations = 100, seed = CellMembrane::GetSeed())
  expect_identical(result1, result2)
})

test_that("PermutationTestForImmPortPathwayOccurrences input is case-insensitive", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  gene_sample <- head(GetImmPortGenes(), 30)
  result_upper <- PermutationTestForImmPortPathwayOccurrences(toupper(gene_sample), nPermutations = 100, seed = CellMembrane::GetSeed())
  result_lower <- PermutationTestForImmPortPathwayOccurrences(tolower(gene_sample), nPermutations = 100, seed = CellMembrane::GetSeed())
  expect_identical(result_upper, result_lower)
})

test_that("PermutationTestForImmPortPathwayOccurrences detects enrichment signal", {
  CellMembrane::SetSeed(CellMembrane::GetSeed())
  # all genes from one pathway should produce a low (significant) p-value for that pathway
  gene_to_pathway <- GetImmPortGeneMapping()
  # find a pathway with at least 10 genes
  pathway_sizes   <- vapply(names(gene_to_pathway), function(g) length(gene_to_pathway[[g]]), integer(1))
  large_pathways  <- names(which(table(unlist(gene_to_pathway)) >= 10))
  target_pathway  <- large_pathways[[1]] #grab "Activation of innate immune response"
  focused_genes   <- names(which(vapply(gene_to_pathway, function(pws) target_pathway %in% pws, logical(1))))
  result <- PermutationTestForImmPortPathwayOccurrences(focused_genes, nPermutations = 500, seed = CellMembrane::GetSeed())
  target_row <- result[result$Pathway == target_pathway, ]
  expect_true(nrow(target_row) == 1)
  expect_true(target_row$PValue < 0.05)
})
