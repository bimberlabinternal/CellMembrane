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
  assayObj <- GetAssay(seuratObj, assay = 'RNA')
  expect_equal(nrow(seuratObj), length(slot(assayObj, GetAssayMetadataSlotName(assayObj))$GeneId))
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


test_that("seurat barcode duplicate code works as expected", {
	seuratObj1 <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj1 <- Seurat::DietSeurat(seuratObj1)
	seuratObj1$DatasetId <- 12345
	seuratObj1$BarcodePrefix <- 12345
	seuratObj1 <- Seurat::RenameCells(object = seuratObj1, add.cell.id = 12345)
	
	seuratObj2 <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj2 <- Seurat::DietSeurat(seuratObj2)
	seuratObj2$DatasetId <- 23456
	seuratObj2$BarcodePrefix <- 23456
	seuratObj2 <- Seurat::RenameCells(object = seuratObj2, add.cell.id = 23456)

	nn <- colnames(seuratObj2)
	nn[1:10] <- colnames(seuratObj1)[1:10]
	seuratObj2 <- Seurat::RenameCells(object = seuratObj2, new.names = nn)
	
	toMerge <- list('12345' = seuratObj1, '23456' = seuratObj2)
	sm <- MergeSeuratObjs(toMerge, projectName = 'TestMerge')

	expect_equal(ncol(sm), ncol(seuratObj1) +  ncol(seuratObj2) - 10)
	
	# Repeat where barcode prefix added by merge:
	seuratObj1 <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj1 <- Seurat::DietSeurat(seuratObj1)

	seuratObj2 <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj2 <- Seurat::DietSeurat(seuratObj2)

	toMerge <- list('12345' = seuratObj1, '23456' = seuratObj2)
	sm <- MergeSeuratObjs(toMerge, projectName = 'TestMerge')
	expect_equal(ncol(sm), ncol(seuratObj1) +  ncol(seuratObj2))	
})

test_that("Assumptions about Seurat layers are true", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))

  assayA <- SeuratObject::CreateAssayObject(Seurat::GetAssayData(seuratObj[['RNA']], layer = 'counts'))
  colnames(assayA) <- paste0('A_', colnames(assayA))
  expect_true(inherits(assayA, 'Assay'))
  expect_false(inherits(assayA, 'Assay5'))

  assayB <- SeuratObject::CreateAssayObject(Seurat::GetAssayData(seuratObj[['RNA']], layer = 'counts'))
  colnames(assayB) <- paste0('B_', colnames(assayB))
  
  assayC <- SeuratObject::CreateAssayObject(Seurat::GetAssayData(seuratObj[['RNA']], layer = 'counts'))
  colnames(assayC) <- paste0('C_', colnames(assayC))

  assay5A <- SeuratObject::CreateAssay5Object(Seurat::GetAssayData(assayA, layer = 'counts'))
  assay5B <- SeuratObject::CreateAssay5Object(Seurat::GetAssayData(assayB, layer = 'counts'))
  assay5C <- SeuratObject::CreateAssay5Object(Seurat::GetAssayData(assayC, layer = 'counts'))
  expect_true(inherits(assay5A, 'Assay5'))
  expect_false(inherits(assay5A, 'Assay'))
  
  merge1 <- merge(assayA, assayB)
  expect_equal(c('counts', 'data'), SeuratObject::Layers(merge1))
  expect_true(inherits(merge1, 'Assay'))
  expect_false(inherits(merge1, 'Assay5'))
  
  merge5 <- merge(assay5A, assay5B, collapse = FALSE)
  expect_equal(c('counts.1', 'counts.2'), SeuratObject::Layers(merge5))
  merge5 <- SeuratObject::JoinLayers(merge5)
  expect_equal(c('counts'), SeuratObject::Layers(merge5))
  expect_true(inherits(merge5, 'Assay5'))
  expect_false(inherits(merge5, 'Assay'))

  seurat5a <- SeuratObject::CreateSeuratObject(project = 'Proj1', counts = SeuratObject::SetAssayData(assay5A, layer = 'layer', new.data = Seurat::GetAssayData(assay5A, layer = 'counts')))
  seurat5b <- SeuratObject::CreateSeuratObject(project = 'Proj2', counts = SeuratObject::SetAssayData(assay5B, layer = 'layer', new.data = Seurat::GetAssayData(assay5B, layer = 'counts')))
  seurat5c <- SeuratObject::CreateSeuratObject(project = 'Proj3', counts = SeuratObject::SetAssayData(assay5C, layer = 'layer', new.data = Seurat::GetAssayData(assay5C, layer = 'counts')))
  
  merge5v2 <- merge(seurat5a, list(Proj2 = seurat5b, Proj3 = seurat5c), collapse = FALSE)
   
  expect_equal(c('counts.Proj1', 'counts.Proj2', 'counts.Proj3', 'layer.Proj1', 'layer.Proj2', 'layer.Proj3'), SeuratObject::Layers(merge5v2))
  
  merge5v2Joined <- SeuratObject::JoinLayers(merge5v2)
  expect_equal(c('counts', 'layer.Proj1', 'layer.Proj2', 'layer.Proj3'), SeuratObject::Layers(merge5v2Joined))
  
  merge5v2Joined <- SeuratObject::JoinLayers(merge5v2, layers = .FindLayersToJoin(merge5v2, 'RNA'))
  expect_equal(c('counts', 'layer'), sort(SeuratObject::Layers(merge5v2Joined)))
})


