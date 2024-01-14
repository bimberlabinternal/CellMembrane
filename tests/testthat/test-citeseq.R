context("scRNAseq")

test_that("Cite-Seq Normalization Works", {
	adts <- read.table('../testdata/raw_feature_bc_matrix/validADTS.32851.metadata.txt', sep = '\t', header = T, fill = TRUE)
	adts$rowname <- adts$tagname
	mat <- .LoadCiteSeqData('../testdata/raw_feature_bc_matrix', adtWhitelist = adts$tagname, featureMetadata = adts)

	df <- data.frame(cellbarcode = colnames(mat), count = colSums(as.matrix(Seurat::GetAssayData(mat, slot = 'counts'))))
	df <- dplyr::arrange(df, desc(count))
	
	seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	sc <- Seurat::GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')[,1:1000]
	colnames(sc) <- df$cellbarcode[1:1000]
	seuratObj <- Seurat::CreateSeuratObject(counts = sc)
	
	ret <- .NormalizeDsbWithEmptyDrops(seuratObj, unfilteredAdtAssay = mat, emptyDropNIters = 1000)
	expect_equal(nrow(Seurat::GetAssayData(ret, slot = 'counts')), nrow(adts))
})

test_that("Cite-Seq Append Works", {
	#Reduce size, touch up data:
	seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj <- seuratObj[1:1000,1:100]
	
	set.seed(1234)
	dummyCounts <- as.matrix(DropletUtils:::simCounts(ngenes = 8, nempty = 20))
	
	citeseqData1 <- dummyCounts[1:4,1:50]
	colnames(citeseqData1) <- colnames(seuratObj)[1:50]
	colnames(citeseqData1) <- sample(colnames(citeseqData1), length(colnames(citeseqData1)), replace = FALSE)
	rownames(citeseqData1) <- paste0('ADT', 1:4)
	inputPath1 <- './input1'
	if (dir.exists(inputPath1)) {
	  unlink(inputPath1, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath1, Seurat::as.sparse(citeseqData1))
	
	citeseqData2 <- dummyCounts[3:8,41:100]
	colnames(citeseqData2) <- colnames(seuratObj)[41:100]
	colnames(citeseqData2) <- sample(colnames(citeseqData2), length(colnames(citeseqData2)), replace = FALSE)
	rownames(citeseqData2) <- paste0('ADT', 3:8)
	inputPath2 <- './input2'
	if (dir.exists(inputPath2)) {
	  unlink(inputPath2, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath2, Seurat::as.sparse(citeseqData2))

	citeseqData3 <- dummyCounts[1:2,61:100]
	colnames(citeseqData3) <- colnames(seuratObj)[61:100]
	colnames(citeseqData3) <- sample(colnames(citeseqData3), length(colnames(citeseqData3)), replace = FALSE)
	rownames(citeseqData3) <- c('NewMarkerName', 'unmapped')
	inputPath3 <- './input3'
	if (dir.exists(inputPath3)) {
	  unlink(inputPath3, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath3, Seurat::as.sparse(citeseqData3))

	# #First w/o prefix:
	seuratObj$DatasetId <- 'Set1'
	seuratObjCite <- AppendCiteSeq(seuratObj = seuratObj, unfilteredMatrixDir = inputPath1, datasetId = NULL, normalizeMethod = NULL)
	expect_equal(colnames(seuratObjCite@assays$ADT), colnames(seuratObj@assays$RNA))
	expect_equal(rownames(seuratObjCite@assays$ADT), rownames(citeseqData1)[rownames(citeseqData1) != 'unmapped'])
	data <- Seurat::GetAssayData(seuratObjCite, assay = 'ADT', slot = 'counts')
	expect_equal(max(data[,51:100]), 0) #These have no data
	expect_equal(sum(data[,1:50] == 0), 32)
	
	#Rename cells with a barcodeprefix:
	seuratObj@meta.data$DatasetId <- c(rep('12345', 50), rep('67890', 50))
	seuratObj <- Seurat::RenameCells(object = seuratObj, new.names = paste0(seuratObj@meta.data$DatasetId, '_', colnames(seuratObj)))
	
	#With prefix, no assay present
	seuratObjCite <- AppendCiteSeq(seuratObj = seuratObj, unfilteredMatrixDir = inputPath1, datasetId = '12345', normalizeMethod = NULL)
	expect_equal(colnames(seuratObjCite@assays$ADT), colnames(seuratObj@assays$RNA))
	expect_equal(rownames(seuratObjCite@assays$ADT), rownames(citeseqData1)[rownames(citeseqData1) != 'unmapped'])
	data <- Seurat::GetAssayData(seuratObjCite, assay = 'ADT', slot = 'counts')
	expect_equal(max(data[,51:100]), 0) #These have no data
	expect_equal(sum(data[,1:50] == 0), 32)
	
	# # #Now add again, existing assay present:
	seuratObjCite2 <- AppendCiteSeq(seuratObj = seuratObjCite, unfilteredMatrixDir = inputPath2, datasetId  = '67890', minRowSum = 0, normalizeMethod = NULL)
	d1 <- Seurat::GetAssayData(seuratObjCite, 'ADT', slot = 'counts')
	d2 <- Seurat::GetAssayData(seuratObjCite2, 'ADT', slot = 'counts')
	expect_equal(d2[1:4,1:50], d1[1:4,1:50]) #The original data
	expect_equal(max(d2[5:8,1:50]), 0)
	expect_equal(nrow(d2), 8)
	expect_equal(max(d2[1:2,51:100]), 0) #The original data
	
	# #This time with: differing sets of features:
	seuratObjCite3 <- AppendCiteSeq(seuratObj = seuratObjCite, unfilteredMatrixDir = inputPath3, datasetId = '67890', minRowSum = 0, normalizeMethod = NULL)
	d1 <- Seurat::GetAssayData(seuratObjCite, 'ADT', slot = 'counts')
	d2 <- Seurat::GetAssayData(seuratObjCite3, 'ADT', slot = 'counts')
	expect_equal(nrow(d2), 5) #our new marker is added
	expect_equal(d2[1:4,1:50], d1[1:4,1:50]) #The original data
	expect_equal(max(d2[1:4,51:100]), 0) #not present in original
	expect_equal(max(d2[5,51:100]), 59)
	
	unlink(inputPath1, recursive = T)
	unlink(inputPath2, recursive = T)
	unlink(inputPath3, recursive = T)
})

test_that("ADT Rename Works", {
	#Reduce size, touch up data:
	seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS('../testdata/seuratOutput.rds')))
	seuratObj <- seuratObj[1:1000,1:100]
	
	set.seed(1234)
	dummyCounts <- as.matrix(DropletUtils:::simCounts(ngenes = 8, nempty = 20))
	
	citeseqData1 <- dummyCounts[1:4,1:50]
	colnames(citeseqData1) <- colnames(seuratObj)[1:50]
	rownames(citeseqData1) <- paste0('ADT', 1:4)
	inputPath1 <- './input1'
	if (dir.exists(inputPath1)) {
		unlink(inputPath1, recursive = T)
	}
	DropletUtils::write10xCounts(inputPath1, Seurat::as.sparse(citeseqData1))
	
	#Test feature metadata add:
	meta <- data.frame(rowname = c('ADT1'), sequence = c('CTCCTCTGCAATTAC'), markername = c('RenamedFeature'), otherfield = c('Value1'))

	seuratObj$DatasetId <- 'Set1'
	seuratObjCiteMeta <- AppendCiteSeq(seuratObj = seuratObj, unfilteredMatrixDir = inputPath1, datasetId = NULL, featureMetadata = meta, normalizeMethod = NULL)
	
	expect_equal(colnames(seuratObjCiteMeta@assays$ADT), colnames(seuratObj@assays$RNA))
	expect_true('RenamedFeature' %in% rownames(seuratObjCiteMeta@assays$ADT))
	expect_equal(Seurat::GetAssay(seuratObjCiteMeta, assay = 'ADT')@meta.features$otherfield, c('Value1', NA, NA, NA))

	unlink(inputPath1, recursive = T)
})
	