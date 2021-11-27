#' @title PlotImmuneMarkers
#'
#' @description Generate a set of Seurat FeaturePlots for common immune cell markers
#' @param seuratObj A seurat object
#' @param reductions Vector of reduction(s) to use
#' @export
#' @import Seurat
PlotImmuneMarkers <- function(seuratObj, reductions = c('tsne', 'umap')) {
	reductions <- intersect(reductions, names(seuratObj@reductions))
	if (length(reductions) == 0) {
		stop('None of the requested reductions are present!')
	}

	#ENSMMUG00000003532=CD8b
	PlotMarkerSeries(seuratObj, reductions, c('CD8A', 'CD8B', 'ENSMMUG00000003532', 'CD4', 'IL7R', 'CD3D', 'CD3E','CD3G'), 'CD8/CD4 Markers')

	#Eff v. Mem:
	#IL7R = CD127
	#IL2RA = CD25
	#PTPRC = CD45
	#SELL = CD62-L / CD-197
	PlotMarkerSeries(seuratObj, reductions, c('CCR7', 'SELL', 'GZMB', 'CCR5', 'IL2RA', 'PTPRC', 'IL7R', 'CTLA4', 'FAS', 'CD28', 'CD27'), 'Effector vs. Memory')

	#CD8 Activation
	# XCL1 = ENSMMUG00000060218
	# CCL4 = ENSMMUG00000008111
	# LOC100423131 = XCL1, ENSMMUG00000013779, Lymphotactin
	# LOC100430627 = CCL4, ENSMMUG00000008111
	# LOC100426632 = C-C motif chemokine 4
	# LOC100426537 = C-C motif chemokine 3-like
	# CD154 = CD40L
	# LAMP1 = CD107a, cytotoxic capacity
	PlotMarkerSeries(seuratObj, reductions, c('IFNG', 'CD69', 'TNF', 'NFKBID', 'LTB', 'TNFRSF9', 'CCL4L1', 'NR4A3', 'TNFSF14', 'CD82', 'PIGT', 'IRF8', 'IRF4', 'RGCC', 'PD1', 'PDCD1', 'TNFAIP3', 'LOC100423131', 'LOC100430627', 'ENSMMUG00000013779', 'XCL1', 'ENSMMUG00000060218', 'CCL4', 'ENSMMUG00000008111', 'CCL3', 'PLEK', 'NR4A2', 'LOC100426537', 'LOC114673087', 'KLF10', 'GADD45B', 'CD154', 'LAMP1'), 'CD8 Activation Markers')

	PlotMarkerSeries(seuratObj, reductions, c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'), 'Cytotoxicity')

	PlotMarkerSet(seuratObj, reductions, 'B-cell Markers', c('MS4A1', 'CD79A', 'CD74', 'DRA'))

	PlotMarkerSet(seuratObj, reductions, 'Monocyte', c('LYZ', 'CST3', 'S100A6', 'VIM'))

	# ITGB3 = CD61
	PlotMarkerSet(seuratObj, reductions, 'Platelet', c('ITGB3'))

	# MRC1 = CD206
	# ITGAM = CD11b, AM=CD11b-, Non-AM=CD11b+
	# Non-AMs: CD16+/CD206-/HLA-DR+/CD11b+
	# M1/M2: CCR7 (M1), CD163 (M2), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7050058/
	# CD163, CD68 = macrophage markers
	PlotMarkerSeries(seuratObj, reductions, c('CD14', 'FCGR3A', 'S100A8', 'S100A6', 'MARCO', 'MRC1', 'CD163', 'CHIT1', 'APOBEC3A', 'ITGAM', 'HLA-DRB1', 'CCR7', 'CD68'), 'Myeloid')

	# IL3RA = CD123 / pDC
	# CLEC4C = CD303 / pDC
	# CD1c = mDC
	# THBD = CD141 / mDC
	# CD80, CD86 = co-stimulatory. low expression = tolerogenic
	PlotMarkerSeries(seuratObj, reductions, c('CD14', 'FCGR3A', 'IL3RA', 'CLEC4C', 'CD1c', 'THBD', 'CD80', 'CD86', 'TGFB1'), 'DCs')

	PlotMarkerSet(seuratObj, reductions, 'Regulatory Markers', c('CD4', 'FOXP3', 'IL2RA'))

	PlotMarkerSet(seuratObj, reductions, 'Treg17', c('RORA', 'RORB', 'RORC', 'IL4', 'STAT3'))

	PlotMarkerSeries(seuratObj, reductions, c('ZBTB16', 'DPP4'), 'MAIT')

	PlotMarkerSeries(seuratObj, reductions, c('TBX21', 'GATA3', 'RORC', 'FOXP3', 'BCL6', 'EOMES', 'TOX'), 'Transcription Factors')

	PlotMarkerSeries(seuratObj, reductions, c('TIGIT', 'CTLA4', 'BTLA', 'PDCD1', 'CD274'), 'Inhibitory Markers')

	#DAP10/12
	PlotMarkerSet(seuratObj, reductions, 'Signaling', c('HCST', 'TYROBP', 'SYK', 'ZAP70'))

	#LILR/KIR:
	PlotMarkerSeries(seuratObj, reductions, c('LILRA5','LILRA6','LILRB4','LILRB5','KIR2DL4','KIR3DX1', 'MAMU-KIR', 'KIR2DL4', 'KIR3DL2'), 'LILR/KIR')

	PlotMarkerSeries(seuratObj, reductions, c('FCGR1A','FCGR2A','FCGR2B','FCGR3'), 'FCGR')

	#Cytokines
	cytokines <- c('IL1A','IL1B','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL2','IL2RA','IL2RB','IL2RG','IL3','IL3RA','IL4','IL4I1','IL4R','IL5','IL5RA','IL6','IL6R','IL6ST','IL7','IL7R','IL9','IL10','IL10RA','IL11','IL12A','IL12B','IL12RB1','IL12RB2','IL13','IL13RA2','IL15','IL15RA','IL16','IL17A','IL17B','IL17C','IL17D','IL17F','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL18BP','IL18R1','IL18RAP','IL19','IL20','IL20RA','IL20RB','IL21','IL21R','IL22','IL22RA2','IL23A','IL24','IL25','IL26','IL27','IL27RA','IL31','IL31RA','IL33','IL34','IL36A','IL36B','IL36G','IL37','ILDR1','ILDR2','ILF2','ILF3','ILK','ILKAP','ILVBL')
	PlotMarkerSeries(seuratObj, reductions, cytokines, 'Cytokines/Receptors')

	# KLRC2 = ENSMMUG00000050862
	klrs <- c('KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRF2', 'KLRG1', 'KLRG2', 'KLRC2', 'KLRC3', 'ENSMMUG00000050862')
	PlotMarkerSeries(seuratObj, reductions, klrs, 'KLRs')

	PlotMarkerSet(seuratObj, reductions, 'Resident Memory', c('ITGAE', 'ITGB7', 'CD69', 'CXCR6'))

	#chemokines
	chemokines <- c('CCL1','CCL11','CCL13','CCL16','CCL17','CCL18','CCL19','CCL2','CCL20','CCL21','CCL22','CCL24','CCL25','CCL26','CCL27','CCL28','CCL5','CCL7','CCL8')
	chemokines <- c(chemokines, c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CCRL2'))
	chemokines <- c(chemokines, c('CXCL1','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17','CXCL5','CXCL6','CXCL8','CXCL9','CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','XCR1'))

	PlotMarkerSeries(seuratObj, reductions, chemokines, 'Chemokines/Receptors')

	PlotMarkerSeries(seuratObj, reductions, c('MKI67'), 'Cell Proliferation')

	#PlotMarkerSeries(seuratObj, reductions, c('EPCAM'), 'Epithelial Cells')

	# LOC710951 = TRAC
	# LOC114677140 = TRBC1
	# LOC711031 = TRDC
	# LOC720538 = TRGC1
	# LOC705095 = TRGC2
	PlotMarkerSeries(seuratObj, reductions, c('LOC710951', 'LOC114677140', 'LOC711031', 'LOC720538', 'LOC705095'), 'TCR Constant Region')
}

PlotMarkerSeries <- function(seuratObj, reductions, features, title, setSize = 4) {
	featuresToPlot <- unique(intersect(features, row.names(seuratObj)))
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * 4) + 1
		end <- min((start + 3), length(featuresToPlot))
		genes <- featuresToPlot[start:end]


		PlotMarkerSet(seuratObj, reductions, title, genes)
	}
}

.RemoveUnchangedOrZero <- function(seuratObj, reduction, features) {
	ret <- c()
	#Remove zeros or unchanged:
	dims <- paste0(Key(object = seuratObj[[reduction]]), c(1,2))
	data <- FetchData(object = seuratObj, vars = c(dims, features), cells = colnames(x = seuratObj))
	for (feature in features) {
		if (sum(data[,feature] > 0) > 1 && length(unique(data[, feature])) > 1) {
			ret <- c(ret, feature)
		}
	}

	return(ret)
}

#' @import Seurat
#' @import patchwork
PlotMarkerSet <- function(seuratObj, reductions, title, features) {
	P <- NULL
	missingFeats <- c()
	for (reduction in reductions) {
		featuresToPlot <- intersect(features, row.names(seuratObj))
		featuresToPlot <- .RemoveUnchangedOrZero(seuratObj, reduction, featuresToPlot)

		if (length(features) != length(featuresToPlot)){
			missing <- features[!(features %in% featuresToPlot)]
			missingFeats <- unique(c(missingFeats, missing))
		}

		if (length(featuresToPlot) == 0){
			print('None of the requested features were present, skipping')
			next
		}

		P1 <- FeaturePlot(seuratObj, features = featuresToPlot, reduction = reduction, min.cutoff = 'q05', max.cutoff = 'q95')
		if (all(is.null(P))) {
			P <- P1
		} else {
			P <- P | P1
		}
	}

	if (length(missingFeats) > 0) {
		print(paste0('The following features were requested, but not present: ', paste0(missingFeats, collapse = ',')))
	}

	if (!all(is.null(P))) {
		P <- P + patchwork::plot_annotation(title = title)
		print(P)
	}
}
