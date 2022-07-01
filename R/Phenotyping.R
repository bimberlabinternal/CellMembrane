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
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CD8A', 'CD8B', 'ENSMMUG00000003532', 'CD4', 'IL7R', 'CD3D', 'CD3E','CD3G'), 'CD8/CD4 Markers')

	#Eff v. Mem:
	#IL7R = CD127
	#IL2RA = CD25
	#PTPRC = CD45
	#SELL = CD62-L / CD-197
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CCR7', 'SELL', 'GZMB', 'CCR5', 'IL2RA', 'PTPRC', 'IL7R', 'CTLA4', 'FAS', 'CD28', 'CD27', 'ITGA4', 'ITGB7', 'ITGB1'), 'Effector vs. Memory')

	#CD8 Activation
	# XCL1 = ENSMMUG00000060218
	# CCL4 = ENSMMUG00000008111
	# LOC100423131 = XCL1, ENSMMUG00000013779, Lymphotactin
	# LOC100430627 = CCL4, ENSMMUG00000008111
	# LOC100426632 = C-C motif chemokine 4
	# LOC100426537 = C-C motif chemokine 3-like
	# LOC114673087 = C-C motif chemokine 3-like 1
	# LOC100429751 = C-X-C chemokine receptor type 1-like
	# LOC701946 = C-X-C motif chemokine 5
	# LOC703222 = C-X-C motif chemokine 5
	# LOC100423954 = C-C motif chemokine 3-like 1
	# CD154 = CD40L
	# LAMP1 = CD107a, cytotoxic capacity
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('IFNG', 'CD69', 'TNF', 'NFKBID', 'LTB', 'TNFRSF9', 'CCL4L1', 'NR4A3', 'TNFSF14', 'CD82', 'PIGT', 'IRF8', 'IRF4', 'IRF2', 'RGCC', 'PD1', 'PDCD1', 'TNFAIP3', 'LOC100423131', 'LOC100430627', 'ENSMMUG00000013779', 'XCL1', 'ENSMMUG00000060218', 'CCL4', 'ENSMMUG00000008111', 'CCL3', 'PLEK', 'NR4A2', 'LOC100426537', 'LOC114673087', 'KLF10', 'GADD45B', 'CD154', 'LAMP1', 'LOC703222', 'LOC701946', 'LOC100429751', 'LOC100423954'), 'CD8 Activation Markers')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'), 'Cytotoxicity')

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'B-cell Markers', c('MS4A1', 'CD79A', 'CD74', 'DRA'))

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Monocyte', c('LYZ', 'CST3', 'S100A6', 'VIM'))

	# ITGB3 = CD61
	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Platelet/MK', c('ITGB3', 'PPBP'))

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Stemness', c('CD34'))

	# MRC1 = CD206
	# ITGAM = CD11b, AM=CD11b-, Non-AM=CD11b+
	# Non-AMs: CD16+/CD206-/HLA-DR+/CD11b+
	# M1/M2: CCR7 (M1), CD163 (M2), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7050058/
	# CD163, CD68 = macrophage markers
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CD14', 'FCGR3A', 'S100A8', 'S100A6', 'MARCO', 'MRC1', 'CD163', 'CHIT1', 'APOBEC3A', 'ITGAM', 'HLA-DRB1', 'CCR7', 'CD68'), 'Myeloid')

	# IL3RA = CD123 / pDC
	# CLEC4C = CD303 / pDC
	# CD1c = mDC
	# THBD = CD141 / mDC
	# CD80, CD86 = co-stimulatory. low expression = tolerogenic
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CD14', 'FCGR3A', 'IL3RA', 'CLEC4C', 'CD1c', 'THBD', 'CD80', 'CD86', 'TGFB1'), 'DCs')

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Regulatory Markers', c('CD4', 'FOXP3', 'IL2RA', 'NR4A1'))

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Treg17', c('RORA', 'RORB', 'RORC', 'IL4', 'STAT3'))

	# ZBTB16 = PLZF
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('HAVCR2'), 'Th1')

	# ZBTB16 = PLZF
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('ZBTB16', 'DPP4'), 'MAIT')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('TBX21', 'GATA3', 'RORC', 'FOXP3', 'BCL6', 'EOMES', 'TOX', 'GATA2', 'TCF7', 'KLF2', 'NR4A1', 'LEF1', 'PRDM1', 'ID2', 'ID3'), 'Transcription Factors')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('TIGIT', 'CTLA4', 'BTLA', 'PDCD1', 'CD274'), 'Inhibitory Markers')

	# https://www.frontiersin.org/articles/10.3389/fimmu.2016.00076/full#:~:text=The%20three%20main%20pathways%20activated,activation%20(1%2C%202).

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'NFATs', c('STAT1', 'STAT2', 'STAT3', 'STAT4'))

	PlotMarkerSet(seuratObj, reductions = reductions, features =  'NFATs', c('NFATC1', 'NFATC2', 'NFAT3', 'NFATC4', 'NFAT5'))

	#DAP10/12
	# LOC707555 = ZAP70
	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Signaling', c('HCST', 'TYROBP', 'SYK', 'ZAP70', 'LOC707555', 'ITK'))

	#LILR/KIR:
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('LILRA5','LILRA6','LILRB4','LILRB5','KIR2DL4','KIR3DX1', 'MAMU-KIR', 'KIR2DL4', 'KIR3DL2'), 'LILR/KIR')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('FCGR1A','FCGR2A','FCGR2B','FCGR3', 'FCGR3A'), 'FCGR')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  GetGeneSet('EffectorCytokines'), 'Effector Cytokines')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  GetGeneSet('ExhaustionOrInhibitory'), 'Exhaustion Or Inhibitory')

	#Cytokines
	cytokines <- c('IL1A','IL1B','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL2','IL2RA','IL2RB','IL2RG','IL3','IL3RA','IL4','IL4I1','IL4R','IL5','IL5RA','IL6','IL6R','IL6ST','IL7','IL7R','IL9','IL10','IL10RA','IL11','IL12A','IL12B','IL12RB1','IL12RB2','IL13','IL13RA2','IL15','IL15RA','IL16','IL17A','IL17B','IL17C','IL17D','IL17F','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL18BP','IL18R1','IL18RAP','IL19','IL20','IL20RA','IL20RB','IL21','IL21R','IL22','IL22RA2','IL23A','IL24','IL25','IL26','IL27','IL27RA','IL31','IL31RA','IL33','IL34','IL36A','IL36B','IL36G','IL37','ILDR1','ILDR2','ILF2','ILF3','ILK','ILKAP','ILVBL')
	# LOC114675338	interleukin-17 receptor B-like
	# LOC114675718	interleukin-3 receptor subunit alpha-like
	# LOC100426853	interleukin-9 receptor-like
	cytokines <- c(cytokines, c('LOC114675338', 'LOC114675718', 'LOC100426853'))
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  cytokines, 'Cytokines/Receptors')

	# KLRC2 = ENSMMUG00000050862
	klrs <- c('KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRF2', 'KLRG1', 'KLRG2', 'KLRC2', 'KLRC3', 'KLRK1', 'ENSMMUG00000050862')
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  klrs, 'KLRs')

	# NCR1 = NKp46
	# NCR2 = NKp44
	# NCR3 = NKp30
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('NCR1', 'NCR2', 'NCR3'), 'NCRs')

	# ITGAE = CD103
	PlotMarkerSet(seuratObj, reductions = reductions, features =  'Resident Memory', c('ITGAE', 'ITGB7', 'CD69', 'CXCR6'))

	#chemokines
	chemokines <- c('CCL1','CCL11','CCL13','CCL16','CCL17','CCL18','CCL19','CCL2','CCL20','CCL21','CCL22','CCL24','CCL25','CCL26','CCL27','CCL28','CCL5','CCL7','CCL8')
	chemokines <- c(chemokines, c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CCRL2'))
	chemokines <- c(chemokines, c('CXCL1','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17','CXCL5','CXCL6','CXCL8','CXCL9','CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','XCR1'))

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  chemokines, 'Chemokines/Receptors')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('MKI67', 'TOP2A'), 'Cell Proliferation')

	#PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('EPCAM'), 'Epithelial Cells')

	# LOC710951 = TRAC
	# LOC114677140 = TRBC1
	# LOC711031 = TRDC
	# LOC720538 = TRGC1
	# LOC705095 = TRGC2
	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('LOC710951', 'LOC114677140', 'LOC711031', 'LOC720538', 'LOC705095'), 'TCR Constant Region')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('MMP2','COL1A1','COL1A2','COL5A1','LUM','PDGFRA'), 'Stromal')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CCR9','LILRA4'), 'pDC')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CDH1','FLT1'), 'Epithelial')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('LYZ','CSF1R','MSR1','MAFB','CD300E'), 'MoMacDC')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('CSF3R','FCGR3B'), 'Neutrophils')

	PlotMarkerSeries(seuratObj, reductions = reductions, features =  c('HBB','HBA2','HBA1'), 'Erythrocyte')

}


#' @title PlotMarkerSeries
#'
#' @description Iteratively plots a set of markers
#' @param seuratObj The seurat object
#' @param features A vector of feature names
#' @param reductions The reductions to plot
#' @param title An optional title of this plot series
#' @param setSize The maximum number of features to include per FeaturePlot
#' @export
#' @import Seurat
PlotMarkerSeries <- function(seuratObj, features, reductions = c('umap'), title = NULL, setSize = 4) {
	featuresToPlot <- unique(intersect(features, row.names(seuratObj)))
	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * 4) + 1
		end <- min((start + 3), length(featuresToPlot))
		genes <- featuresToPlot[start:end]


		PlotMarkerSet(seuratObj, reductions = reductions, title = title, features = genes)
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
		if (!is.null(title)) {
			P <- P + patchwork::plot_annotation(title = title)
		}

		print(P)
	}
}


pkg.env$GENE_SETS <- list()

.RegisterGeneSet <- function(name, genes) {
	if (name %in% names(pkg.env$GENE_SETS)) {
		stop(paste0('Name has already been registered: ', name))
	}

	pkg.env$GENE_SETS[[name]] <- genes
}

#' @title GetGeneSet
#'
#' @description Returns a vector with the set of genes registered under the provided name.
#' @param name The name of the gene set
#' @export
#' @import Seurat
GetGeneSet <- function(name) {
	if (! (name %in% names(pkg.env$GENE_SETS))) {
		warning(paste0('Unknown gene set: ', name))
	}

	return(pkg.env$GENE_SETS[[name]])
}

.RegisterGeneSet('MMul10TcrGenes', c('LOC705095','LOC703029','LOC696306','LOC710951','LOC106999340','LOC106996262','LOC106999345','LOC106997707','LOC106997706','LOC710149','LOC700771','LOC699427','LOC711871','LOC709081','LOC698785','LOC114677052','LOC114676933','LOC106999353','LOC106999351','LOC106999350','LOC106999349','LOC106999348','LOC106999347','LOC106999346','LOC106999343','LOC106999341','LOC106999339','LOC106999337','LOC106999336','LOC106999335','LOC106999312','LOC106997705','LOC106997704','LOC106997703','LOC106997702','LOC106997697','LOC106997453','LOC106997452','LOC106997451','LOC106995765','LOC106992460','LOC106992446','LOC106992434','LOC106992433','LOC720538','LOC720456','LOC716949','LOC716866','LOC711537','LOC711386','LOC711194','LOC711141','LOC711066','LOC711031','LOC710821','LOC710627','LOC710455','LOC710361','LOC710183','LOC710093','LOC709531','LOC708581','LOC708328','LOC704883','LOC703153','LOC702904','LOC702550','LOC702113','LOC701992','LOC701875','LOC701745','LOC701395','LOC701262','LOC701152','LOC700224','LOC700154','LOC700105','LOC699912','LOC699790','LOC699543','LOC699298','LOC699162','LOC698913','LOC698543','LOC698289','LOC698161','LOC697792','LOC697466','LOC697234','LOC697054','LOC696752','LOC696684','LOC696557','LOC696075','LOC695943','LOC114679533','LOC114679531','LOC114677140','LOC114677139','LOC114677137','LOC114677136','LOC114677055','LOC114677054','LOC114677050','LOC114677049','LOC114677047','LOC114675766'))
.RegisterGeneSet('CD8_Activation.1', c('CCL4L1','MIR155HG','RGCC','NFKBIA','IFNG','NR4A3','TNFSF14','CCL3'))
.RegisterGeneSet('Cytotoxicity', c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'))
.RegisterGeneSet('EffectorCytokines', c('GZMA','GZMB','GZMH','GZMK', 'GZMM', 'PRF1', 'NKG7', 'GNLY', 'IFNG', 'FASLG', 'TNF', 'IL17A', 'IL2'))
.RegisterGeneSet('ExhaustionOrInhibitory', c('PDCD1', 'TIGIT', 'HAVCR2', 'LAG3', 'CTLA4', 'VTCN1', 'CD244', 'KLRG1', 'TNFRSF14', 'BTLA', 'CD160'))

.RegisterGeneSet('Myelocytes', c("OLFM4", "LTF", "CAMP", "LCN2"))
.RegisterGeneSet('Pro_Myelocytes', c("BPI", "MPO", "ELANE", "RTD1A", "RTD1B", "DEFA1B", "AZU1"))
.RegisterGeneSet("MMul10_MHC", c("MAMU-A", "MAMU-A3", "MAMU-AG", "LOC719250", "LOC100426197", "LOC714466", "LOC694372", "LOC114677644", "LOC114675360", "LOC106997902", "LOC698738", "LOC106996077", "LOC106992378", "LOC100424348", "LOC106995519", "LOC106992468", "LOC114676051", "LOC106996627", "LOC106997893", "LOC106997885", "LOC114675646", "LOC114669738", "LOC720132", "LOC719702", "LOC106996676", "LOC714964", "LOC114669810", "LOC114675357", "LOC100428435", "LOC100429195", "LOC699987", "LOC114677642", "LOC723473", "LOC106995461", "LOC106995452"))
.RegisterGeneSet("MMul10_KIR", c("KIR3DL0", "KIR3DL2", "KIR2DS4", "KIR2DL4", "KIR3DS05", "KIR3DL12", "KIR3DL1", "KIR3DH", "KIR3DH5", "KIR3DL21", "KIR3DL11", "MAMU-KIR", "LOC106994859", "LOC114669735", "LOC100424026", "LOC100125572"))

# Generated using: sort(rownames(seuratObj@assays$RNA)[grepl(rownames(seuratObj@assays$RNA), pattern = '^RP[SL]')])
ribosomalGenes <- c("RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36AL", "RPL37", "RPL37A", "RPL37A-1", "RPL38", "RPL39", "RPL39L", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPLP0", "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS19BP1", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27A-1", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2", "RPS6KC1", "RPS6KL1", "RPS7", "RPS8", "RPS9", "RPSA")

# These were identified by grep on the descriptions of all LOC* genes for ribosomal
ribosomalGenes <- c(ribosomalGenes, c("LOC717349", "LOC718596", "LOC705988", "LOC703571", "LOC100429589", "LOC114672439", "LOC711644", "LOC705840", "LOC114674534", "LOC697510", "LOC100426562", "LOC100430397", "LOC705494", "LOC718356", "LOC708673", "LOC708221", "LOC693728", "LOC701224", "LOC695949", "LOC100430766", "LOC708539", "LOC114670541", "LOC694687", "LOC696921", "LOC106999910", "LOC705621", "LOC706532", "LOC701369", "LOC106994001", "LOC114675582", "LOC706087", "LOC711318", "LOC695920", "LOC697723", "LOC700257", "LOC710331", "LOC700807", "LOC106992793", "LOC717846", "LOC694712", "LOC708815", "LOC702606", "LOC707204", "LOC106996615", "LOC711714", "LOC697734", "LOC100425632", "LOC709255", "LOC703037", "LOC710558", "LOC114670812", "LOC721996", "LOC716997", "LOC698732", "LOC702174", "LOC114679384", "LOC697431", "LOC704429", "LOC106992303", "LOC114676316", "LOC100430733", "LOC708154", "LOC114673092", "LOC697548", "LOC695895", "LOC708833", "LOC701909", "LOC715154", "LOC717289", "LOC106998255", "LOC701084", "LOC707369", "LOC697035", "LOC720167", "LOC700616", "LOC700401", "LOC106996191", "LOC106998695", "LOC707090", "LOC717450", "LOC693713", "LOC713382", "LOC716142", "LOC696119", "LOC714106", "LOC695421", "LOC709128", "LOC714782", "LOC693486", "LOC106992469", "LOC699246", "LOC708791", "LOC699630", "LOC697141", "LOC710726", "LOC107000432", "LOC100428100", "LOC718153", "LOC100427366", "LOC693343", "LOC114669718", "LOC713514", "LOC698010", "LOC715166", "LOC695327", "LOC702713", "LOC704640", "LOC114674200", "LOC700134", "LOC715041", "LOC722876", "LOC710644", "LOC705596", "LOC702961", "LOC704054", "LOC715668", "LOC717801", "LOC707085", "LOC114679603", "LOC106999518", "LOC718267", "LOC705327", "LOC106992281", "LOC695218", "LOC702957", "LOC708147", "LOC719947", "LOC701392", "LOC699166", "LOC106999718", "LOC702760", "LOC702328", "LOC704267", "LOC719289", "LOC106997304", "LOC100423207", "LOC100428069", "LOC106999080", "LOC100427696", "LOC712160", "LOC712267", "LOC100425072", "LOC100427254", "LOC696068", "LOC694967", "LOC708492", "LOC710583", "LOC717763", "LOC693820", "LOC114678655", "LOC114678366", "LOC114678384", "LOC114671497", "LOC114675033", "LOC697719", "LOC694071", "LOC695523", "LOC703455", "LOC715814", "LOC106998837", "LOC713060", "LOC693904", "LOC711043", "LOC709241", "LOC693584", "LOC698301", "LOC703631", "LOC701292", "LOC706910", "LOC712630", "LOC714801", "LOC717838", "LOC698099", "LOC694067", "LOC706082", "LOC702514", "LOC711324", "LOC703616", "LOC704396", "LOC712274", "LOC698273", "LOC718556", "LOC715043", "LOC696015", "LOC698297", "LOC700262", "LOC695670", "LOC698383", "LOC708603", "LOC701421", "LOC697065", "LOC114677128", "LOC714656", "LOC106996346", "LOC703683", "LOC701691", "LOC698792", "LOC107000151", "LOC697561", "LOC710477", "LOC100427278", "LOC720748", "LOC114672187", "LOC708140", "LOC698768", "LOC716735", "LOC114673268", "LOC114674609", "LOC114674676", "LOC114674882", "LOC114675055", "LOC114675456", "LOC114674615", "LOC114674626", "LOC114674635", "LOC114674636", "LOC114674637", "LOC114674638", "LOC114674640", "LOC114674642", "LOC114674650", "LOC114674666", "LOC114674683", "LOC114673870", "LOC114675629", "LOC114675635", "LOC114675637", "LOC114675638", "LOC114675641", "LOC114675642", "LOC114675643", "LOC114675636", "LOC114675644", "LOC114675630", "LOC114675633", "LOC114675631", "LOC114675632", "LOC114675739", "LOC114675743", "LOC114675745", "LOC114675746", "LOC114675747", "LOC114675740", "LOC114675741", "LOC114675856", "LOC114675854", "LOC114675853", "LOC114675857", "LOC114675858", "LOC114675859", "LOC114675860", "LOC114675913", "LOC114675919", "LOC114675922", "LOC114675923", "LOC114675924", "LOC114675925", "LOC114675914", "LOC114675915", "LOC114675916", "LOC114675918", "LOC114675141", "LOC718979", "LOC106993215", "LOC106995250", "LOC709583", "LOC114679608", "LOC717779", "LOC707414", "LOC709376", "LOC710038", "LOC706793", "LOC106999723", "LOC696958", "LOC106993579", "LOC106995074", "LOC694799", "LOC106995394", "LOC698530", "LOC107000213", "LOC705400", "LOC720019", "LOC698143", "LOC697585", "LOC713714", "LOC702620", "LOC715037", "LOC699544", "LOC100427498", "LOC699681", "LOC701302", "LOC700008", "LOC698713", "LOC702871", "LOC693856", "LOC701386", "LOC716126", "LOC696423", "LOC114669857", "LOC697476", "LOC100429843", "LOC696901", "LOC693272", "LOC704083", "LOC702723", "LOC703433", "LOC722545", "LOC696877", "LOC106998431", "LOC106999185", "LOC706492", "LOC699695", "LOC713984", "LOC703435", "LOC709778", "LOC694122", "LOC711449", "LOC701763", "LOC719242", "LOC715711", "LOC700872", "LOC717246", "LOC698182", "LOC701156", "LOC106996680", "LOC710859", "LOC106997591", "LOC707834", "LOC100428649", "LOC695340", "LOC718357", "LOC702542", "LOC718020", "LOC699687", "LOC711174", "LOC709462", "LOC712987", "LOC106995957", "LOC708800", "LOC717286", "LOC698999", "LOC696145", "LOC106999371", "LOC698377", "LOC709045", "LOC711760", "LOC114680233", "LOC107000909", "LOC106994077", "LOC701710", "LOC114674062", "LOC699398", "LOC693578", "LOC106992566", "LOC106997285", "LOC718737", "LOC706606", "LOC704012", "LOC702297", "LOC114680376", "LOC696154", "LOC721751", "LOC704238", "LOC709820", "LOC114676029", "LOC708995", "LOC710502", "LOC695715", "LOC698942", "LOC696134", "LOC696785", "LOC698684", "LOC695740", "LOC706177", "LOC717053", "LOC698197", "LOC714598", "LOC106998794", "LOC704510", "LOC698130", "LOC107000483", "LOC716320", "LOC106994271", "LOC708074", "LOC100425638", "LOC709201", "LOC707358", "LOC706712", "LOC701457", "LOC114676221", "LOC114677736", "LOC100423145", "LOC708782", "LOC714171", "LOC707437", "LOC701466", "LOC719029", "LOC697987", "LOC714120", "LOC708085", "LOC694937", "LOC106993116", "LOC706340", "LOC697214", "LOC106995063", "LOC699344", "LOC114676734", "LOC702875", "LOC695122", "LOC717674", "LOC712555", "LOC711491", "LOC702677", "LOC704365", "LOC708118", "LOC709678", "LOC106994904", "LOC703784", "LOC702847", "LOC709086", "LOC114674318", "LOC693848", "LOC719244", "LOC114680572", "LOC719770", "LOC710317", "LOC714620", "LOC703648", "LOC703365", "LOC706181", "LOC700425", "LOC708644", "LOC709650", "LOC705067", "LOC713497", "LOC106994132", "LOC703919", "LOC114671097", "LOC703853", "LOC114674108", "LOC114676998", "LOC712599", "LOC114678374", "LOC114678397", "LOC114679023", "LOC699392", "LOC114678167", "LOC694228", "LOC710665", "LOC107000401", "LOC100423661", "LOC106999480", "LOC106997918", "LOC107000809", "LOC114673683", "LOC106998280", "LOC106998813", "LOC100429752", "LOC100423423", "LOC106994361", "LOC106994914", "LOC705174", "LOC707072", "LOC700862", "LOC693472", "LOC106998156", "LOC721192", "LOC707117", "LOC705979", "LOC705677", "LOC710287", "LOC704345", "LOC699239", "LOC702116", "LOC706754", "LOC703797", "LOC710877", "LOC696151", "LOC702932", "LOC695743", "LOC693287", "LOC703024", "LOC698657", "LOC100429815", "LOC707342", "LOC695463", "LOC699805", "LOC722113", "LOC694734", "LOC106998116", "LOC716075", "LOC703563", "LOC713469", "LOC714185", "LOC703251", "LOC703283", "LOC698417", "LOC704062", "LOC698226", "LOC114677434", "LOC693480", "LOC114676944", "LOC702978", "LOC100427953", "LOC100427759", "LOC704256", "LOC700045", "LOC702525", "LOC710595", "LOC698574", "LOC707182", "LOC700578", "LOC698448", "LOC709444", "LOC698086", "LOC714855", "LOC708405", "LOC114675505", "LOC697178", "LOC703216", "LOC716530", "LOC716150", "LOC694093", "LOC718478", "LOC707861", "LOC106998326", "LOC695279", "LOC106995044", "LOC100423364", "LOC100423295", "LOC704394", "LOC703794", "LOC693573", "LOC699375", "LOC701250", "LOC706496", "LOC693947", "LOC713991", "LOC710642", "LOC109910387", "LOC109910386", "LOC106997745", "LOC714700", "LOC100426636", "LOC106998409"))
.RegisterGeneSet("MMul10_Ribosomal", ribosomalGenes)

# Identified by grep on the descriptions of all LOC* genes for 'mitochond'
mitochondrialGenes <- c("LOC703370", "LOC721815", "LOC720015", "LOC703345", "LOC708979", "LOC704978", "LOC717349", "LOC718596", "LOC705988", "LOC703571", "LOC100429589", "LOC114672439", "LOC711644", "LOC705840", "LOC697510", "LOC100426562", "LOC100430397", "LOC705494", "LOC718356", "LOC708673", "LOC708221", "LOC693728", "LOC701224", "LOC695949", "LOC716130", "LOC100430580", "LOC100423136", "LOC106998857", "LOC701346", "LOC695671", "LOC696244", "LOC100427705", "LOC114673694", "LOC106997368", "LOC107000220", "LOC696031", "LOC708562", "LOC716946", "LOC715928", "LOC713542", "LOC719688", "LOC699117", "LOC694182", "LOC696890", "LOC707114", "LOC704076", "LOC704440", "LOC709495", "LOC702581", "LOC711636", "LOC698075", "LOC711183", "LOC698715", "LOC694260", "LOC707476", "LOC700815", "LOC711417", "LOC706691", "LOC696144", "LOC695041", "LOC701131", "LOC114676443", "LOC712886", "LOC114669858", "LOC106996230", "LOC701964", "LOC705178", "LOC700087", "LOC695277", "LOC695732", "LOC114678889", "LOC106995294", "LOC100425230", "LOC100425700", "LOC697409", "LOC710274", "LOC707227", "LOC701661", "LOC696956", "LOC106996920", "LOC701236", "LOC106992283", "LOC717879", "LOC114674848", "LOC106998783", "LOC106998836", "LOC114678228", "LOC114673956", "LOC114677410", "LOC114678786", "LOC114679364", "LOC114679428", "LOC114680400", "LOC114671073", "LOC100423866", "LOC106993150", "LOC106996000", "LOC107000216", "LOC699137", "LOC106996051", "LOC713467", "LOC114675539", "LOC693843", "LOC697678", "LOC702427", "LOC114675509", "LOC698520", "LOC106997494", "LOC100425012", "LOC106999579", "LOC107000931", "LOC106997075", "LOC717923", "LOC106998133", "LOC712330", "LOC706686", "LOC715327", "LOC693799", "LOC693766", "LOC698522", "LOC114675319", "LOC718745", "LOC718390", "LOC114678805", "LOC706366", "LOC114674124", "LOC106996824", "LOC106994877", "LOC106992611", "LOC107000653", "LOC107000976", "LOC713028", "LOC693813", "LOC711899", "LOC718711", "LOC114669720", "LOC721639", "LOC698302", "LOC718712", "LOC698305", "LOC706710", "LOC705282", "LOC114670386", "LOC699779", "LOC100427337", "LOC100429934", "LOC704917", "LOC715730", "LOC705066", "LOC709894", "LOC711634", "LOC723069", "LOC705183", "LOC702309", "LOC701966", "LOC704336", "LOC697416", "LOC106993325", "LOC716730", "LOC717607", "LOC106997886", "LOC107000665", "LOC106999488", "LOC106998010", "LOC705476", "LOC114673096", "LOC696696", "LOC715179")
.RegisterGeneSet("MMul10_Mitochondrial", mitochondrialGenes)

igHeavyVariable <- c("LOC715358", "LOC114669771", "LOC106999620", "LOC106999619", "LOC114679735", "LOC114679700", "LOC114679702", "LOC114679730", "LOC114679696", "LOC106999621", "LOC114679699", "LOC721017", "LOC106995637", "LOC721043", "LOC721139", "LOC114679710", "LOC715165", "LOC720899", "LOC714894", "LOC106999623", "LOC106999617", "LOC714310", "LOC723533", "LOC715400", "LOC701407", "LOC720918", "LOC106999616", "LOC721119", "LOC114679704", "LOC715260", "LOC719378", "LOC106999631", "LOC721089", "LOC106999608", "LOC721132", "LOC720890", "LOC106999629", "LOC114679692", "LOC100430012", "LOC106999615", "LOC114669802", "LOC114679711", "LOC714939", "LOC715543", "LOC106999624", "LOC721104", "LOC114679713", "LOC721116", "LOC114679694", "LOC114679729", "LOC106999604", "LOC715594", "LOC106995900", "LOC106999610", "LOC106999609", "LOC721083", "LOC106999633", "LOC106996003", "LOC114679703", "LOC721002", "LOC106999628", "LOC114679693", "LOC722278", "LOC720985", "LOC720452", "LOC720904", "LOC114679728", "LOC114679707", "LOC720935", "LOC721108", "LOC720974", "LOC106999618", "LOC114679733", "LOC715733")
igKappaVariable <- c("LOC106993058", "LOC721353", "LOC703375", "LOC701600", "LOC701068", "LOC106993063", "LOC707007", "LOC709162", "LOC709703", "LOC709066", "LOC106993071", "LOC709793", "LOC707181", "LOC106993142", "LOC106992444", "LOC106993055", "LOC703139", "LOC114671809", "LOC106993143", "LOC723660", "LOC703610", "LOC106993056", "LOC114671797", "LOC708976", "LOC106993144", "LOC106998680", "LOC106993140", "LOC106993064", "LOC106993141", "LOC706869", "LOC707609", "LOC708249", "LOC709890", "LOC106993048", "LOC700696")
igLambdaVariable <-c("LOC107000570", "LOC106992425", "LOC701976", "LOC708547", "LOC107000578", "LOC114670549", "LOC107000585", "LOC107000586", "LOC107000581", "LOC114670550", "LOC107000576", "LOC107000577", "LOC702421", "LOC701772", "LOC107000579", "LOC701240")
.RegisterGeneSet('MMul10_Ig_Variable', unique(c(igHeavyVariable, igKappaVariable, igLambdaVariable)))

.RegisterGeneSet('VariableGenes_Exclusion.1', unique(c(
	GetGeneSet('MMul10TcrGenes'),
	GetGeneSet('MMul10_MHC'),
	GetGeneSet('MMul10_KIR'),
	GetGeneSet('MMul10_Ribosomal'),
	GetGeneSet('MMul10_Mitochondrial'),
	GetGeneSet('MMul10_Ig_Variable')
)))

#' @title GetMMul10TcrGenes
#'
#' @description Returns a vector with MMul10 gene IDs (NCBI build) for TCR genes.
#' @export
#' @import Seurat
GetMMul10TcrGenes <- function(){
	return(GetGeneSet('MMul10TcrGenes'))
}

#' @title GetMMul10IgGenes
#'
#' @description Returns a vector with MMul10 gene IDs (NCBI build) for TCR genes.
#' @export
#' @import Seurat
GetMMul10IgGenes <- function(){
	return(c(
		'LOC720839', #immunoglobulin heavy constant alpha 1-like
		'LOC708891', #immunoglobulin heavy constant gamma 1-like
		'LOC114679689', #immunoglobulin heavy constant epsilon-like
		'LOC114679691', #immunoglobulin heavy constant gamma 2-like
		'LOC114679690', #immunoglobulin heavy constant gamma 4-like
		'LOC710905', #immunoglobulin heavy constant gamma 4-like
		'LOC711872', # immunoglobulin heavy constant mu-like
		'LOC698810', #immunoglobulin iota chain
		'LOC114679695', #immunoglobulin mu heavy chain-like
		'LOC701504', # immunoglobulin kappa light chain
		'LOC107000555', # immunoglobulin lambda constant 6-like
		'LOC106996055', #immunoglobulin lambda-1 light chain-like
		'LOC708771', #immunoglobulin lambda constant 6-like
		'LOC106992418', #immunoglobulin lambda constant 6-like
		'LOC114670307', #immunoglobulin lambda constant 6-like
		'LOC106999340' # immunoglobulin kappa light chain
	))
}

#' @title ExpandGeneList
#'
#' @description Takes an input gene list and identifies any entries matching registered gene sets. Those will be expanded to the full gene list.
#' @param genes A vector of genes or gene set names
#' @param verbose Whether to log information about matches
#' @export
ExpandGeneList <- function(genes, verbose = TRUE) {
	genesMatchingSets <- genes[genes %in% names(pkg.env$GENE_SETS)]
	if (verbose && length(genesMatchingSets) > 0) {
		print(paste0('The following symbols match gene sets and will be expanded: ', paste0(genesMatchingSets, collapse = ',')))
	}

	ret <- genes[!genes %in% names(pkg.env$GENE_SETS)]
	for (geneSet in genesMatchingSets) {
		ret <- unique(c(ret, pkg.env$GENE_SETS[[geneSet]]))
	}

	return(ret)
}