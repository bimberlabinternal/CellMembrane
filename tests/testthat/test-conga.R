context("scRNAseq")

test_that("RunConga works", {
  if (!reticulate::py_module_available('conga')) {
    print('conga module not found, debugging:')
    print(reticulate::py_list_packages())
    if ('conga' %in% reticulate::py_list_packages()$package) {
      tryCatch({
        reticulate::import('conga')
      }, error = function(e){
        print("Error with reticulate::import('conga')")
        print(conditionMessage(e))
        traceback()
      })
    }

    warning('The python conga module has not been installed!')
    return()
  }

  #read data
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))
  tcr_df <- read.table("../testdata/tcr_df.csv", header = T, sep = ",") #this is spoofed TCR data from 438-21

  # limit to cells present in the TCR data:
  hasBothChains <- intersect(tcr_df$barcode[tcr_df$chain == 'TRA' & tcr_df$is_cell == 'true'], tcr_df$barcode[tcr_df$chain == 'TRB' & tcr_df$is_cell == 'true'])
  seuratObj <- subset(seuratObj, cells = hasBothChains)
  print(seuratObj)
  
  #The tcr_df was originally downloaded via Rdiscvr, but altered like so:
  #iterate through the colnames of the seurat object, replacing the values of tcr_df$barcodes and then drop the rest (the goal is to fake the tcr data using the gex barcodes)
  #i = 0
  #for(cell_barcode in unique(colnames(seuratObj))){
  #  barcode_to_replace <- unique(tcr_df$barcode)[i]
  #  tcr_df[tcr_df$barcode == barcode_to_replace,"barcode"] <- cell_barcode
  #  i=i+1
  #}
  #tcr_df <- tcr_df[tcr_df$barcode %in% colnames(seuratObj), ]
  #write.table(tcr_df, "../testdata/tcr_df.csv", col.names = T, sep = ",")

  congaSeuratObj <- RunCoNGA(seuratObj = seuratObj,
                             assayName = "RNA",
                             tcrClonesFile = "../testdata/tcr_df.csv",
                             seuratToCongaDir = "../testdata/tmpoutput",
                             organism = "rhesus",
                             runCongaOutputFilePrefix = "conga_output",
                             gexDatatype = "10x_h5",
                             runCongaOutfilePrefixForQcPlots = "qc_plots",
                             runCongaOutputDirectory = "../testdata/tmpoutput", 
                             congaMetadataPrefix = "conga_")
  #test that the GEX file exists (i.e that SeuratToConga worked).
  testthat::expect_true(file.exists("../testdata/tmpoutput/GEX.h5"))
  #test that conga ran successfully
  testthat::expect_true(file.exists("../testdata/tmpoutput/conga_output_results_summary.html"))
  #test that clustering worked and was appended to the seurat object.
  testthat::expect_true(1 %in% congaSeuratObj@meta.data[,"conga_clusters_gex"])
  unlink("./tmpoutput", recursive = TRUE)
})

test_that("CalculateTcrDiversity works", {
  if (!reticulate::py_module_available('tcrdist')) {
    print('tcrdist3 module not found, debugging:')
    print(reticulate::py_list_packages())
    if ('tcrdist' %in% reticulate::py_list_packages()$package) {
      tryCatch({
        reticulate::import('tcrdist')
      }, error = function(e){
        print("Error with reticulate::import('tcrdist')")
        print(conditionMessage(e))
        traceback()
      })
    }
    
    warning('The python tcrdist3 module has not been installed!')
    expect_true(reticulate::py_module_available('conga'))
  }

  dat <- read.table("../testdata/clones_file.txt", sep = '\t', header = TRUE)
  dat <- dat[c('clone_id', 'clone_size', 'va_gene', 'vb_gene', 'cdr3a', 'cdr3b')]
  names(dat) <- c('sampleId', 'clone_size', 'v_a_gene', 'v_b_gene', 'cdr3_a_aa', 'cdr3_b_aa')
  dat$sampleId <- unlist(sapply(dat$sampleId, function(x){
    return(unlist(strsplit(x = x, split = '_'))[1])
  }))
  dat <- dat[!is.na(dat$v_a_gene) & !is.na(dat$v_b_gene) & !is.na(dat$cdr3_a_aa) & !is.na(dat$cdr3_b_aa),]
  df <- CalculateTcrDiversity(dat,
                order1 = 1,
                order2 = 200)

  testthat::expect_equal(nrow(df), 199)
})

test_that("QuantifyTcrClones  works", {
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(readRDS("../testdata/seuratOutput.rds")))
  seuratObj <- QuantifyTcrClones(seuratObj, "../testdata/tcr_df.csv", groupingFields = 'ClusterNames_0.2')

  expect_equal(length(unique(seuratObj$clonotypeID)), 400)
  expect_equal(length(unique(seuratObj$cloneSize)), 7)
  expect_equal(max(seuratObj$cloneProportion, na.rm = TRUE), 0.151, tolerance = 0.001)
})

test_that("CalculateTcrRepertioreStats  works", {
  df <- read.csv("../testdata/RepertoireData.csv")
  dfout <- CalculateTcrRepertioreStats(df, "cDNA_ID")
  
  expect_equal(nrow(dfout), 38)
  expect_equal(dfout[dfout$MetricName == "TRA_0.2",]$Value, 14)
  expect_equal(dfout[dfout$MetricName == "TRA_TopFrac",]$Value, 0.0423, tolerance = 0.001)
})

