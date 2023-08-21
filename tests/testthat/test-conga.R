context("scRNAseq")

test_that("RunConga works", {
  #temp solution, delete before pushing
  setwd("/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tests/testdata/")
  #read data
  seuratObj <- readRDS("../testdata/seuratOutput.rds")
  tcr_df <- read.table("../testdata/tcr_df.csv", header = T, sep = ",") #this is spoofed TCR data from 438-21
  
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
  #create a temporary directory to store the output from RunCoNGA
  seuratToCongaDirectory <- "../testdata/tmpoutput"
  SeuratToCoNGA(seuratObj, 
                tcr_file_from_rdiscvr = "../testdata/tcr_df.csv", 
                seuratToCongaDirectory = seuratToCongaDirectory, 
                overwrite = T)
  #test that the GEX file exists (i.e that SeuratToConga worked).
  testthat::expect_true(file.exists(paste0(seuratToCongaDirectory, "/GEX.h5")))
  #create a temporary directory to store the output from RunCoNGA
  setwd(file.path("../testdata", "tmpoutput"))
  RunCoNGA(seuratToCongaDirectory = seuratToCongaDirectory, 
           gex_datatype = "10x_h5",
           organism = "rhesus", 
           outfile_prefix = "testing_conga",
           outfile_prefix_for_qc_plots = "plots",
           clones_file = "/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tests/testdata/tmpoutput/tmp_clonesfile.txt", 
           working_directory = getwd())
  
})
  