context("scRNAseq")

test_that("RunConga works", {
 #this technically already has the TCR data appended, but that should be fine.
 seuratObj <- readRDS("../testdata/seuratOutputWithTCR.rds")
 SeuratToCoNGA(seuratObj, clonotypesFile = "../tcr_df.csv", 
               outputDir = "./tmpoutput")
 #outfile_prefix
 RunCoNGA(variable_features_file = "/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tmpoutput/varfeats.csv", 
          gex_datafile = "/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tmpoutput/GEX.h5", 
          gex_datatype = "10x_h5",
          tcr_datafile = "/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tests/testdata/tcr_df.csv", 
          organism = "rhesus", 
          outfile_prefix = "testing_conga",
          outfile_prefix_for_qc_plots = "plots",
          clones_file = "/Users/mcelfreg/scRNASeq/bimber-lab/CellMembrane/tests/testdata/tmp_clonesfile.txt")
 
})