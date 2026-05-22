
#' @title Run StarCAT
#'
#' @description Runs starCAT on a Seurat object and appends per-cell scores and program usage to its metadata.
#' @param seuratObj The Seurat object containing the data to be run through starCAT.
#' @param reference The starCAT reference to use. Either a built-in reference name (e.g. 'TCAT.V1', 'MYELOID.GLIOMA.V1', 'BONEMARROW.CD34POS.HSPC.V1') or a path to a custom reference .tsv/.txt file.
#' @param outputDirectory The directory where starCAT outputs (<name>.rf_usage_normalized.txt and, when the reference ships add-on scores, <name>.scores.txt) will be written.
#' @param assayName The name of the assay containing the counts matrix that will be passed to starCAT.
#' @param outputPrefix Prefix for starCAT output files.
#' @param starcatMetadataPrefix A prefix prepended to columns added to seuratObj@meta.data. Scores are added as <prefix>score_*, usage columns as <prefix>usage_*.
#' @param addUsageToMetadata If TRUE, append the cNMF program usage matrix to seuratObj@meta.data. If FALSE, usage is left as a file on disk in outputDirectory.
#' @param cacheDirectory Directory where starCAT downloads / caches built-in reference files. Defaults to <outputDirectory>/cache. Use a shared persistent path to avoid re-downloading the reference across runs.
#' @param cleanUpIntermediateFiles If TRUE (default), delete everything RunStarCAT wrote once results are on the Seurat object. Set FALSE to keep the starCAT outputs and reference cache in outputDirectory.
#' @return The input Seurat object with starCAT scores (and optionally program usage) appended to meta.data.
#' @examples
#' \dontrun{
#'   seuratObj <- RunStarCAT(seuratObj = seuratObj,
#'                           reference = 'TCAT.V1')
#' }
#' @export
RunStarCAT <- function(seuratObj,
                       reference                = "TCAT.V1",
                       outputDirectory          = "./starcat_output",
                       assayName                = "RNA",
                       outputPrefix             = "starcat_run",
                       starcatMetadataPrefix    = "starcat_",
                       addUsageToMetadata       = TRUE,
                       cacheDirectory           = NULL,
                       cleanUpIntermediateFiles = TRUE) {

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('starcat')) {
    stop('The starcat python package has not been installed! If you believe it has been installed, run reticulate::import("starcat") to get more information and debug')
  }

  if (is.null(seuratObj) || !inherits(seuratObj, 'Seurat')) {
    stop('seuratObj must be a Seurat object.')
  }

  if (!is.character(reference) || length(reference) != 1L || !nzchar(reference)) {
    stop("`reference` must be a single non-empty string (e.g. 'TCAT.V1' or a path to a reference TSV).")
  }

  if (grepl("/$", outputDirectory)) {
    outputDirectory <- gsub("/$", "", outputDirectory)
  }

  if (!dir.exists(outputDirectory)) {
    dir.create(outputDirectory, recursive = TRUE)
  }

  outputDirectory <- gsub(R.utils::getAbsolutePath(outputDirectory), pattern = '\\\\', replacement = '/')

  if (is.null(cacheDirectory)) {
    cacheDirectory <- paste0(outputDirectory, '/cache')
  }
  if (!dir.exists(cacheDirectory)) {
    dir.create(cacheDirectory, recursive = TRUE)
  }
  cacheDirectory <- gsub(R.utils::getAbsolutePath(cacheDirectory), pattern = '\\\\', replacement = '/')

  gexOutFile <- R.utils::getAbsolutePath(paste0(outputDirectory, "/counts.h5"), mustWork = FALSE)
  gexOutFile <- gsub(gexOutFile, pattern = '\\\\', replacement = '/')

  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'counts'),
                               path = gexOutFile,
                               overwrite = TRUE)

  pyTemplate <- readr::read_file(system.file("scripts/run_StarCAT.py", package = "CellMembrane"))
  if (!nzchar(pyTemplate)) {
    stop('Could not locate inst/scripts/run_StarCAT.py from the installed CellMembrane package. Reinstall the package.')
  }

  script <- tempfile(fileext = ".py")
  readr::write_file(pyTemplate, script)

  call <- paste0("run_StarCAT(gex_datafile='", gexOutFile,
                 "', reference='", reference,
                 "', output_dir='", outputDirectory,
                 "', name='", outputPrefix,
                 "', cachedir='", cacheDirectory,
                 "')")

  readr::write_file('\n', script, append = TRUE)
  readr::write_file(call, script, append = TRUE)

  status <- system2(reticulate::py_exe(), script)
  if (is.integer(status) && status != 0L) {
    stop(sprintf("starcat python script exited with status %d.", status))
  }

  usageFile  <- paste0(outputDirectory, "/", outputPrefix, ".rf_usage_normalized.txt")
  scoresFile <- paste0(outputDirectory, "/", outputPrefix, ".scores.txt")

  if (!file.exists(usageFile)) {
    stop("starCAT finished but the expected usage output was not found at ", usageFile)
  }

  usage <- utils::read.table(usageFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

  if (!setequal(rownames(usage), colnames(seuratObj))) {
    stop("The cell barcodes in the starCAT usage output do not match the cells in the Seurat object.")
  }

  scores <- NULL
  if (file.exists(scoresFile)) {
    scores <- utils::read.table(scoresFile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    if (!setequal(rownames(scores), colnames(seuratObj))) {
      stop("The cell barcodes in the starCAT scores output do not match the cells in the Seurat object.")
    }
    colnames(scores) <- paste0(starcatMetadataPrefix, "score_", colnames(scores))
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = scores)
  } else {
    warning("No scores file found at ", scoresFile, ". The selected reference may not ship add-on score definitions; only usage will be attached.")
  }

  if (addUsageToMetadata) {
    colnames(usage) <- paste0(starcatMetadataPrefix, "usage_", colnames(usage))
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = usage)
  }

  if (cleanUpIntermediateFiles) {
    unlink(c(gexOutFile, script, usageFile, scoresFile))
    unlink(cacheDirectory, recursive = TRUE)
    if (length(list.files(outputDirectory, all.files = TRUE, no.. = TRUE)) == 0L) {
      unlink(outputDirectory, recursive = TRUE)
    }
  }

  return(seuratObj)
}
