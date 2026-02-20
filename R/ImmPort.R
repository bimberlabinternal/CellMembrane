# functions for working with ImmPort gene/pathway data.
# the gene-to-pathway map is stored as inst/extdata/immport_gene_pathway_map.rds
# This was last run on Feb 20, 2026. To update, see below to pull the latest GMT from ImmPort and rebuild the RDS using the following code: 
# library(GSEABase)
# raw_lines <- readLines("all_gene_lists.gmt") # gene list from immport (link will download the file): https://s3.immport.org/release/genelists/current/all_gene_lists.gmt?download=true
# clean_lines <- raw_lines[nzchar(trimws(raw_lines))]
# clean_lines <- clean_lines[grepl("\t", clean_lines)]
# tmp_clean <- tempfile(fileext = ".gmt")
# writeLines(clean_lines, tmp_clean)
# gsets <- getGmt(tmp_clean, sep = "\t")
# pathway_list <- geneIds(gsets)
# pathwayRDS <- saveRDS(pathway_list, file = "inst/extdata/immport_gene_pathway_map.rds")

###############
### Helpers ###
###############

#general function for pulling this mapping into memory, caching to avoid repeated file reads.
.LoadImmPortGeneMapping <- function() {
  if (!is.null(pkg.env$immport_gene_to_pathway)) {
    return(pkg.env$immport_gene_to_pathway)
  }

  rds_path <- system.file("extdata", "immport_gene_pathway_map.rds", package = "CellMembrane")
  if (!nzchar(rds_path)) {
    stop("immport_gene_pathway_map.rds not found. Run inst/scripts/BuildImmPortMapping.R to generate it.")
  }

  gene_to_pathway <- readRDS(rds_path)
  pkg.env$immport_gene_to_pathway <- gene_to_pathway

  return(gene_to_pathway)
}

  #helper for pathway counting in the permutation test. Takes a vector of gene symbols and returns a table of pathway counts.
  .CountPathways <- function(genes) {
    gene_to_pathway <- .LoadImmPortGeneMapping()
    all_pathway_hits <- unlist(gene_to_pathway[genes], use.names = FALSE)
    if (length(all_pathway_hits) == 0) return(integer(0))
    return(table(all_pathway_hits))
  }

#' @title Get ImmPort Gene Mapping
#' @description Returns the full named list mapping each gene symbol to its ImmPort pathways.
#' @return A named list where each name is an uppercase gene symbol and each value is a character vector of pathway names.
#' @export
GetImmPortGeneMapping <- function() {
  return(.LoadImmPortGeneMapping())
}


#' @title Get ImmPort Genes
#' @description Returns all unique gene symbols present in the ImmPort gene list.
#' @return A character vector of uppercase gene symbols.
#' @export
GetImmPortGenes <- function() {
  gene_to_pathway <- .LoadImmPortGeneMapping()
  return(names(gene_to_pathway))
}


#' @title Gene To ImmPort Pathway
#' @description Looks up the ImmPort pathway(s) associated with a given gene symbol.
#' @param geneSymbol A single gene symbol string (case-insensitive).
#' @return A character vector of pathway names, or NA if the gene is not found.
#' @export
GeneToImmPortPathway <- function(geneSymbol) {
  if (length(geneSymbol) != 1 || !is.character(geneSymbol)) {
    stop("geneSymbol must be a single character string.")
  }

  normalized_symbol <- toupper(trimws(geneSymbol))
  gene_to_pathway <- .LoadImmPortGeneMapping()
  matched_pathways <- gene_to_pathway[[normalized_symbol]]

  if (is.null(matched_pathways)) {
    return(NA_character_)
  }

  return(matched_pathways)
}

###########################
### Seurat Manipulation ###
###########################

#' @title Subset Seurat To ImmPort Genes
#' @description Subsets the features of a Seurat assay to only those present in the ImmPort gene list.
#' @param seuratObj A Seurat object.
#' @param assayName The assay to subset. Defaults to the active default assay.
#' @return The Seurat object with features reduced to the ImmPort gene intersection.
#' @export
SubsetSeuratToImmPortGenes <- function(seuratObj, assayName = NULL) {
  if (is.null(assayName)) {
    assayName <- Seurat::DefaultAssay(seuratObj)
  }

  if (!assayName %in% names(seuratObj@assays)) {
    stop("Assay not found in seuratObj: ", assayName)
  }

  all_features <- rownames(seuratObj[[assayName]])
  immport_genes <- GetImmPortGenes()

  # uppercase comparison to be consistent with how the mapping is keyed
  features_upper <- toupper(all_features)
  immport_set <- immport_genes

  matched_indices <- which(features_upper %in% immport_set)
  retained_features <- all_features[matched_indices]

  n_original <- length(all_features)
  n_retained <- length(retained_features)
  
  #handle feature intersection cases
  #null case
  if (n_retained == n_original) {
    print(paste0("All ", n_original, " features in assay '", assayName, "' are present in ImmPort. No subsetting needed."))
    return(seuratObj)
  } else if (n_retained == 0) {
    stop("No features in assay '", assayName, "' matched the ImmPort gene list. Check that gene symbols are HGNC-standard.")
  }
  #otherwise, report how many features are retained and subset the assay
  print(paste0("Subsetting assay '", assayName, "': retaining ", n_retained, " of ", n_original, " features present in ImmPort."))

  seuratObj <- subset(seuratObj, features = retained_features)

  return(seuratObj)
}

########################
### Statistical Test ###
########################

#' @title Permutation Test For ImmPort Pathway Occurrences
#' @description Tests whether ImmPort pathways are over-represented in a gene set relative to random draws
#' from the ImmPort gene universe. For each pathway observed in the input set, the observed gene count is
#' compared against a null distribution built by resampling the same number of genes from the full ImmPort
#' gene space.
#' @param geneSet A character vector of gene symbols (e.g., DE result genes). Case-insensitive.
#' @param nPermutations Number of random resamples used to build the null distribution. Default is 1000.
#' @param pAdjustMethod Multiple testing correction method passed to \code{p.adjust}. Default is \code{"BH"}.
#' @param seed Integer random seed for reproducibility. Default is \code{GetSeed()}.
#' @return A data.frame ordered by ascending PValue with columns:
#'   \itemize{
#'     \item \code{Pathway}: ImmPort pathway name.
#'     \item \code{ObservedCount}: Number of input genes mapping to the pathway.
#'     \item \code{MeanNullCount}: Mean count across all permutations.
#'     \item \code{PValue}: Empirical one-tailed p-value (proportion of permutations where null count >= observed).
#'     \item \code{PValueAdj}: Multiple-testing-adjusted p-value.
#'   }
#' @export
PermutationTestForImmPortPathwayOccurrences <- function(geneSet, nPermutations = 1000, pAdjustMethod = "BH", seed = GetSeed()) {
  if (!is.character(geneSet) || length(geneSet) == 0) {
    stop("geneSet must be a non-empty character vector.")
  }

  SetSeed(seed)
  
  gene_to_pathway <- .LoadImmPortGeneMapping()
  immport_universe <- names(gene_to_pathway)

  # normalize and intersect input gene set with ImmPort universe
  query_genes <- intersect(toupper(trimws(geneSet)), immport_universe)
  n_query <- length(query_genes)

  if (n_query == 0) {
    stop("None of the genes in geneSet were found in the ImmPort gene universe. Check that gene symbols are HGNC-standard.")
  }
  print(paste0("Retained ", n_query, " of ", length(geneSet), " input genes found in the ImmPort universe (", length(immport_universe), " total ImmPort genes)."))

  # observed pathway counts for the input gene set
  observed_table <- .CountPathways(query_genes)
  observed_pathways <- names(observed_table)

  if (length(observed_pathways) == 0) {
    stop("No pathway annotations were found for the intersected gene set.")
  }

  # build a binary membership matrix (rows = immport_universe, cols = observed_pathways)
  # so we can avoid per-gene loops for speed.
  all_genes_flat    <- rep(immport_universe, lengths(gene_to_pathway[immport_universe]))
  all_pathways_flat <- unlist(gene_to_pathway[immport_universe], use.names = FALSE)
  keep              <- all_pathways_flat %in% observed_pathways
  membership_matrix <- matrix(FALSE,
                              nrow = length(immport_universe),
                              ncol = length(observed_pathways))
  membership_matrix[cbind(match(all_genes_flat[keep], immport_universe),
                          match(all_pathways_flat[keep], observed_pathways))] <- TRUE

  # build null distribution: each permutation is a matrix row-subset + colSums
  null_counts <- matrix(0L, nrow = length(observed_pathways), ncol = nPermutations,
                        dimnames = list(observed_pathways, NULL))

  for (i in seq_len(nPermutations)) {
    sampled_idx      <- sample(length(immport_universe), size = n_query, replace = FALSE)
    null_counts[, i] <- colSums(membership_matrix[sampled_idx, , drop = FALSE])
  }

  # empirical one-tailed p-value: proportion of permutations where null >= observed.
  # R column-major recycling broadcasts observed_vec correctly across all permutation columns.
  observed_vec <- as.integer(observed_table[observed_pathways])
  p_values <- rowMeans(null_counts >= observed_vec)

  result_df <- data.frame(
    Pathway      = observed_pathways,
    ObservedCount = observed_vec,
    MeanNullCount = round(rowMeans(null_counts), 3),
    PValue       = p_values,
    PValueAdj    = p.adjust(p_values, method = pAdjustMethod),
    stringsAsFactors = FALSE
  )

  result_df <- result_df[order(result_df$PValue), ]
  rownames(result_df) <- NULL

  return(result_df)
}





