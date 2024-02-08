#' @import ggplot2


#' @title PathwayEnrichment
#'
#' @description Runs GSEA analysis on a seurat object
#' @param seuratObj The Seurat object containing the data. 
#' @param groupField Field on which to group to calculate AUC
#' @param msigdbSpecies MSigDBSpecies name to retrieve gene sets from. Run msigdbr_species() to see available species.
#' @param msigdbCategory MSigDB collection abbreviation, such as H or C1. Run msigdbr_collections() to see available collections.
#' @param scoreType This parameter defines the GSEA score type.Possible options are("std","pos", "neg")
#' @param selectedPathways A set of pathways to test.
#' @param msigdbSubcategory MSigDB sub-collection abbreviation, such as CGP or BP
#' @return A seurat object with NES values for each pathway appended to each cell 
#' @export

PathwayEnrichment <- function(seuratObj,
                              groupField,
                              msigdbSpecies = "Macaca mulatta",
                              msigdbCategory = "H",
                              scoreType = "pos",
                              selectedPathways = c(
                                "HALLMARK_TGF_BETA_SIGNALING",
                                "HALLMARK_GLYCOLYSIS",
                                "HALLMARK_IL2_STAT5_SIGNALING",
                                "HALLMARK_INFLAMMATORY_RESPONSE",
                                "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                                "HALLMARK_G2M_CHECKPOINT",
                                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                "HALLMARK_DNA_REPAIR",
                                "HALLMARK_APOPTOSIS",
                                "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                                "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                                "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
                              ),
                              msigdbSubcategory = NULL) {
  gene_ranks <- presto::wilcoxauc(seuratObj, group_by = groupField)
  
  msigdb_df <- msigdbr::msigdbr(species = msigdbSpecies,
                                category = msigdbCategory,
                                subcategory = msigdbSubcategory) %>%
    dplyr::select(gs_name, entrez_gene, gene_symbol) %>%
    dplyr::distinct()
  
  if (!is.null(selectedPathways)) {
    msigdb_df <- msigdb_df %>% dplyr::filter(gs_name %in% selectedPathways)
  }
  
  fgsea_sets <- msigdb_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  gsea_result <- data.frame()
  
  for (path in unique(msigdb_df$gs_name)) {
    for (groupName in unique(seuratObj@meta.data[, groupField])) {
      group_genes <- gene_ranks %>%
        dplyr::filter(group == groupName) %>%
        dplyr::arrange(desc(auc)) %>%
        dplyr::select(feature, auc)
      
      ranks <- tibble::deframe(group_genes)
      
      fgsea_result <-
        fgsea::fgsea(fgsea_sets, stats = ranks, scoreType = scoreType)
      
      fgsea_result <- fgsea_result %>%
        tibble::as_tibble() %>%
        dplyr::filter(pathway == path) %>%
        dplyr::mutate(groupName = groupName)
      
      gsea_result <-
        rbind(gsea_result, fgsea_result) %>%  dplyr::arrange(desc(NES))
    }
  }
  
  return(gsea_result)
}