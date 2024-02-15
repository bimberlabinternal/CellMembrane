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
#' @param gseaPackage This parameter specifies the R package to utilize for conducting GSEA analysis.Possible options are("clusterProfiler", "fgsea")
#' @param msigdbSubcategory MSigDB sub-collection abbreviation, such as CGP or BP
#' @param genesToExclude A vector of genes to exclude from the pathway gene sets.
#' @return A Seurat object with the result of the GSEA analysis stored in the 'Misc' slot under 'NES_' followed by gseaPackage.
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
                              gseaPackage = "clusterProfiler",
                              msigdbSubcategory = NULL,
                              genesToExclude = NULL) {
  
  if (gseaPackage %in% c("fgsea", "clusterProfiler")) {
    print(paste("Using", gseaPackage, "R package to run GSEA analysis"))
  
    } else{
    stop(
      "Invalid value for gseaPackage parameter. Please choose from: ",
      paste(c("fgsea", "clusterProfiler"), collapse = ", ")
    )
  }
  gene_ranks <- presto::wilcoxauc(seuratObj, group_by = groupField)
  
  msigdb_df <- msigdbr::msigdbr(species = msigdbSpecies,
                                category = msigdbCategory,
                                subcategory = msigdbSubcategory) %>%
    dplyr::select(gs_name, entrez_gene, gene_symbol) %>%
    dplyr::distinct()
  
  if (!is.null(selectedPathways)) {
    msigdb_df <-
      msigdb_df %>% dplyr::filter(gs_name %in% selectedPathways)
  }
  
  if (!is.null(genesToExclude)) {
    msigdb_df <- msigdb_df %>% dplyr::filter(!gene_symbol %in% genesToExclude)
  }
  
  gsea_result <- data.frame()
  
  if (gseaPackage == "fgsea") {
    fgsea_sets <- msigdb_df %>% split(x = .$gene_symbol, f = .$gs_name)
  } else if (gseaPackage == "clusterProfiler") {
    msigdb_df <- msigdb_df %>%
      dplyr::select(gs_name, entrez_gene) %>%
      dplyr::distinct()
  }
  
  for (groupName in unique(seuratObj@meta.data[, groupField])) {
    group_genes <- gene_ranks %>%
      dplyr::filter(group == groupName) %>%
      dplyr::arrange(desc(auc)) %>%
      dplyr::select(feature, auc)
    
    ranks <- tibble::deframe(group_genes)
    
    if (gseaPackage == "fgsea") {
      temp_result <-
        fgsea::fgsea(fgsea_sets, stats = ranks, scoreType = scoreType) %>%
        mutate(pvalue = pval,
               setSize = size) %>%
        dplyr::select(-pval, -size, -leadingEdge)
      
    } else if (gseaPackage == "clusterProfiler") {
      ortholog_genes <-
        babelgene::orthologs(names(ranks), species = msigdbSpecies, human = FALSE)
      group_genes <-
        group_genes %>% 
        left_join(ortholog_genes, by = c("feature" = "symbol")) %>% 
        filter(!is.na(entrez))
      
      geneList <- group_genes[, "auc"]
      names(geneList) <- as.character(group_genes[, "entrez"])
      geneList <- sort(geneList, decreasing = TRUE)
      
      temp_result <-
        clusterProfiler::GSEA(
          geneList,
          TERM2GENE = msigdb_df,
          scoreType = 'pos',
          pvalueCutoff = 1
        )
      temp_result <- temp_result@result %>%
        dplyr::mutate(pathway = ID,
                      ES = enrichmentScore,
                      padj = p.adjust) %>%
        dplyr::select(-ID, -Description,-enrichmentScore,-p.adjust)
    }
    
    temp_result <- temp_result %>%
      tibble::as_tibble() %>%
      dplyr::mutate(groupName = groupName)
    
    gsea_result <-
      rbind(gsea_result, temp_result) %>%
      dplyr::group_by(pathway) %>%
      dplyr::arrange(desc(NES), .by_group = TRUE)
  }
  
  # Plot NES for each pathway - Lollipop plots
  
  for (path in unique(gsea_result$pathway)) {
    nes_data <-
      gsea_result %>% 
      filter(pathway == path) %>% 
      arrange(desc(NES))
    
    p <- ggplot(nes_data, aes(reorder(groupName, NES), NES)) +
      geom_segment(aes(
        y = -Inf,
        yend = NES,
        x = reorder(groupName, NES),
        xend = reorder(groupName, NES)
      )) +
      geom_point(aes(
        y = NES,
        x = reorder(groupName, NES),
        fill = -log10(pvalue)
      ),
      shape = 21,
      size = 5) +
      coord_flip() +
      labs(x = "Group", y = "Normalized Enrichment Score",
           title = path) +
      scale_fill_gradientn(colors = c("white", "red")) +
      egg::theme_presentation(base_size = 12) +
      theme(axis.ticks.y = element_blank())
    
    print(p) 
  }
  
  Misc(seuratObj, slot = paste("NES", gseaPackage, sep = "_")) <- gsea_result
  return(seuratObj)
  
}