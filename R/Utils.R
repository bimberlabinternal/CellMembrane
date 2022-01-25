#' @importFrom graphics abline boxplot legend lines plot points segments
#' @importFrom methods new
#' @importFrom stats approxfun cmdscale dist kmeans lm na.omit prcomp quantile sd setNames wilcox.test
#' @importFrom utils head read.csv read.table tail write.csv write.table

utils::globalVariables(
  names = c('X', 'Y', 'SortOrder'),
  package = 'CellMembrane',
  add = TRUE
)

pkg.env <- new.env(parent=emptyenv());

pkg.env$RANDOM_SEED <- 1234
set.seed(pkg.env$RANDOM_SEED)

.GetCCGenes <- function(){
  # Cell cycle genes were obtained from the Seurat example (See regev_lab_cell_cycle_genes.txt)
  # and stored using use_data(internal = T) (https://github.com/r-lib/usethis and use_data)
  # cc.genes
  # g2m.genes.orig

  return(cc.genes)
}

.GetSPhaseGenes <- function(){
  return (.GetCCGenes()[1:43])
}

.GetG2MGenes <- function() {
  return(unique(c(g2m.genes.orig, .GetCCGenes()[44:97])))
}

.WriteLogMsg <- function(msg, prefixTime = TRUE, file = 'CellMembrane.log.txt') {
  if (prefixTime) {
    msg <- paste0(Sys.time(), ' ', msg)
  }

  write(msg, file = file, append = T)
}

.InferPerplexityFromSeuratObj <- function(seuratObj, perplexity = 30) {
  return(.InferPerplexity(ncol(seuratObj), perplexity))
}

.InferPerplexity <- function(sampleNumber, perplexity = 30) {
  if (sampleNumber - 1 < 3 * perplexity) {
    print(paste0('Perplexity is too large for the number of samples: ', sampleNumber))
    perplexityNew <- floor((sampleNumber - 1) / 3)
    print(paste0('lowering from ', perplexity, ' to: ', perplexityNew))
    perplexity <- perplexityNew
  }

  return(perplexity)
}

.PossiblyAddBarcodePrefix <- function(seuratObj, datasetId, datasetName = NULL) {
  if (!('BarcodePrefix' %in% names(seuratObj@meta.data))) {
    print(paste0('Adding barcode prefix: ', datasetId))
    seuratObj <- RenameCells(object = seuratObj, add.cell.id = datasetId)
    seuratObj[['BarcodePrefix']] <- c(datasetId)
    seuratObj[['DatasetId']] <- c(datasetId)
    if (!is.null(datasetName)) {
      seuratObj[['DatasetName']] <- datasetName
    }
  } else {
    print('Barcode prefix already added')
  }

  return(seuratObj)
}


#' @title Update Gene Model
#'
#' @description Substitutes LOC genes to more common gene IDs
#' @param features a vector of features to be updated to more common gene IDS
#' @param predictions Logical describing whether or not to include lower quality/speculative updates to the gene model
#' @param automatic Logical describing whether or not to attempt to automatically replace LOC genes
#' @param ambiguity_suffix String to be concatenated at the end of genes with non-conclusive 1:many LOC -> symbol mappings
#' @importFrom biomaRt useMart getLDS
#' @importFrom Biobase lcPrefix
#' @importFrom dplyr full_join
#' @import org.Mmu.eg.db 

.UpdateGeneModel <- function(features, predictions = T, automatic = T, ambiguity_suffix = ".PERHAPS"){
# TRAC = LOC710951 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJ95TNA1016&from=84561541&to=84565299)
features <- gsub("LOC710951", "TRAC", features)

# TRBC1 = LOC114677140 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJ9S5K0Z01R&from=169372796&to=169374268)
# TRBC2 maps to LOC114677140 as well
features <- gsub("LOC114677140", "TRBC1", features)
#features <- gsub("LOC114677140", ,"TRBC2", features)

#TRDC = LOC711031 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJA2XTBW01R&from=84480217&to=84482757)
features <- gsub("LOC711031", "TRDC", features)

# TRGC1 = LOC720538 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJAHHBHJ016&from=76878214&to=76881635)
features <- gsub("LOC720538", "TRGC1", features)

#TRGC2 = LOC705095 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041756.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJASX36T016&from=76909393&to=76914320)
features <- gsub("LOC705095", "TRGC2", features)

#TRDV1 = LOC720456 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=2&RID=UJB64E6Z013&from=84055100&to=84055663)
#TRDV1 also maps to TRDC (LOC710951) so use that if the variable region seems odd to include
features <- gsub("LOC720456", "TRDV1", features)

#IGKC = LOC701504 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041766.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBE65YT013&from=18130540&to=18130862)
features <- gsub("LOC701504", "IGKC", features)

#IGHG3 = LOC114679691 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBJKUU9016&from=168001652&to=168005818) 
#These are gamma-2 like instead of gamma-3 like, but they map well. 
features <- gsub("LOC114679691", "IGHG3", features)

#IGHG1 = LOC708891 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJBVCMEA013&from=168070075&to=168071681)
features <- gsub("LOC708891", "IGHG1", features)

#IGHA1 = LOC720839 (source: https://www.ncbi.nlm.nih.gov/nucleotide/NC_041760.1?report=genbank&log$=nuclalign&blast_rank=1&RID=UJC2E4S8013&from=167910270&to=167911722)
features <- gsub("LOC720839", "IGHA1", features)

#HLA.DP1 = LOC114669810 (source: https://www.ncbi.nlm.nih.gov/nuccore/XM_028842593.1)
features <- gsub("LOC114669810", "HLA.DPA1", features)

#FCN1 = LOC712405 (source: https://www.ncbi.nlm.nih.gov/nuccore/XM_015116399.2)
features <- gsub("LOC712405", "FCN1", features)

#XCL1 = LOC100423131 (source: https://www.ncbi.nlm.nih.gov/gene/100423131)
features <- gsub("LOC100423131", "XCL1.2", features)

#CCL4 = LOC100430627 (source: https://www.ncbi.nlm.nih.gov/gene/100430627)
#LOC100426632 also maps to CCL4
features <- gsub("LOC100430627", "CCL4", features)

#CCL3L1 = LOC100426537 (source: https://www.ncbi.nlm.nih.gov/gene/100426537)
#alternative name: CCL3L3
features <- gsub("LOC100426537", "CCL3L1", features)

  #these are predicted proteins/isoforms that seem less annotated
  if(predictions){
    #TYMP = LOC716161 (source: https://www.ncbi.nlm.nih.gov/protein/1622837354)
    features <- gsub("LOC716161", "TYMP", features)
    
    #Automatically pull in gene names from ensembl and entrez based annotations
    if(automatic){
      #Set up Ensembl translations
      rhesus <- useMart('ensembl', dataset="mmulatta_gene_ensembl")
      human <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
      
      table_merged <- getLDS(
        mart = rhesus,
        attributes = c('ensembl_gene_id', 'vgnc', 'external_gene_name', 'chromosome_name'),
        martL = human,
        attributesL = c('ensembl_gene_id','hgnc_symbol','gene_biotype', 'chromosome_name'))
      
      #Set up Entrez translations
      symbols <- mapIds(org.Mmu.eg.db, keys = features, keytype = "SYMBOL", column = "ENSEMBL")
      entrez_mapping <- as.data.frame(symbols) 
      entrez_mapping$gene_name <- rownames(entrez_mapping)
      
      #Join the tables based on common ensembl ID
      entrez_and_ensembl <- full_join(table_merged, entrez_mapping, by = c("Gene.stable.ID" = "symbols"))
      
      #Gene.name is mapped by Ensembl and gene_name is mapped by Entrez
      translation_table <- entrez_and_ensembl[rownames(seuratObj) %in% entrez_and_ensembl$Gene.name | rownames(seuratObj) %in% entrez_and_ensembl$gene_name, c("HGNC.symbol", "Gene.name", "Chromosome.scaffold.name", "Gene.type",  "gene_name")]
      #remove genes that aren't mapped by symbol by either entrez or ensembl
      translation_table<- translation_table[!is.na(translation_table$gene_name) & !is.na(translation_table$Gene.name),]
      
      
      #Iterate through features and attempt to map LOCs
      features_before_automatic_replacement <- features
      for (gene in features){
        if (gene %in% translation_table$gene_name){
          #regex match LOC genes
          if(grepl("^LOC", gene)){
            #create (potentially a vector) of NCBI LOC -> HGNC + ENTREZ translations
            vector_of_translations <- unique(translation_table[translation_table$gene_name == gene, c("HGNC.symbol", "Gene.name")])
            vector_of_translations <- vector_of_translations[!is.na(vector_of_translations)]
            vector_of_translations <- unique(vector_of_translations[!(vector_of_translations == "")])
            if(!identical(vector_of_translations ,character(0))){
            #Check if LOC genes have 1 -> many mapping
            if (length(vector_of_translations) == 1){
              #If the mapping is one to one, simply accept
              features[features == gene] <- vector_of_translations
            } else if (length(vector_of_translations) > 1){
              #If there are still multiple gene names, then attempt to find the best common prefix/gene name
              if (length(vector_of_translations) >1){
                common_prefix <- Biobase::lcPrefix(vector_of_translations)
                #if there's no common substring, abort and revert to using the LOC gene name
                if (common_prefix == ""){
                  common_prefix <- gene
                  features[features == gene] <- common_prefix
                } else{
                  #suffix denoting confusion for manual analysis and to omit from automatic downstream analysis
                  common_prefix <- paste0(common_prefix, ambiguity_suffix)
                  #replace gene name with common_prefix
                  features[features == gene] <- common_prefix
                }
              }
            } 
          }
        }
      }
      }
      #Print summary of replaced features and ensure suffixed features are unique
      print(paste0(length(features_before_automatic_replacement[grepl("^LOC",features_before_automatic_replacement)]) - length(features[grepl("^LOC",features)]), " LOC genes replaced"))
      print(paste0(length(features[grepl(paste0(ambiguity_suffix), features)]), " low confidence gene-names replaced"))
      features <- make.unique(features)
  }
}
return(features)
}

#' @title Set random seed
#'
#' @description Sets the seed used for R‘s random number generator, which should be used in all internal functions
#' @param seed The random seed
#' @export
SetSeed <- function(seed) {
  pkg.env$RANDOM_SEED <- seed
  set.seed(pkg.env$RANDOM_SEED)
}

#' @title Get random seed
#'
#' @description Sets a random seed, which should be used in all internal functions
#' @param seed The random seed
#' @export
GetSeed <- function(seed) {
  return(pkg.env$RANDOM_SEED)
}

#' @title Rename a vector of genes using CD nomenclature
#'
#' @description This compares a vector of genes to CD nomenclature, based on https://www.genenames.org/data/genegroup/#!/group/471. Any gene symbol that matches will have the CD name appended to the end (i.e. KLRB1 becomes "KLRB1 (CD161)")
#' @param inputGenes A vector of genes to be aliased.
#' @export
RenameUsingCD <- function(inputGenes) {
  vals <- CellMembrane::cdGenes[c('GeneSymbol', 'PreviousSymbols', 'Synonyms')]
  vals$Aliases <- sapply(paste0(vals$PreviousSymbols, ',', vals$Synonyms), function(x){
    x <- unlist(strsplit(x, split = ','))
    x <- x[x != '']
    x <- x[grep(x, pattern = '^CD')]
    x <- sort(x)
    x <- unique(toupper(x))

    return(paste0(x, collapse = ','))
  })
  vals <- vals[!is.na(vals$Aliases) & vals$Aliases != '',]
  vals <- vals[c('GeneSymbol','Aliases')]
  vals

  toMerge <- data.frame(Name = inputGenes, SortOrder = 1:length(inputGenes))
  toMerge <- merge(toMerge, vals, by.x = 'Name', by.y = 'GeneSymbol', all.x = T)
  sel <- !is.na(toMerge$Aliases)
  toMerge$Name[sel] <- paste0(toMerge$Name[sel], ' (', toMerge$Aliases[sel], ')')
  toMerge <- arrange(toMerge, SortOrder)

  return(toMerge$Name)
}

#' @title Update Seurat Object Gene Model From Monkey LOCs to Human Gene IDs
#'
#' @description Uses UpdateGeneModel to change LOC genes to more common human gene IDs in the RNA assay
#' @param seuratObj The Seurat Object to be updated
#' @param verbose This prints the slot and assay number, primarily for debugging once seurat updates.
#' @export

UpdateMacaqueMmul10NcbiGeneSymbols <- function(seuratObj, verbose = T){
  if (class(seuratObj)[[1]] != "Seurat"){
    stop(paste0('Please provide a Seurat Object'))
  }
  for (assay in names(seuratObj@assays)){
    for (slot in names(attributes(seuratObj@assays[[assay]]))){
      if(verbose){
        print(paste("Updating Gene Names in Assay:", assay, "Slot:", slot))
      }
      #updating common gene expression slots
      # @data could be either S4 or a matrix depending on normalization of @counts
      if (any(grepl(slot, c("counts", "data", "scale.data")))){
        if (typeof(attr(x = seuratObj@assays[[assay]], which = slot)) == "S4"){
          attr(x = seuratObj@assays[[assay]], which = slot)@Dimnames[[1]] <- .UpdateGeneModel(attr(x = seuratObj@assays[[assay]], which = slot)@Dimnames[[1]])
        }
        if (typeof(attr(x = seuratObj@assays[[assay]], which = slot)) == "double" | typeof(attr(x = seuratObj@assays[[assay]], which = slot)) == "integer"){
          rownames(attr(x = seuratObj@assays[[assay]], which = slot)) <- .UpdateGeneModel(rownames(attr(x = seuratObj@assays[[assay]], which = slot)))
        }
        
      }
      #updating variable features (expected to be a vector)
      if (slot == "var.features"){
        attr(x = seuratObj@assays[[assay]], which = slot) <- .UpdateGeneModel(attr(x = seuratObj@assays[[assay]], which = slot))
      }
      #updating meta.features (expected to be a matrix)
      if (slot == "meta.features"){
        rownames(attr(x = seuratObj@assays[[assay]], which = slot)) <- .UpdateGeneModel(rownames(attr(x = seuratObj@assays[[assay]], which = slot)))
        #feature metadata is tracked by assay, so this matrix could be empty.
        #This is an update that seems to be CellMembrane specific
        if(dim(attr(x = seuratObj@assays[[assay]], which = slot))[[2]] != 0 & any(grepl("GeneId", names(attr(x = seuratObj@assays[[assay]], which = slot)))))
          attr(x = seuratObj@assays[[assay]], which = slot)$GeneId <- .UpdateGeneModel(attr(x = seuratObj@assays[[assay]], which = slot)$GeneId)
      }
    }
  }
  #this block attempts to update the features in PCA reductions. (expected to be two matrices)
  if (any(grepl("pca", names(seuratObj@reductions)))){
    rownames(seuratObj@reductions$pca@feature.loadings) <- .UpdateGeneModel(rownames(seuratObj@reductions$pca@feature.loadings))
    rownames(seuratObj@reductions$pca@feature.loadings.projected) <- .UpdateGeneModel(rownames(seuratObj@reductions$pca@feature.loadings.projected))
  }
  return(seuratObj)
}