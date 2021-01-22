#' Table mapping gene symbol to CD gene name
#'
#'
#' @format A data frame with 394 rows and 4 variables:
#' \describe{
#'   \item{GeneSymbol}{Official GeneSymbol}
#'   \item{ApprovedName}{A description}
#'   \item{PreviousSymbols}{Previous symbols.  This is the field typically used for aliasing gene symbols}
#'   \item{Synonyms}{Additional aliases}
#'   ...
#' }
#' @source \url{https://www.genenames.org/data/genegroup/#!/group/471}
"cdGenes"