% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseudoBulk.R
\name{FilterPseudobulkContrasts}
\alias{FilterPseudobulkContrasts}
\title{FilterPseudobulkContrasts}
\usage{
FilterPseudobulkContrasts(
  logicList = NULL,
  design = NULL,
  useRequireIdenticalLogic = T,
  requireIdenticalFields = NULL,
  filteredContrastsOutputFile = "./filtered_contrasts.tsv"
)
}
\arguments{
\item{logicList}{The list that defines the study design. Please see examples for the setup and use. This list be a list of lists where each sub-list has three entries. The first defines the metadata field (e.g. Tissue or Timepoint) to which the logic gate will be applied. The second defines the logic gate that will be used for that field (one of: all (AND), any (OR), xor, nand, nor, or xnor. The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.}

\item{design}{a design/model matrix returned by DesignModelMatrix().}

\item{useRequireIdenticalLogic}{Whether or not to apply require_identical logic using requireIdenticalFields. It is possible to use logic gates to keep metadata fields constant, but using this saves time & effort/iteration.}

\item{requireIdenticalFields}{The metadata columns of the SeuratObj that you wish to keep constant. For example: defining "Tissue" in this vector would filter all contrasts comparing Liver to PBMC, as their tissue is not identical.}

\item{filteredContrastsOutputFile}{The file to write the list of filtered contrasts to (one pair of samples per row, separated into their metadata fields).}
}
\value{
A dataframe of pairwise contrasts, with column names indicating the directionality (positive vs negative) and associated metadata column (inferred by contrast_columns).
}
\description{
This is designed to accept a study design defined by a series of logical gates applied to metadata fields, test them for equivalence to some criterion, and then filter pairwise contrasts accordingly.
}
\details{
\tabular{cccccccccc}{
   \strong{Filtered} \tab \emph{Reason} \tab \strong{positive_contrast_Challenge} \tab \emph{positive_contrast_Tissue} \tab \strong{positive_contrast_cell_type} \tab \emph{positive_contrast_SampleType} \tab \strong{negative_contrast_Challenge} \tab \emph{negative_contrast_Tissue} \tab \strong{negative_contrast_cell_type} \tab \emph{negative_contrast_SampleType} \cr
   yes \tab non-identical cell_types \tab Mtb \tab Spleen \tab \emph{\strong{Myeloid}} \tab Necropsy \tab Mock.challenged \tab Spleen \tab \emph{\strong{Bcell}} \tab Necropsy \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
   no \tab N/A \tab Mtb \tab Lung \tab T.NK \tab Necropsy \tab Mock.challenged \tab Lung \tab T.NK \tab Necropsy \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
   yes \tab fails Challenge xor gate \tab \emph{\strong{Mock.challenged}} \tab Lung \tab T.NK \tab Necropsy \tab \emph{\strong{Mock.challenged}} \tab Lung \tab T.NK \tab Necropsy \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
    \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \tab  \cr
   yes \tab fails SampleType AND gate \tab Mock.challenged \tab MesLN \tab Bcell \tab \emph{\strong{Baseline}} \tab Mock.challenged \tab MesLN \tab Bcell \tab \emph{\strong{Necropsy}} \cr
}
}
\examples{
\dontrun{
#Set up design matrix
design <- DesignModelMatrix(pseudobulked_seuratObj, 
                            contrast_columns = c('Challenge', 'Tissue', 'cell_type', 'SampleType'), 
                            sampleIdCol = 'cDNA_ID')


#Form the study design. This should be a list of gates (each specified as a list itself) with three elements. 
#The first element defines the metadata column or field that will be used for filtering. 
#The second defines the logic gate (AND -> any, OR -> any, xor, xnor, nand, nor). 
#The third defines the specific value of the metadata field that will be tested for equivalence against the positive and negative contrasts within the gate.
logicList <- list(
                    list('Challenge', 'xor', 'Mock.challenged'), 
                    list('SampleType', 'all', 'Necropsy')
            )


#Finally, enumerate all possible contrasts, then filter and return them according to the study design defined by logicList.
#If you want to perform "Like-vs-Like" contrasts (i.e. T cell vs T cell within a cell_type metadata field), specify those columns as required to be identical and set useRequireIdenticalLogic to TRUE.
filtered_contrasts <- FilterPseudobulkContrasts(logicList = logicList, 
                                                design = design, 
                                                useRequireIdenticalLogic = T, 
                                                requireIdenticalFields = c('Tissue','cell_type', 'SampleType'), 
                                                filteredContrastsOutputFile = './filtered_contrasts.tsv')
                                                
}
}
