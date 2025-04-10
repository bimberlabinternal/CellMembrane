% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spatial.R
\name{DetectCellStructuresBasedOnCellType}
\alias{DetectCellStructuresBasedOnCellType}
\title{DetectCellStructuresBasedOnCellType}
\usage{
DetectCellStructuresBasedOnCellType(
  seuratObjectMetadata,
  cellTypeField = "cell_type",
  minimumClusterSizeCoefficient = 0.05,
  fovField = "fov",
  fovWhitelist = NULL,
  cellTypeWhiteList = NULL,
  xCoordinateField = "x_FOV_px",
  yCoordinateField = "y_FOV_px",
  substructureMetaDataFieldName = "Substructure",
  summarizeLocalResults = TRUE
)
}
\arguments{
\item{seuratObjectMetadata}{Metadata dataframe storing, at minimum, a cell type or "transcriptomic clustering" type field, a field denoting FOV, and two spatial coordinates.}

\item{cellTypeField}{A character field storing cell type annotations, but could be any discrete cell clustering assignment.}

\item{minimumClusterSizeCoefficient}{The percentage of the input data (after cell type splitting) to be considered as a feasible minimum cluster size. Determines how many "noise data points" will be added by DBSCAN.}

\item{fovField}{the metadata column that stores the Field of View information.}

\item{fovWhitelist}{An optional whitelist of FOVs. By default, the function will loop over all FOVs, which could be time consuming.}

\item{cellTypeWhiteList}{A vector of genes that constitute your substructure (e.g. c("Bcell", "BCell") for B cell follicles).}

\item{xCoordinateField}{The metadata column that stores the x coordinate information within the Field of View.}

\item{yCoordinateField}{The metadata column that stores the y coordinate information within the Field of View.}

\item{substructureMetaDataFieldName}{An annotation that will be concatenated during the results. "Local" FOV information will be concatenated using "fov + substructureMetaDataFieldName + a substructure index" within the columns of the metadata.}

\item{summarizeLocalResults}{An optional boolean that will wrap up the various substructureMetaDataFieldName columns into two single columns. One, which determines if a cell is within ANY of the defined substructures, stored in the output column "Within_Local + substructureMetaDataFieldName". The second is a metadata column that displays which of the local substructures the cell belongs in, concatenated as "Local + substructureMetaDataFieldName". "Local + substructureMetaDataFieldName + 0" is always the noise designation.}
}
\value{
Returns a dataframe containing columns related to the substructures found within the images at varying scopes. With summarizeLocalResults = FALSE, (number of FOVs) x (number of subtructures + 1) columns will be added. summarizeLocalResults rolls these high resolution results into two additional columns relative to the fovField.
}
\description{
A DBSCAN-based spatial cluster detection method to find dense cellular structures within spatial data. This function is currently implemented with FOV-based, iterative, classification in mind, but is technically extendable to global analysis of a single frame.
}
\examples{
\dontrun{
#Perform Cell Structure detection for B cell follicles. 

metadata <- DetectCellStructuresBasedOnCellType(seuratObjectMetadata, 
cellTypeField = "cell_type", 
minimumClusterSizeCoefficient = 0.05,
fovField = "fov",
fovWhitelist = 1,
cellTypeWhiteList = c("Bcell", "B_cell", "B.cell"),
xCoordinateField = "x_FOV_px", 
yCoordinateField = "y_FOV_px", 
substructureMetaDataFieldName = "BCF",
summarizeLocalResults = TRUE
)

#load/install packages for plotting
library(pacman)
p_load(ggplot2, dplyr, egg, patchwork)

#define plotting layout
layout <- "
#AAAA#
BBBCCC
"

#Plot results
ggplot(metadata \%>\% filter(fov == 1), 
aes(x = x_FOV_px, y = y_FOV_px, color = simple_cellType)) + 
geom_point() + 
egg::theme_article() + 
ggtitle('Cell Type assignment')
ggplot(metadata \%>\% filter(fov == 1), aes(x = x_FOV_px, y = y_FOV_px, color = factor(Local_BCF))) + 
geom_point() + 
egg::theme_article() +  
ggtitle('Specific substructure cell assignment') + 
ggplot(metadata \%>\% filter(fov == 1), 
aes(x = x_FOV_px, y = y_FOV_px, color = Within_Local_BCF)) + 
geom_point() + 
egg::theme_article() + 
ggtitle('Non-specific substructure cell assignment') + 
plot_layout(design = layout, guides = "collect")
}
}
