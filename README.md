![R Build and Checks](https://github.com/BimberLabInternal/CellMembrane/workflows/R%20Build%20and%20Checks/badge.svg)

# CellMembrane
An R package with wrappers and pipelines for single cell RNA-seq analysis

## Table of Contents
* [Overview](#overview)
* [Example Usage](#example)
* [Installation](#installation)

### <a name = "overview">Overview</a>



### <a name="example">Example Usage</a>


### <a name="installation">Installation</a>

```{r}
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlabinternal/cellmembrane', ref = 'master', dependencies = TRUE, upgrade = 'always')
```
    
Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/BimberLabInternal/packages/container/package/CellMembrane). We recommend using a specific release, which you can do using tags: 

```
docker pull ghcr.io/bimberlabinternal/cellmembrane:latest
```
