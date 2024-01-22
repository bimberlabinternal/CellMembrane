[![R Build and Checks](https://github.com/bimberlabinternal/CellMembrane/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bimberlabinternal/CellMembrane/actions/workflows/R-CMD-check.yaml)

# CellMembrane
An R package with wrappers and pipelines for single cell RNA-seq analysis

## Table of Contents
* [Overview](#overview)
* [Example Usage](#example)
* [Installation](#installation)

### <a name = "overview">Overview</a>

This package was created to run the Bimber Lab's single-cell RNA-seq pipelines. It primarily provides high-level wrappers around tools (hence the name). The general idea is to facilitate pipelining tools, make those pipelines fault-tolerant across input data, and to provide QC/visualizations when running those pipelines to help interpret results or issues.

### <a name="installation">Installation</a>

```{r}
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlabinternal/cellmembrane', dependencies = TRUE, upgrade = 'always')
```
    
Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/BimberLabInternal/packages/container/package/CellMembrane). We recommend using a specific release, which you can do using tags: 

```
docker pull ghcr.io/bimberlabinternal/cellmembrane:latest
```
