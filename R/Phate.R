# Based on Seurat v2 code: https://www.rdocumentation.org/packages/Seurat/versions/2.3.4/topics/RunPHATE

.RunPHATE <- function(
  object,
  assay = NULL,
  n.components = 2L,
  knn = 5L,
  decay = 40L,
  n.landmark=2000L,
  gamma = 1,
  t = "auto",
  mds.solver = "sgd",
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max=100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = NA,
  reduction.key = "PHATE_",
  k = NULL,
  alpha = NULL,
  ...
) {
  if (!is.null(k)) {
    message("Argument k is deprecated. Using knn instead.")
    knn <- k
  }
  if (!is.null(alpha)) {
    message("Argument alpha is deprecated. Using decay instead.")
    decay <- alpha
  }

  phate_output <- phateR::phate(
    object,
    ndim = n.components,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    knn.dist.method = knn.dist.method,
    init = NULL,
    mds.solver = mds.solver,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed.use,
    ...
  )
  phate_output <- as.matrix(phate_output)
  colnames(x = phate_output) <- paste0(reduction.key, 1:ncol(x = phate_output))
  rownames(x = phate_output) <- rownames(object)
  phate.reduction <- CreateDimReducObject(
    embeddings = phate_output,
    key = reduction.key,
    assay = assay
  )
  return(phate.reduction)
}

#' Run PHATE
#'
#' PHATE is a data reduction method specifically designed for visualizing
#' **high** dimensional data in **low** dimensional spaces.
#' To run, you must first install the `phate` python
#' package (e.g. via pip install phate). Details on this package can be
#' found here: \url{https://github.com/KrishnaswamyLab/PHATE}. For a more in depth
#' discussion of the mathematics underlying PHATE, see the bioRxiv paper here:
#' \url{https://www.biorxiv.org/content/early/2017/12/01/120378}.
#'
#' @param dims Which dimensions to use as input features, used only if
#' \code{features} is NULL
#' @param source This can either be PCA or counts (meaning it will be run on the raw counts)
#' @param features If set, run PHATE on this subset of features (instead of running on a
#' set of reduced dimensions). Not set (NULL) by default; \code{dims} must be NULL to run
#' on features
#' @param assay Assay to pull data for when using \code{features}
#' @param n.components Total number of dimensions to embed in PHATE.
#' @param knn int, optional, default: 5
#' number of nearest neighbors on which to build kernel
#' @param decay int, optional, default: 40
#' sets decay rate of kernel tails.
#' If NA, alpha decaying kernel is not used
#' @param n.landmark int, optional, default: 2000
#' number of landmarks to use in fast PHATE
#' @param gamma float, optional, default: 1
#' Informational distance constant between -1 and 1.
#' `gamma=1` gives the PHATE log potential, `gamma=0` gives
#' a square root potential.
#' @param t int, optional, default: 'auto'
#' power to which the diffusion operator is powered
#' sets the level of diffusion
#' @param mds.solver {'sgd', 'smacof'}, optional, default: 'sgd'
#' which solver to use for metric MDS. SGD is substantially faster,
#' but produces slightly less optimal results. Note that SMACOF was used
#' for all figures in the PHATE paper.
#' @param knn.dist.method string, optional, default: 'euclidean'.
#' recommended values: 'euclidean', 'cosine', 'precomputed'
#' Any metric from `scipy.spatial.distance` can be used
#' distance metric for building kNN graph. If 'precomputed',
#' `data` should be an n_samples x n_samples distance or
#' affinity matrix. Distance matrices are assumed to have zeros
#' down the diagonal, while affinity matrices are assumed to have
#' non-zero values down the diagonal. This is detected automatically using
#' `data[0,0]`. You can override this detection with
#' `knn.dist.method='precomputed_distance'` or
#' `knn.dist.method='precomputed_affinity'`.
#' @param mds.method string, optional, default: 'metric'
#' choose from 'classic', 'metric', and 'nonmetric'
#' which MDS algorithm is used for dimensionality reduction
#' @param mds.dist.method string, optional, default: 'euclidean'
#' recommended values: 'euclidean' and 'cosine'
#' @param t.max int, optional, default: 100.
#' Maximum value of t to test for automatic t selection.
#' @param npca int, optional, default: 100
#' Number of principal components to use for calculating
#' neighborhoods. For extremely large datasets, using
#' n_pca < 20 allows neighborhoods to be calculated in
#' log(n_samples) time.
#' @param plot.optimal.t boolean, optional, default: FALSE
#' If TRUE, produce a plot showing the Von Neumann Entropy
#' curve for automatic t selection.
#' @param verbose `int` or `boolean`, optional (default : 1)
#' If `TRUE` or `> 0`, print verbose updates.
#' @param n.jobs `int`, optional (default: 1)
#' The number of jobs to use for the computation.
#' If -1 all CPUs are used. If 1 is given, no parallel computing code is
#' used at all, which is useful for debugging.
#' For n_jobs below -1, (n.cpus + 1 + n.jobs) are used. Thus for
#' n_jobs = -2, all CPUs but one are used
#' @param seed.use int or `NA`, random state (default: `NA`)
#' @param reduction.name dimensional reduction name, specifies the position in
#' the object$dr list. phate by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. PHATE_ by default
#' @param k Deprecated. Use `knn`.
#' @param alpha Deprecated. Use `decay`.
#'
#' @return Returns a Seurat object containing a PHATE representation
#'
#' @references Moon K, van Dijk D, Wang Z, Gigante S,
#' Burkhardt D, Chen W, van den Elzen A,
#' Hirn M, Coifman R, Ivanova N, Wolf G and Krishnaswamy S (2017).
#' "Visualizing Transitions and Structure for High Dimensional Data
#' Exploration." _bioRxiv_, pp. 120378. doi: 10.1101/120378
#' (URL: http://doi.org/10.1101/120378),
#' <URL: https://www.biorxiv.org/content/early/2017/12/01/120378>.
#'
#' @export
#'
RunPHATE <- function(
  object,
  dims = NULL,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  n.components = 2L,
  knn = 5L,
  decay = 40L,
  n.landmark=2000L,
  gamma = 1,
  t = "auto",
  mds.solver = 'sgd',
  knn.dist.method = "euclidean",
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed.use = NA,
  reduction.name = "phate",
  reduction.key = "PHATE_",
  k = NULL,
  alpha = NULL,
  ...
) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('GMM_Demux')) {
    stop('GMM_Demux has not been installed!')
  }

  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = dims) && !is.null(x = features)) {
    stop("Please specify only one of the following arguments: dims or features")
  }

  if (!is.null(x = features)) {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- assay %||% DefaultAssay(object = object[[reduction]])
  } else {
    data.use <- Embeddings(object[[reduction]])
    assay <- assay %||% DefaultAssay(object = object[[reduction]])
  }
  object[[reduction.name]] <- .RunPHATE(
    object = data.use,
    assay = assay,
    n.components = n.components,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    knn.dist.method = knn.dist.method,
    mds.solver = mds.solver,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed.use = seed.use,
    reduction.key = reduction.key,
    k = k,
    alpha = alpha,
    ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}