#' @include Utils.R
#' @include Preprocessing.R
#' @import Seurat

utils::globalVariables(
  names = c('Topics'),
  package = 'CellMembrane',
  add = TRUE
)


#' @title DoLdaParameterScan
#'
#' @description This will run LDA on the target assay
#' @param seuratObj A Seurat object.
#' @param ntopics Passed to runLDA.
#' @param normalizationMethod The method used for Seurat::NormalizeData()
#' @param varFeatures The number of variable features to use in the LDA model. The more features that are used, the slower the model will run and the more noise that will be introduced, but the model will be more complete in representing your entire dataset.
#' @param randomSeed Passed to runLDA seed.number argument
#' @param assayName The name of the source assay
#' @param nCores The number of cores to use
#' @export
DoLdaParameterScan <- function(seuratObj, ntopics = seq(5, 50, by=5), normalizationMethod = "CLR", varFeatures = 5000, randomSeed = GetSeed(), assayName = Seurat::DefaultAssay(seuratObj), nCores = 1) {
  # Perform normalization once:
  seuratObj <- .DoNormalizationLda(seuratObj, normalizationMethod = normalizationMethod, varFeatures = varFeatures, assayName = assayName)

  ldaResults <- runLDA(seuratObj, ntopics = ntopics, normalizationMethod = normalizationMethod, seed.number = randomSeed, varFeatures = varFeatures, cores = nCores, assayName = assayName, skipNormalization = TRUE)
  if (length(ntopics) > 1) {
    print(LDAelbowPlot(ldaResults, seuratObj))
  }

  return(ldaResults)
}

.DoNormalizationLda <- function(seuratObj, normalizationMethod, varFeatures, assayName) {
  seuratObj <- Seurat::NormalizeData(seuratObj, assay = assayName, normalization.method = normalizationMethod)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, assay = assayName, nfeatures = varFeatures)

  return(seuratObj)
}

#' @title Runs LDA Model
#'
#' @description This function runs an LDA model on scRNA-seq expression data
#'
#' @param seuratObj Seurat object containing the data the model was created with.
#' @param ntopics Number of topics to be used in the model. If parallel == TRUE, a vector of topics to run should be inputted
#' @param alpha the value for alpha in the LDA model
#' @param beta the value for beta in the LDA model
#' @param varFeatures the number of variable features to use in the LDA model. The more features that are used, the slower the model will run and the more noise that will be introduced, but the model will be more complete in representing your entire dataset.
#' @param iterations the number of iterations used when learning the LDA model.
#' @param burnin number of iterations to run to allow the model to learn before calculating certain statistics. Models start at random points, so this allows model to get closer to the fit before certain statistics are calculated.
#' @param cores Number of cores to use, only applicable if parallel = TRUE
#' @param normalizationMethod Normalization method used by Seurat NormalizeData. Options are CLR, LogNormalize and RC.
#' @param assayName The name of the assay holding the source data
#' @param seed.number random integer to set seed
#' @param skipNormalization If true, the data are assumed to be pre-normalized. Both normalization and Seurat::FindVarialeFeatures() are skipped. Therefore the arguments normalizationMethod and varFeatures are ignored.
#'
#' @author TITAN
#' @references https://github.com/ohsu-cedar-comp-hub/TITAN
#'
#' @return LDA Model
#' @export
#'
#' @import Seurat
runLDA <- function(seuratObj,
                   ntopics,
                   alpha = 50,
                   beta = 0.1,
                   varFeatures = 5000,
                   iterations = 500,
                   burnin = 250,
                   seed.number = GetSeed(),
                   cores = 1,
                   normalizationMethod = "CLR",
                   assayName = Seurat::DefaultAssay(seuratObj),
                   skipNormalization = FALSE) {

  ## Set seed
  set.seed(seed.number)

  #Normalize and extract the gene expression data from the Seurat Object
  if (!skipNormalization) {
    seuratObj <- .DoNormalizationLda(seuratObj, normalizationMethod = normalizationMethod, varFeatures = varFeatures, assayName = assayName)
  }

  Object.sparse <- Seurat::GetAssayData(seuratObj, slot = "data",assay = assayName)
  Object.sparse <- Object.sparse[Seurat::VariableFeatures(seuratObj, assay = assayName),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
  data.use <- Matrix::Matrix(Object.sparse, sparse = T)

  data.use      <- data.use * 10
  data.use      <- round(data.use)
  data.use      <- Matrix::Matrix(data.use, sparse = T)
  sumMat        <- Matrix::summary(data.use)
  cellList      <- split(as.integer(data.use@i), sumMat$j)
  ValueList     <- split(as.integer(sumMat$x), sumMat$j)
  cellList      <- mapply(rbind, cellList, ValueList, SIMPLIFY=F)
  Genes         <- rownames(data.use)
  cellList      <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1,]+1];x})

  #Run model
  model_maker <- function(topics) {
    print(paste0('Running LDA with ', topics, ' topics'))
    selected.Model <- lda::lda.collapsed.gibbs.sampler(
      cellList,
      topics,
      Genes,
      num.iterations=iterations,
      alpha=alpha,
      eta=beta,
      compute.log.likelihood=TRUE,
      burnin=burnin)[-1]

    return(selected.Model)
  }

  if (length(ntopics) > 1) {
    models <- parallel::mclapply(ntopics, model_maker, mc.cores = cores)
    return(models)
  } else {
    m <- model_maker(ntopics)
    return(m)
  }
}


LDAelbowPlot <- function(ldaResults, seuratObj, assayName = "RNA") {
  Object.sparse <- GetAssayData(seuratObj, slot = "data",assay = assayName)
  Object.sparse <- Object.sparse[VariableFeatures(seuratObj, assay = assayName),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
  data.use <- Matrix::Matrix(Object.sparse, sparse = T)
  data.use <- data.use * 10
  data.use <- round(data.use)

  #initialize necessary variables
  perp_list     <- NULL
  topic_numbers <- NULL
  RPC           <- NULL

  for (topic_num in names(ldaResults)) {
    model <- ldaResults[[topic_num]]
    topic_numbers <- c(topic_numbers, topic_num)

    #extract document-term matrix
    docterMat     <- t(as.matrix(data.use))
    docterMat     <- methods::as(docterMat, "sparseMatrix")

    #calculate topic word distribution
    topworddist   <- .normalize(model$topics, byrow = T)

    #calculate document topic distribution
    doctopdist    <- .normalize(t(model$document_sums), byrow = T)

    #calculate perpelexity
    perp          <- text2vec::perplexity(docterMat, topworddist, doctopdist)
    perp_list     <- c(perp_list, perp)

    #calculate RPC (rate of perplexity change)
    if (length(perp_list) > 1) {
      RPC_temp <- abs((perp_list[length(perp_list)] - perp_list[length(perp_list) - 1]) / (topic_numbers[length(topic_numbers)] - topic_numbers[length(topic_numbers) - 1]))
      RPC <- c(RPC, RPC_temp)
    }
  }

  #build plot dataframe and create ggplot object
  plot_df           <- as.data.frame(cbind(topic_numbers[-1], RPC))
  colnames(plot_df) <- c("Topics", "RPC")
  p                 <- ggplot(data = plot_df, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point()

  return(p)
}

.normalize <- function(x, byrow = TRUE, tol = 1e-6) {
  object <- x

  if (!byrow){
    object.new <- t(.normalize( t(object), tol = tol, byrow = TRUE) )
  } else {
    if (is.matrix(object) || any(is(object) == "Matrix")) {
      object.new <- .threshold(object, min = 0)
      max.pos <- integer(0)
      all.zeros <- which(apply(object.new, 1, function(u) all(u == 0)))
      if (any(all.zeros)) {
        if (length(all.zeros) > 1){
          max.pos <- apply(object[all.zeros,], 1, which.max)
        } else if (length(all.zeros) == 1) {
          max.pos <- which.max(object[all.zeros,])
        }
        object.new[cbind(all.zeros, max.pos)] <- 1
      }
      if (any(is(object.new) == "Matrix")) {
        # normalize rows for sparse matrices
        object.new <- Diagonal(x = 1 / rowSums(object.new)) %*% object.new
      } else {
        object.new <- sweep(object.new, 1, rowSums(object.new), "/")
      }
      if (tol > 0) {
        object.new[object.new < tol] <- 0
        if (any(is(object.new) == "Matrix")) {
          # normalize rows for sparse matrices
          object.new <- Diagonal(x = 1 / rowSums(object.new)) %*% object.new
        } else {
          object.new <- sweep(object.new, 1, rowSums(object.new), "/")
        }
      }
    } else if (is.vector(object)) {
      max.pos <- which.max(object)
      object.new <- .threshold(object, min = 0)
      if (all(object.new == 0)) {
        object.new[max.pos] <- 1
      } else {
        object.new <- object.new / sum(object.new)
      }
    }
  }

  return(object.new)
}

.threshold <- function(x, min = -Inf, max = Inf){
  if (min > -Inf){
    x[x < min] <- min
  }
  if (max < Inf){
    x[x > max] <- max
  }
  invisible(x)
}