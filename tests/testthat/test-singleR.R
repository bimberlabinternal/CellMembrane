context("scRNAseq")

test_that("SingleR works as expected", {
    set.seed(CellMembrane::GetSeed())

    seuratObj <- readRDS('../testdata/seuratOutput.rds')

    results <- 'singleR.txt'
    rawDataFile <- 'singleR.results.txt'

    nGene <- nrow(seuratObj)
    nCell <- ncol(seuratObj)
    datasets <- c('hpca', 'blueprint', 'dice', 'monaco')
    seuratObj <- RunSingleR(seuratObj = seuratObj, resultTableFile = results, rawDataFile = rawDataFile, datasets = datasets)

    nGene2 <- nrow(seuratObj)
    nCell2 <- ncol(seuratObj)

    expect_equal(nGene, nGene2)
    expect_equal(nCell, nCell2)

    print(table(seuratObj$hpca.label))
    print(table(seuratObj$hpca.label.fine))

    print(table(seuratObj$SingleRConsensus))
    expect_equal(1540, sum(seuratObj$SingleRConsensus == 'NK/T_cell', na.rm = T))

    allData <- read.table(rawDataFile, header = T, sep = '\t')
    print(nrow(allData))
    
    #TODO:
    unlink(rawDataFile)

    expect_equal(99, sum(seuratObj$hpca.label == 'NK_cell'))
    expect_equal(1404, sum(seuratObj$hpca.label == 'T_cells'))
    expect_equal(0, sum(seuratObj$hpca.label == 'B_cell'))
    expect_equal(0, sum(seuratObj$hpca.label == 'Neutrophils'))
    expect_equal(54, sum(seuratObj$hpca.label == 'Unknown'))

    expect_equal(150, sum(seuratObj$hpca.label.fine == 'T_cell:CD8+_Central_memory'))
    expect_equal(60, sum(seuratObj$hpca.label.fine == 'T_cell:CD8+_naive'))

    expect_equal(ncol(seuratObj), nrow(read.table(results, sep = '\t', header = T)))

    unlink(results)

    expect_equal(1252, sum(seuratObj$dice.label == 'NK cells'))
    expect_equal(228, sum(seuratObj$dice.label == 'T cells, CD4+'))
    expect_equal(3, sum(seuratObj$dice.label == 'Unknown'))

    for (dataset in datasets) {
        expect_true(!is.null(seuratObj@meta.data[[paste0(dataset,'.label')]]))
        expect_true(!is.null(seuratObj@meta.data[[paste0(dataset,'.label.fine')]]))
    }

    #saveRDS(seuratObj, file = '../testdata/singleR.rds')
})

