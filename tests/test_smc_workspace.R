require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.rescale.tree <- function() {
    # NB: root branch length is ignored!
    tree <- read.tree(text='(A:2,B:6):0.1;')
    config <- list(norm.mode="MEAN")
    result <- parse.input.tree(tree, config)
    expected <- read.tree(text='(A:0.5,B:1.5):0.1;')
    checkEquals(expected, result)
}

test.parse.input.tree <- function() {
    # check ladderization
    config <- list(norm.mode="NONE")
    result <- parse.input.tree('(s4:2,(s3:2,(s1:2,s2:2):2):2):2;', config)
    expected <- read.tree(text='(((s1:2,s2:2):2,s3:2):2,s4:2):2;')
    checkEquals(expected, result)
}

test.get.tip.heights <- function() {
    result <- get.tip.heights(t2)
    # ((A:0.1,B:0.2):0.1,C:0.3):0;
    expected <- c(0.1, 0, 0)
    checkEquals(expected, result)
}

test.parse.labels <- function() {
    labels <- c('1_A', '2_B', '3_A', '4_')
    regex <- '_([A-Z]*)$'
    result <- parse.labels(labels, regex)
    expected <- c('A', 'B', 'A', '')
    checkEquals(expected, result)
}

test.init.workspace <- function() {
    config <- list(
        params=c('N'),
        nparticle=10,
        nsample=3,
        norm.mode='NONE'
    )
    class(config) <- 'smc.config'
    result <- init.workspace(t1, config)
    checkEquals(result$n.tips, 2)
    checkEquals(result$tip.heights,  c(0.1, 0.0))
    checkEquals(result$tip.labels, NA)
    checkEquals(nrow(result$particles), 10)
    checkEquals(ncol(result$particles), 1)
}

