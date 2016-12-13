require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')

test.rescale.tree <- function() {
    expected <- 0.15
    result <- mean(temp$edge.length)
    checkEquals(expected, result)

    temp <- rescale.tree(t1, "MEAN")
    expected <- 1.0
    result <- mean(temp$edge.length)
    checkEquals(expected, result)
}


test.kernel.trivial <- function() {
    expected <- 1
    result <- tree.kernel(t1, t1, lambda=0.1, sigma=2.0, rho=1, normalize=1, rescale.mode="NONE")
    # normalized kernel score of a tree against itself should always be 1
    checkEquals(expected, result)
}

test.kernel.unnormalized <- function() {
    # t1 has one internal node (n1) with two tips
    # delta(n1,n1) = lambda * rho * [(1 + delta(cn1, cn1)) * (1 + delta(cn2, cn2))]
    # delta(cn1, cn1) = lambda (tips)
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 1 * 1 * ((1+1) * (1+1))
    checkEquals(expected, result)

    result <- tree.kernel(t1, t1, lambda=0.5, sigma=1, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 0.5 * 1 * ((1+0.5) * (1+0.5))
    checkEquals(expected, result)

    # RBF = 1 for all sigma because branch lengths are identical
    result <- tree.kernel(t1, t1, lambda=1.0, sigma=0.5, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)

    # t1 is nested within t2
    result <- tree.kernel(t1, t2, lambda=1.0, sigma=Inf, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)

    # t3 has the same topology as t1 but different branch lengths
    result <- tree.kernel(t1, t3, lambda=1.0, sigma=Inf, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4  # no RBF penalty when sigma -> Infinity
    checkEquals(expected, result)

    # note sigma parameter for RBF *includes* the factor of two!
    result <- tree.kernel(t1, t3, lambda=1.0, sigma=1.0, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 1 * exp(-((0.3-0.1)^2 + (0.4-0.2)^2) / (1.0)) * ((1+1) * (1+1))
    checkEquals(expected, result)
}


test.kernel.normalized <- function() {
    result <- tree.kernel(t1, t1, lambda=1, sigma=0.5, rho=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
}


test.kernel.label <- function() {
	# essentially unlabelled
	cat("essentially unlabelled\n")
	result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=t1$tip.label, label2=t1$tip.label, gamma=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
        
    # same tree labelled
    cat("same tree labelled\n")
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=t1$tip.label, label2=t1$tip.label, gamma=0, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
    
    # same tree labels switched
    cat("same tree labels switched\n")
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=t1$tip.label, label2=t1$tip.label[2:1], gamma=0, normalize=0, rescale.mode='NONE')
    expected <- exp(-(.1^2+.1^2)/(1.0)) * ((1+1) * (1+1))
    checkEquals(expected, result)
    
    # same tree different labels
    cat("same tree different labels\n")
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=c(1, 2), label2=c(3, 4), gamma=0, normalize=0, rescale.mode='NONE')
    expected <- 1
    checkEquals(expected, result)
    
    # half gamma
    cat("half gamma\n")
	result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=c(1, 2), label2=c(3, 4), gamma=.5, normalize=0, rescale.mode='NONE')
    expected <- 1.5 * 1.5
    checkEquals(expected, result)
    
    # same tree different labels, ignore mismatch
    cat("same tree different labels, ignore mismatch\n")
    result <- tree.kernel(t1, t1, lambda=1, sigma=1, rho=1, label1=c(1, 2), label2=c(3, 4), gamma=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
    
    # subset tree
    cat("subset tree\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=t1$tip.label, label2=t2$tip.label, gamma=0, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
    
    # subset tree, half gamma
    cat("subset tree, half gamma\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=t1$tip.label, label2=t2$tip.label, gamma=.5, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
    
    # subset tree, labels different
    cat("subset tree, labels different\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=c(1, 1), label2=c(2, 2, 2), gamma=0, normalize=0, rescale.mode='NONE')
    expected <- 1
    checkEquals(expected, result)
    
    # subset tree, ignore mismatch
    cat("subset tree, ignore mismatch\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=c(1, 1), label2=c(2, 2, 2), gamma=0.5, normalize=0, rescale.mode='NONE')
    expected <- 1.5 * 1.5
    checkEquals(expected, result)
    
    # subset tree, ignore mismatch
    cat("subset tree, ignore mismatch\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=c(1, 1), label2=c(2, 2, 2), gamma=1, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
    
    # subset tree, outgroup a
    cat("subset tree, outgroup a\n")
    result <- tree.kernel(t1, t2, lambda=1, sigma=1, rho=1, label1=c(1, 2), label2=c(1, 2, 1), gamma=0, normalize=0, rescale.mode='NONE')
    expected <- 4
    checkEquals(expected, result)
}