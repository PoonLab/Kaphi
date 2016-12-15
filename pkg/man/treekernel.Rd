\name{tree.kernel}
\alias{tree.kernel}
\alias{rescale.tree}
\alias{preprocess.tree}
\alias{to.newick}
\alias{utk}

\title{Tree shape kernel functions}
\description{
    A kernel function maps objects into a high dimensional feature space where the inner product between objects can be computed efficiently.

    Kaphi uses a
    \code{tree.kernel} is a wrapper function that calls the associated function in the C extension \code{src/kernel.c} to compute a tree shape kernel function.

}
