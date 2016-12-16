# This file is part of Kaphi.

# Kaphi is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Kaphi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Kaphi.  If not, see <http://www.gnu.org/licenses/>.

require(Kaphi, quietly=TRUE)
require(RUnit, quietly=TRUE)
source('tests/fixtures/simple-trees.R')


test.rescale.tree <- function() {
    # NB: root branch length is ignored!
    tree <- read.tree(text='(A:2,B:6):0.1;')
    config <- list(decay.factor=0.1, rbf.variance=2.0, sst.control=1, norm.mode="MEAN")
    result <- parse.input.tree(tree, config)
    expected <- read.tree(text='(A:0.5,B:1.5):0.1;')

    # compare serializations to ignore attributes produced by parser
    checkEquals(expected, result)
}

test.parse.input.tree <- function() {
    # check ladderization
    config <- list(decay.factor=0.1, rbf.variance=2.0, sst.control=1, norm.mode="NONE")
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
        decay.factor=1.0,
        rbf.variance=2.0,
        sst.control=1,
        norm.mode='NONE'
    )
    class(config) <- 'smc.config'
    result <- init.workspace(t1, config)

    checkEquals(result$n.tips, 2)
    checkEquals(result$tip.heights,  c(0.1, 0.0))
    checkEquals(result$tip.labels, NA)
    checkEquals(nrow(result$particles), 10)
    checkEquals(ncol(result$particles), 1)

    checkEquals(Ntip(result$obs.tree), 2)
    checkEquals(result$obs.tree$kernel, 4)
}

