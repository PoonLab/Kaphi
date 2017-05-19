require(ape, quietly=TRUE)

# create test fixtures
t1 <- read.tree(text="(A:0.1,B:0.2):0;")
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")
t3 <- read.tree(text="(C:0.3,D:0.4):0;")

# ultrametric trees
t4 <- read.tree(text="(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0;")
t5 <- read.tree(text="(((A:0.1,B:0.1):0.15,C:0.25):0.05,D:0.3):0;")

