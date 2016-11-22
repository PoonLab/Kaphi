require(ape)

# create test fixtures
t1 <- read.tree(text="(A:0.1,B:0.2):0;")
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")
t3 <- read.tree(text="(C:0.3,D:0.4):0;")
