require(ape, quietly=TRUE)

# create test fixtures
t1 <- read.tree(text="(A:0.1,B:0.2):0;")
t2 <- read.tree(text="((A:0.1,B:0.2):0.1,C:0.3):0;")
t3 <- read.tree(text="(C:0.3,D:0.4):0;")

# ultrametric trees
t4 <- read.tree(text="(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0;")
t5 <- read.tree(text="(((A:0.1,B:0.1):0.15,C:0.25):0.05,D:0.3):0;")

# Trees used in Kuhner & Yamato, 2014
t6 <- read.tree(text="((A:3, B:3):7, ((C:2, D:2):5, E:7):3):0;")
t7 <- read.tree(text="(((A:2.5, B:2.5):2, C:4.5):3, (D:1, E:1):6.5):0;")
