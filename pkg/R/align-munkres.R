## munkres (hungarian) algorithm

## based off of https://www.r-bloggers.com/munkres-assignment-algorithm-with-rcpparmadillo/

## STEP ONE: subtract minimal element of each row from each element in this row
step.one <- function(step, cost, N) {
  sapply(length(N), function(r) { 
    min <- min(cost[r,])
    sapply(cost[r,], function(y) {
      y <- y - min
    })
  })
  step <- 2
}


## STEP TWO: search for a zero in the modified cost matrix from STEP ONE
# only taking the first zero in a row and column
# indM = indicator matrix indicates this zero by setting the corresponding element at (x,y) to 1
step.two <- function(step, cost, indM, rcov, ccov, N) {
  sapply(length(N), function(r) {
    sapply(length(N), function(c) {
      if (cost[r,c] == 0.0 && rcov[r] == 0 && ccov[c] == 0) {
        indM[r,c] <- 1
        rcov[r] <- 1
        ccov[c] <- 1
        break
      }
    })
  })
  # for reuse later (wipe clean)
  rcov <- as.vector(rep(c(0), length=length(rcov)))
  ccov <- as.vector(rep(c(0), length=length(ccov)))
  step <- 3
}


## STEP THREE: cover each column w/ a *starred zero*
# do this by looking for 1's in indM
# if N columns are covered, all *starred zero*s describe a complete assignment --> go to STEP SEVEN and finish
# otherwise, go to STEP FOUR
step.three <- function(step, indM, ccov, N) {
  colcount <- 0
  sapply(length(N), function(r) {
    sapply(length(N), function(c) {
      if (indM[r,c] == 1) {
        ccov[c] = 1
        colcount <- colcount + 1
      }
    })
  })
  if (colcount == N) {
    step = 7
  } else {
    step = 4
  }
}


## STEP FOUR: find non-covered zeros and prime them
# if there are zeros in a row and none of them are *starred*, /prime/ them
# requires 3 helper functions: find.noncovered.zero(), star.in.row(), find.star.in.row()
step.four <- function(stpe, cost, indM, rcov, ccov, rpath_0, cpath_0, N) {
  row = col <- -1
  done <- FALSE
  while (!done) {
    search <- find.noncovered.zero(row, col, cost, rcov, ccov, N)
    
    if (search[1] == -1) {
      done <- TRUE
      step <- 6
    } else {
      # uncovered zero
      indM[row,col] <- 2
      if (star.in.row(row, indM, N)) {
        col.num <- find.star.in.row(row, col, indM, N)
        # cover the row w/ *starred zero* and uncover column w/ *starred zero*
        rcov[row] <- 1
        ccov[col] <- 0
      } else {
        # no *starred zero* in row w/ uncovered zero
        done <- TRUE
        step <- 5
        rpath_0 <- row
        cpath_0 <- col
      }
    }
  }
}

# helper function to search for non-covered zeros
find.noncovered.zero <- function(row, col, cost, rcov, ccov, N) {
  r = 0
  done = FALSE
  row = col = -1
  while (!done) {
    c = 0
    while (TRUE) {
      if (cost[r,c] == 0.0 && rcov[r] == 0 && ccov[c] == 0) {
        row = r
        col = c
        done = TRUE
      }
      c <- c + 1
      if (c == N || done) {
        break
      }
    }
    r <- r + 1
    if (r == N) {
      done = TRUE
    }
  }
  return(c(row, col)) #inserted
}

# if no uncovered zero found in STEP FOUR --> go to STEP SIX
# if uncovered zero found, we set position in indM to 2
# search for a *starred zero* in row w/ the uncovered zero
# uncover column w/ the starred zero
# cover row w/ the starred zero
# two helper functions are required to do this
star.in.row <- function(row, indM, N) {
  tmp <- FALSE
  sapply(length(N), function(c) {
    if (indM[row,c] == 1) {
      tmp <- TRUE
      break
    }
  })
  return(tmp)
}

find.star.in.row <- function(row, col, indM, N) {
  col <- -1
  sapply(length(N), function(c) {
    if (indM[row,c] == 1) {
      col = c
    }
  })
  return(col) #inserted
}







































