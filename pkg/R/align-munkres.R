## munkres (hungarian) algorithm

## based off of https://www.r-bloggers.com/munkres-assignment-algorithm-with-rcpparmadillo/

## STEP ONE: subtract minimal element of each row from each element in this row
step.one <- function(step, cost, N) {
  for (r in 1:N) { 
    min <- min(cost[r,])
    for (c in 1:N) {
      cost[r,c] <- cost[r,c] - min
    }
  }
  step <- 2
  return(list(step, cost))
}


## STEP TWO: search for a zero in the modified cost matrix from STEP ONE
# only taking the first zero in a row and column
# indM = indicator matrix indicates this zero by setting the corresponding element at (x,y) to 1
step.two <- function(step, cost, indM, rcov, ccov, N) {
  for (r in 1:N) {
    for (c in 1:N) {
      if (cost[r,c] == 0.0 && rcov[r] == 0 && ccov[c] == 0) {
        indM[r,c] <- 1
        rcov[r] <- 1
        ccov[c] <- 1
        break
      }
    }
  }
  step <- 3
  result <- list(step, indM)
  return(result)
}


## STEP THREE: cover each column w/ a *starred zero*
# do this by looking for 1's in indM
# if N columns are covered, all *starred zero*s describe a complete assignment --> go to STEP SEVEN and finish
# otherwise, go to STEP FOUR
step.three <- function(step, indM, ccov, N) {
  colcount <- 0
  for (r in 1:N) {
    for (c in 1:N) {
      if (indM[r,c] == 1) {
        ccov[c] = 1
        colcount <- colcount + 1
      }
    }
  }
  if (colcount == N) {
    step = 7
  } else {
    step = 4
  }
  return(list(step, ccov))
}


## STEP FOUR: find non-covered zeros and prime them
# if there are zeros in a row and none of them are *starred*, /prime/ them
# requires 3 helper functions: .find.noncovered.zero(), .star.in.row(), .find.star.in.row()
step.four <- function(step, cost, indM, rcov, ccov, rpath_0, cpath_0, N) {
  row = col <- -1
  done <- FALSE
  while (!done) {
    search <- .find.noncovered.zero(row, col, cost, rcov, ccov, N)
    row <- search[1]
    col <- search[2]
    if (row == -1) {
      done <- TRUE
      step <- 6
    } else {
      # uncovered zero
      indM[row,col] <- 2
      if (.star.in.row(row, indM, N)) {
        col.num <- .find.star.in.row(row, col, indM, N)
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
  return(list(step, indM, rcov, ccov, rpath_0, cpath_0))
}

# helper function to search for non-covered zeros
.find.noncovered.zero <- function(row, col, cost, rcov, ccov, N) {
  r <- 1
  done <- FALSE
  row = col <- -1
  while (!done) {
    c <- 1
    while (TRUE) {
      if (cost[r,c] == 0.0 && rcov[r] == 0 && ccov[c] == 0) {
        row <- r
        col <- c
        done <- TRUE
      }
      c <- c + 1
      if (c > N || done) {
        break
      }
    }
    r <- r + 1
    if (r > N) {
      done <- TRUE
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
.star.in.row <- function(row, indM, N) {
  tmp <- FALSE
  for (c in 1:N) {
    if (indM[row,c] == 1) {
      tmp <- TRUE
      break
    }
  }
  return(tmp)
}

.find.star.in.row <- function(row, col, indM, N) {
  col <- -1
  for (c in 1:N) {
    if (indM[row,c] == 1) {
      col <- c
    }
  }
  return(col) #inserted
}


## STEP FIVE: constructs path beginning at an uncovered, /primed zero/, and alternating btwn *starred* and /primed/ zeros
# path is continued until a /primed zero/ w/ NO *starred zero* in its column is found
# then all *starred zeros* are unstarred
# then all /primed zeros/ are *starred*
# all primes in the indM are erased and all rows are uncovered
# then return to STEP THREE to cover over the columns again
# requires 4 helper functions: ..find.star.in.col(), ..find.prime.in.row(), ..augment.path(), .erase.primes()
step.five <- function(step, indM, rcov, ccov, path, rpath_0, cpath_0, N) {
  done <- FALSE
  row = col <- -1
  path_count <- 2
  path[path_count-1, 1] <- rpath_0
  path[path_count-1, 2] <- cpath_0
  while (!done) {
    row <- .find.star.in.col(path[path_count-1, 2], row, indM, N)
    if (row > -1) {
      # *starred zero* in row "row"
      path_count <- path_count + 1
      path[path_count-1, 1] <- row
      path[path_count-1, 2] <- path[path_count-2, 2]
    } else {
      done <- TRUE
    }
    if (!done) {
      # if there is a *starred zero*, find a /primed zero/ in this row
      # write index to "col"
      col <- .find.prime.in.row(path[path_count-1, 1], col, inM, N)
      path_count <- path_count + 1
      path[path_count-1, 1] <- path[path_count-2, 1]
      path[path_count-1, 2] <- col
    }
  }
  indM <- .augment.path(path_count, indM, path)
  rcov <- as.vector(rep(c(0), length=length(rcov)))    # clear the covers from rows
  ccov <- as.vector(rep(c(0), length=length(ccov)))
  indM <- .erase.primes(indM, N)
  step <- 3
  return(list(step, indM, rpath_0, cpath_0))
}

# helper function to find *starred zeros* in columns
.find.star.in.col <- function(col, row, indM, N) {
  row <- -1
  for (r in 1:N) {
    if (indM[r,col] == 1) {
      row <- r
    }
  }
  return(row)
}

# helper function to find a primed zero in a row
.find.prime.in.row <- function(row, col, indM, N) {
  for (c in 1:N) {
    if (indM[row,c] == 2) {
      col <- c
    }
  }
  return(col)
}

# helper function to augment the path
.augment.path <- function(path_count, indM, path) {
  for(p in 1:(path_count-1)) {
    if (indM[ path[p,1], path[p,2] ] == 1) {
      indM[ path[p,1], path[p,2] ] <- 0  
    } else {
      indM[ path[p,1], path[p,2] ] <- 1
    }
  }
  return(indM)
}

# helper function to erase /primed zeros/ from indM
.erase.primes <- function(indM, N) {
  for (r in 1:N) {
    for (c in 1:N) {
      if (indM[r,c] == 2) {
        indM[r,c] <- 0
      }
    }
  }
  return(indM)
}


# STEP SIX: take cover vectors 'rcov' and 'ccov' and look into uncovered region of cost matrix for smallest value
# subtract this value from each element in an uncovered column and add it to each element in a covered row
# requires 1 helper function: ..find.smallest()
step.six <- function(step, cost, rcov, ccov, N) {
  minval <- .Machine$double.xmax
  smallest <- .find.smallest(minval, cost, rcov, ccov, N)
  for (r in 1:N) {
    for (c in 1:N) {
      if (rcov[r] == 1) {
        cost[r,c] <- cost[r,c] + smallest
      }
      if (ccov[c] == 0) {
        cost[r,c] <- cost[r,c] - smallest
      }
    }
  }
  step <- 4
  return(list(step, cost))
}

# helper function searches for smallest value in uncovered region of the cost matrix
.find.smallest <- function(minval, cost, rcov, ccov, N) {
  for (r in 1:N) {
    for (c in 1:N) {
      if (rcov[r] == 0 && ccov[c] == 0) {
        if (minval > cost[r,c]) {
          minval <- cost[r,c]
        }
      }
    }
  }
  return(minval)
}


#----------MAIN FUNCTION---------#
hungarian.alg <- function(input_cost) {
  N <- nrow(input_cost)
  step <- 1
  cpath_0 = rpath_0 <- 0
  cost <- input_cost
  indM <- matrix(data=0, nrow=N, ncol=N)
  rcov <- as.vector(rep(c(0), length=N))
  ccov <- as.vector(rep(c(0), length=N))
  path <- matrix(nrow=N*2, ncol=2)    #changed..? originally -->  arma::imat path(2 * N, 2);
  
  done <- FALSE
  while (!done) {
    if (step == 1) {
      res1 <- step.one(step, cost, N)
      step <- res1[[1]]
      cost <- res1[[2]]
    } else if (step == 2) {
      res2 <- step.two(step, cost, indM, rcov, ccov, N)
      step <- res2[[1]]
      indM <- res2[[2]]
    } else if (step == 3) {
      res3 <- step.three(step, indM, ccov, N)
      step <- res3[[1]]
      ccov <- res3[[2]] 
    } else if (step == 4) {
      res4 <- step.four(step, cost, indM, rcov, ccov, rpath_0, cpath_0, N)
      step <- res4[[1]]
      indM <- res4[[2]]
      rcov <- res4[[3]]
      ccov <- res4[[4]]
      rpath_0 <- res4[[5]]
      cpath_0 <- res4[[6]]
    } else if (step == 5) {
      res5 <- step.five(step, indM, rcov, ccov, path, rpath_0, cpath_0, N)
      step <- res5[[1]]
      indM <- res5[[2]]
      rpath_0 <- res5[[3]]
      cpath_0 <- res5[[4]]
    } else if (step == 6) {
      res6 <- step.six(step, cost, rcov, ccov, N)
      step <- res6[[1]]
      cost <- res6[[2]]
    } else if (step == 7) {
      done <- TRUE
    }
  }
  return(indM)
}



