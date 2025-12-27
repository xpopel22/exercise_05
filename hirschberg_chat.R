# hirshberg algoritmus od chatgpt

# needleman_wunsh pro male pripady. zakladni pripady hirshberga
NW_align <- function(X, Y, match, mismatch, gap) {
  xlen <- length(X)
  ylen <- length(Y)
  
  score <- matrix(0, nrow = xlen + 1, ncol = ylen + 1)
  trace <- matrix("", nrow = xlen + 1, ncol = ylen + 1)
  
  score[1, ] <- 0:ylen * gap
  score[, 1] <- 0:xlen * gap
  trace[1, ] <- "L"
  trace[, 1] <- "U"
  trace[1,1] <- ""
  
  for (i in 2:(xlen + 1)) {
    for (j in 2:(ylen + 1)) {
      s <- ifelse(X[i-1] == Y[j-1], match, mismatch)
      
      diag <- score[i-1, j-1] + s
      up   <- score[i-1, j] + gap
      left <- score[i, j-1] + gap
      
      score[i, j] <- max(diag, up, left)
      trace[i, j] <- c("D","U","L")[which.max(c(diag, up, left))]
    }
  }
  
  # backtracking
  i <- xlen + 1
  j <- ylen + 1
  ax <- ay <- character()
  
  while (i > 1 || j > 1) {
    if (trace[i, j] == "D") {
      ax <- c(X[i-1], ax)
      ay <- c(Y[j-1], ay)
      i <- i - 1; j <- j - 1
    } else if (trace[i, j] == "U") {
      ax <- c(X[i-1], ax)
      ay <- c("-", ay)
      i <- i - 1
    } else {
      ax <- c("-", ax)
      ay <- c(Y[j-1], ay)
      j <- j - 1
    }
  }
  
  list(ax, ay)
}


# needle mann wunch jen pro poslední řádek
NW_score <- function(X, Y, match, mismatch, gap) {
  xlen <- length(X)
  ylen <- length(Y)
  
  score <- matrix(0, nrow = xlen + 1, ncol = ylen + 1)
  score[1, ] <- 0:ylen * gap
  score[, 1] <- 0:xlen * gap
  
  for (i in 2:(xlen + 1)) {
    for (j in 2:(ylen + 1)) {
      s <- ifelse(X[i-1] == Y[j-1], match, mismatch)
      score[i, j] <- max(
        score[i-1, j-1] + s,
        score[i-1, j] + gap,
        score[i, j-1] + gap
      )
    }
  }
  
  score[xlen + 1, ]
}

# hirshberg pouze se sekvenci s characters
Hirschberg <- function(X, Y, match, mismatch, gap) {
  
  if (length(X) == 0)
    return(list(rep("-", length(Y)), Y))
  
  if (length(Y) == 0)
    return(list(X, rep("-", length(X))))
  
  if (length(X) == 1 || length(Y) == 1)
    return(NW_align(X, Y, match, mismatch, gap))
  
  xmid <- floor(length(X) / 2)
  
  scoreL <- NW_score(X[1:xmid], Y, match, mismatch, gap)
  scoreR <- NW_score(rev(X[(xmid+1):length(X)]),
                     rev(Y), match, mismatch, gap)
  
  ymid <- which.max(scoreL + rev(scoreR)) - 1
  
  left  <- Hirschberg(X[1:xmid], Y[1:ymid], match, mismatch, gap)
  right <- Hirschberg(X[(xmid+1):length(X)], Y[(ymid+1):length(Y)],
                      match, mismatch, gap)
  
  list(
    c(left[[1]], right[[1]]),
    c(left[[2]], right[[2]])
  )
}

library(Biostrings)

toDNAAlignment <- function(alignment) {
  DNAStringSet(c(
    paste(alignment[[1]], collapse = ""),
    paste(alignment[[2]], collapse = "")
  ))
}


X <- strsplit("TACGAGGCA", "")[[1]]
Y <- strsplit("ACGGA", "")[[1]]

res <- Hirschberg(X, Y, match = 1, mismatch = -1, gap = -2)

paste(res[[1]], collapse = "")
paste(res[[2]], collapse = "")

dna_alignment <- toDNAAlignment(res)
dna_alignment


