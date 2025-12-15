#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param X DNAString object representing NT or AA sequence to be aligned.
#' @param Y DNAString object representing NT or AA sequence to be aligned.
#' @param alignment A list of DNAString objects with alignment of input sequences.
#' @param match An integer value of a score for matching bases.
#' @param mismatch An integer value of a score for mismatching bases.
#' @param gap An integer value of a penalty for gap insertion.
#' @returns A list of DNAString objects with alignment of input sequences.
HirschbergTemplate <- function(X, Y, alignment, match, mismatch, gap){
    
    first_alignment_row <- alignment[[1]] # initialize the first row of alignment
    second_alignment_row <- alignment[[2]] # initialize the second row of alignment
  
    if (length(X) == 0) { # length of X is equal to zero
        for (i in 1:length(Y)) { # for each character in Y
            first_alignment_row <- append(first_alignment_row, DNAString("-"))# add gap
            second_alignment_row <- append(second_alignment_row, Y[i]) # add character from Y
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        #print(alignment)
    } else if (length(Y) == 0) { # length of Y is equal to zero
        for (i in 1:length(X)) { # for each character in X
            first_alignment_row <- append(first_alignment_row, X[i])# add character from X
            second_alignment_row <- append(second_alignment_row, DNAString("-"))# add gap
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        #print(alignment)
    } else if (length(Y) == 1 && length(X) == 1) { # length of X and Y is equal to 1
        first_alignment_row <- append(first_alignment_row, X[1])# add character from X
        second_alignment_row <- append(second_alignment_row, Y[1])# add character from Y
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        #print(alignment)
    } else {
        xlen <- length(X)# length of X
        xmid <- ceiling(xlen/2)# half of the length of X
        ylen <- length(Y)# length of Y

        # NW score for the first half of X and the whole Y
        first_score <- Needle_wunch(X[1:xmid], Y, match, mismatch, gap)
        # NW score for the second half of X and the whole Y (both are reversed)
        second_score <- Needle_wunch(rev(X[xmid:xlen]), rev(Y), match, mismatch, gap)
        

        ymid <- which.max(first_score+rev(second_score)) - 1# index of division for Y

        # The first half
        if (ymid == 0) { # index of division for Y is equal to 0
            # call Hirschberg function for the first half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(X[1:xmid], DNAString(""), alignment, match, mismatch, gap)
        } else {
            # call Hirschberg function for the first half of X and for the first part of Y
            alignment <- HirschbergTemplate(X[1:xmid], Y[1:ymid], alignment, match, mismatch, gap)
        }
        
        # The second half
        if ((xmid + 1) > xlen) { # X cannot be further divided
            # call Hirschberg function for an empty DNAString object and the second half of Y
            alignment <- HirschbergTemplate(DNAString(""), Y[ymid:ylen], alignment, match, mismatch, gap)
        } else if ((ymid + 1) > ylen) { # Y cannot be further divided
            # call Hirschberg function for the second half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(X[(xmid+1):xlen], DNAString(""), alignment, match, mismatch, gap)
        } else {
            # call Hirschberg function for the second half of X and the second part of Y
            alignment <- HirschbergTemplate(X[(xmid+1):xlen], Y[(ymid+1):ylen], alignment, match, mismatch, gap)
        }
    }
    #print(length(X))
    #print(length(Y))
    #alignment[[1]] = alignment[[1]][length(X):length(alignment[[1]])]
    #alignment[[2]] = alignment[[2]][length(Y):length(alignment[[2]])]
    return(alignment)
}


Needle_wunch <- function(X, Y, match, missmatch, gap){
  xlen <- length(X)
  ylen <- length(Y)
  alignment <- matrix(c(1:(ylen*xlen)), nrow = xlen, ncol = ylen)
  for (i in 1:xlen){
    for (j in 1:ylen){
      if (i == 1){
        alignment[i, j] <- (j-1)*gap
      } else if (j == 1) {
        alignment[i, j] <- (i-1)*gap
      }else {
        if (X[i] == Y[j]){
          plus <- match
        } else {
          plus <- missmatch
        }
        alignment[i, j] <- max(c((alignment[(i-1), (j-1)] + plus), (alignment[i-1, j] + gap), (alignment[i, j-1] + gap)))
      }
    }
  }
  return(alignment[xlen, 1:ylen])
}

library(Biostrings)

X <- DNAString("AGTACGCA")
Y <- DNAString("TATGC")
match <- 3
missmatch <- 1
gap <- -3
alignment <- Needle_wunch(X, Y, match, missmatch, gap)
#print(which.max(Needle_wunch(X, Y, match, missmatch, gap)-rev(Needle_wunch(X, Y, match, missmatch, gap))))
alignment <- list(DNAString(""), DNAString(""))
print(alignment)
append(alignment[[1]], DNAString("-"))
alignment <- HirschbergTemplate(X, Y, alignment, match, missmatch, gap)
print(alignment)
