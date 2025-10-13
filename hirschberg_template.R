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
        for (x in 1:length(Y)) { # for each character in Y
            first_alignment_row <- append(first_alignment_row, "-")# add gap
            second_alignment_row <- append(second_alignment_row, Y[i]) # add character from Y
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if (length(Y) == 0) { # length of Y is equal to zero
        for (x in 1:length(Y)) { # for each character in X
            first_alignment_row <- append(first_alignment_row, X[i])# add character from X
            second_alignment_row <- append(second_alignment_row, "-")# add gap
        }
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else if (length(Y) == 1 && length(X) == 1) { # length of X and Y is equal to 1
        first_alignment_row <- append(first_alignment_row, X[i])# add character from X
        second_alignment_row <- append(second_alignment_row, Y[i])# add character from Y
        alignment <- c(DNAStringSet(first_alignment_row), DNAStringSet(second_alignment_row))
        print(alignment)
    } else {
        xlen <- length(X)# length of X
        xmid <- ceiling(xlen/2)# half of the length of X
        ylen <- length(Y)# length of Y

        # NW score for the first half of X and the whole Y
        first_score <- Needle_wunch(X[1:xmid], Y)
        # NW score for the second half of X and the whole Y (both are reversed)
        second_score <- Needle_wunch(rev(X[xmid:xlen]), rev(Y))

        ymid <- which.max(first_score-second_score)# index of division for Y

        # The first half
        if (ymid == 0) { # index of division for Y is equal to 0
            # call Hirschberg function for the first half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(X[1:xmid], DNAString(""), alignment, match, missmatch, gap)
        } else {
            # call Hirschberg function for the first half of X and for the first part of Y
            alignment <- HirschbergTemplate(X[1:xmid], Y[1:ymid], alignment, match, missmatch, gap)
        }
        
        # The second half
        if ((xmid + 1) > xlen) { # X cannot be further divided
            # call Hirschberg function for an empty DNAString object and the second half of Y
            alignment <- HirschbergTemplate(DNAString(""), Y[ymid:ylen], alignment, match, missmatch, gap)
        } else if ((ymid + 1) > ylen) { # Y cannot be further divided
            # call Hirschberg function for the second half of X and for an empty DNAString object
            alignment <- HirschbergTemplate(X[xmid:xlen], DNAString(""), alignment, match, missmatch, gap)
        } else {
            # call Hirschberg function for the second half of X and the second part of Y
            alignment <- HirschbergTemplate(X[xmid:xlen], Y[ymid:ylen], alignment, match, missmatch, gap)
        }
    }
    return(alignment)
}


Needle_wunch <- function(X, Y, match, missmatch, gap){
  xlen <- length(X)
  ylen <- length(Y)
  alignment <- matrix(c(1:ylen*xlen), nrow = xlen, ncol = ylen)
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
        alignment[i, j] <- max(c((alignment[i-1, j-1] + plus), (alignment[i-1, j] + gap), (alignment[i, j-1] + gap)))
      }
    }
  }
  return(alignment[xlen, 1:ylen])
}

X <- DNAString("AGTACGCA")
Y <- DNAString("TATGC")
match <- 3
missmatch <- 1
gap <- -3
alignment <- Needle_wunch(X, Y, match, missmatch, gap)
print(which.max(Needle_wunch(X, Y, match, missmatch, gap)-rev(Needle_wunch(X, Y, match, missmatch, gap))))
alignment <- DNAString("A")
append(DNAString("A"),alignment)
print(alignment)
alignment <- HirschbergTemplate(X, Y, alignment, match, mismatch, gap)