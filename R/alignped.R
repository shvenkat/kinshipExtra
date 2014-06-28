#' Generate a pedigree layout for plotting
#'
#' This function provides convenient access to layout algorithms in this
#' package as well as in kinship2.
#'
#' @inheritParams ggpedigree
#' @param method
#'      one of the strings "classic", "compact" and "1d"
#' @return
#'      pedigree alignment list, such as returned by kinship2::align.pedigree
#' @section Details:
#'      The layout methods available are: classic, compact and 1d. Classic
#'      gives the common layout used in textbooks, with parents centered on
#'      their children below. Compact generates a minimum-width layout, in
#'      which parents may not be directly above their children, leading to the
#'      use of diagonal descending connectors. 1d assigns non-overlapping
#'      horizontal positions to pedigree members across generations, producing
#'      a layout suitable for annotating the axis of another plot.
#'
#'      Note: In a 1d layout, each subtree of the pedigree is compact, in the
#'      sense that no member outside that subtree occupies a position within
#'      the interval spanned by the subtree.
#' @import kinship2
#' @export
alignped <- function(ped, method = "classic") {
    method <- match.arg(method, c("classic", "compact", "1d"))
    pedalign <- switch(method,
        "classic" = align.pedigree(ped, packed = FALSE, align = FALSE),
        "compact" = align.pedigree(ped, packed = TRUE, align = TRUE),
        "1d"      = align.pedigree.1d(ped),
        stop(sprintf("%s is not a valid 'method'", method)))
    return(pedalign)
}

#' Generate a "1-dimensional" layout of a pedigree
#'
#' Given a pedigree, this function returns a pedigree alignment that assigns
#' each person a unique horizontal position, making it suitable for plotting
#' persons in the pedigree along an axis, while depicting their relationships.
#' This is useful in annotating a figure (e.g. haplotype sharing) where the
#' same persons are plotted along one of the axes.
#'
#' In this alignment, each pedigree subtree is compact. A subtree of a pedigree
#' comprises a person, their partner, descendants and descendants' partners.
#' This layout generates compact subtrees, in that each one occupies a
#' horizontal interval _exclusively_. In other words, no persons outside a
#' subtree (e.g. ancestors, siblings, etc.) occupy the interval corresponding
#' to it. This layout is the most "natural" looking 1d layout.
#' @inheritParams alignped
#' @return
#'     pedigree alignment list, such as returned by kinship2::align.pedigree
#' @import kinship2
align.pedigree.1d <- function(ped) {

    pedalign <- alignped(ped, method = "classic")
    # If the pedigree has only one generation, no further work is needed
    if(length(pedalign$n) < 2)
        return(pedalign)

    n      <- pedalign$n
    nid    <- pedalign$nid
    fam    <- pedalign$fam
    spouse <- pedalign$spouse
    pos    <- pedalign$pos
    pos[nid == 0] <- NA         # NA indicates an empty position
    pos[!is.na(pos)] <- 0       # non-empty positions are zero-ed out

    # Process levels one at a time top-down, updating layout positions
    for(i in 1:length(n)) {
        # Process each level left-to-right in extended sibships
        # A extended sibship is a set of full siblings with their partners.
        # Non-founder partners are assumed to be depicted a second time on the
        # right with their ancentors and siblings, with the exception of the
        # rightmost non-founder partner.
        j <- 1
        while(j <= n[i]) {
            # Empty cell, move on.  This case should not occur.
            if(nid[i, j] == 0) {
                warning(sprintf(paste("Unexpected condition:",
                            "Empty position reached at (i, j) = (%i, %i)"),
                        i, j))
                j <- j + 1
                next
            }
            # Person is neither a descendant nor spouse.  This case shouldn't
            # occur.
            if(fam[i, j] == 0 && spouse[i, j] == 0) {
                stop(sprintf(paste("Unsupported pedigree:",
                            "fam = 0 and spouse = 0 for (i, j) = (%i, %i)"),
                        i, j))
                next
            }
            # Pair of founders, insert them at the rightmost position of the
            # current level.
            if(fam[i, j] == 0 && fam[i, j + 1] == 0) {
                pos <- pedpos.insert.right(pos, i, c(j, j + 1))
                j <- j + 2
                next
            }
            # The first of a sibship; gather all sibs and spouses and insert
            # them flanking the parents
            if(fam[i, j] == 0)
                famid <- fam[i, j + 1]
            else
                famid <- fam[i, j]
            sibs <- which(fam[i, ] == famid)
            if(min(sibs) < j)
                stop(sprintf(paste("Unexpected condition:",
                            "starting in the middle of a sibship at (i, j)",
                            "= (%i, %i)"),
                        i, j))
            k <- max(sibs)
            while(spouse[i, k] != 0)
                k <- k + 1

            jmid <- sibs[ceiling((length(sibs) - 1)/2)]
            if(length(jmid) == 0) {
                lxsibs <- integer(0)
                rxsibs <- j:k
            } else {
                while(spouse[i, jmid] != 0)
                    jmid <- jmid + 1
                lxsibs <- j:jmid
                rxsibs <- ((jmid + 1):k)
            }

            dad <- which(nid[i - 1, ] == ped$findex[nid[i, min(sibs)]])
            mom <- which(nid[i - 1, ] == ped$mindex[nid[i, min(sibs)]])
            if(length(dad) != 1 || length(mom) != 1)
                stop(sprintf(paste("Error determining parents:",
                            "(i, j, k, sibs_min, sibs_max, n_dad, n_mom)",
                            "= (%i, %i, %i, %i, %i, %i, %i)"),
                        i, j, k, min(sibs), max(sibs), length(dad),
                        length(mom)))
            if(abs(mom - dad) != 1)
                stop(sprintf(paste("Error determining parents:",
                            "(i, j, k, sibs_min, sibs_max, dad, mom)",
                            "= (%i, %i, %i, %i, %i, %i, %i)"),
                        i, j, k, min(sibs), max(sibs), dad, mom))

            pos <- pedpos.insert.flanking(pos, i, lxsibs, rxsibs, dad, mom)
            j <- k + 1
            next
        }
    }

    if(any(pos[!is.na(pos)] == 0))
        warning("Some persons were not assigned positions, defaulting to zero")
    pos[is.na(pos)] <- 0        # For consistency with kinship2::align.pedigree

    pedalign$pos <- pos
    return(pedalign)
}

#' Insert person(s) at the rightmost position of a given level
#'
#' @param pos
#'     integer matrix, pos element
#' @param i
#'     pedigree level, row of pos
#' @param js
#'     person(s) to insert, col(s) of pos
#' @return
#'     updated pos matrix
pedpos.insert.right <- function(pos, i, js) {
    nj <- length(js)
    m <- max(pos[i, ], na.rm = TRUE)
    pos[!is.na(pos) & pos > m] <- pos[!is.na(pos) & pos > m] + nj
    pos[i, js] <- m + (1:nj)
    return(pos)
}

#' Insert an extended sibship flanking their parents
#'
#' @param pos
#'     integer matrix, pos element
#' @param i
#'     pedigree level, row of pos into which to insert the extended sibship
#' @param lxsibs
#'     subset of extended sibship to be inserted to the left of parents, col(s)
#'     of pos
#' @param rxsibs
#'     subset of extended sibship to be inserted to the right of parents, col(s)
#'     of pos
#' @param dad,mom
#'     parents on level i - 1, col(s) of pos
#' @return
#'     updated pos matrix
pedpos.insert.flanking <- function(pos, i, lxsibs, rxsibs, dad, mom) {
    lparent <- min(dad, mom)
    rparent <- max(dad, mom)
    poslparent <- pos[i - 1, lparent]
    posrparent <- pos[i - 1, rparent]
    ll <- length(lxsibs)
    rr <- length(rxsibs)
    pos[!is.na(pos) & pos > posrparent] <- pos[!is.na(pos) & pos > posrparent] +
        ll + rr
    pos[i, lxsibs] <- poslparent:(poslparent + ll - 1)
    pos[i - 1, lparent] <- poslparent + ll
    pos[i - 1, rparent] <- poslparent + ll + 1
    pos[i, rxsibs] <- (poslparent + ll + 2):(poslparent + ll + rr + 1)
    return(pos)
}

#' Pedigree permutation in (horizontal) plotting order
#'
#' This function returns a permutation that orders the members of a pedigree
#' according to their horizontal plot position. In other words,
#' \code{ped$id[pedOrder(x)]} gives the pedigree member identifiers in the same
#' order as their horizontal position in pedigree alignment \code{x}.
#'
#' @param pedalign
#'      pedigree alignment, such as returned by alignped
#' @return
#'      integer vector that is a permutation of 1:n, where n is the number of
#'      members in the pedigree
#' @export
ped.order <- function(pedalign) {
    ord <- pedalign$nid[order(pedalign$pos)]
    ord <- ord[ord > 0]
    pedSize <- sum(pedalign$n)
    if(length(ord) != pedSize || !setequal(ord, 1:pedSize))
        stop("Internal error: pedOrder result is not a permutation")
    return(ord)
}

#' Horizontal plot position of pedigree members
#'
#' This is a convenience function that returns a named numeric vector, where
#' the names are the pedigree member identifiers and the values are the
#' horizontal plotting coordinates.
#'
#' @inheritParams alignped
#' @inheritParams ped.order
#' @return
#'      named numeric vector, giving the horizontal plot coordinates and named
#'      with pedigree identifier
#' @export
ped.hpos <- function(ped, pedalign) {
    i <- pedalign$nid > 0
    hpos <- pedalign$pos[i]
    names(hpos) <- as.character(ped$id[pedalign$nid[i]])
    return(hpos)
}
