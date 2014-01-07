#' Generate a "1-dimensional" layout of a pedigree
#'
#' Given a pedigree, this function returns a pedigree alignment that assigns
#' each person a unique horizontal position, making it suitable for plotting
#' persons in the pedigree along an axis, while depicting their relationships.
#' This is useful in annotating a figure (e.g. haplotype sharing) where the
#' same persons are plotted along one of the axes.
#'
#' @param ped
#'     pedigree object, as returned by kinship2::pedigree
#' @param ...
#'     additional arguments passed to kinship2::align.pedigree
#' @return
#'     list, as returned by kinship2::align.pedigree
#' @import kinship2
#' @export
align.pedigree.1d <- function(ped,
        method = c("compactsubtree", "preserveorder"),
        ...) {

    align.result <- align.pedigree(ped, ...)
    pedheight <- length(align.result$n)
    if(pedheight < 2)
        return(align.result)

    n <- align.result$n
    nid <- align.result$nid
    fam <- align.result$fam
    spouse <- align.result$spouse
    pos <- align.result$pos

    pos[nid == 0] <- NA
    method <- match.arg(method)
    pos <- switch(method,
        "compactsubtree" = align.pedigree.1d.compactsubtree(pos, n, nid, fam,
            spouse, ped),
        "preserveorder" = align.pedigree.1d.preserveorder(pos),
        stop(sprintf("%s is not a valid method argument to %s", method,
                match.call()[[1]])))
    if(any(pos[!is.na(pos)] == 0))
        warning("Some persons were not assigned positions, defaulting to zero")
    pos[is.na(pos)] <- 0

    align.result$pos <- pos
    return(align.result)
}

#' Generate a 1d layout of a pedigree, preserving the original horizontal order
#'
#' @param pos
#'     integer matrix, pos element from return value of kinship2::align.pedigree
#' @return
#'     integer matrix; updated pos element
align.pedigree.1d.preserveorder <- function(pos) {
    pos[!is.na(pos)] <- rank(pos[!is.na(pos)])
    return(pos)
}

#' Generate a 1d layout of a pedigree where each subtree is compact
#'
#' A subtree of a pedigree comprises a person, their partner, descendants and
#' descendants' partners. This layout generates compact subtrees, in that each
#' one occupies a horizontal interval _exclusively_. In other words, no persons
#' outside a subtree (e.g. ancestors, siblings, etc.) occupy the interval
#' corresponding to it.
#'
#' This layout is the most "natural" looking 1d layout.
#'
#' @param pos
#' @param n
#' @param nid
#' @param fam
#' @param spouse
#'     elements of the return value of kinship2::align.pedigree
#' @param ped
#'     pedigree object, as returned by kinship2::pedigree
#' @return
#'     integer matrix; updated pos element
align.pedigree.1d.compactsubtree <- function(pos, n, nid, fam, spouse, ped) {

    # NA indicates empty position; all non-empty positions are zero-ed out
    pos[!is.na(pos)] <- 0
    # Process levels one at a time, top-down
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
            while(spouse[i, k] != 0 && fam[i, k + 1] == 0)
                k <- k + 1

            jmid <- sibs[ceiling((length(sibs) - 1)/2)]
            while(spouse[i, jmid] != 0 && fam[i, jmid + 1] == 0)
                jmid <- jmid + 1
            lxsibs <- j:jmid
            if(jmid < k)
                rxsibs <- ((jmid + 1):k)
            else
                rxsibs <- integer(0)

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
    return(pos)
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
#' @param dad
#' @param mom
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
