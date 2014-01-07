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
#' @return
#'     list, as returned by kinship2::align.pedigree
#' @import kinship2
#' @export
align.pedigree.1d <- function(ped, ...) {
    align.result <- align.pedigree(ped, ...)

    nlev <- length(align.result$n)
    if(nlev < 2) return(align.result)

    nid <- align.result$nid
    pos <- align.result$pos
    pos[nid == 0] <- NA
    fam <- align.result$fam
    spouse <- align.result$spouse

    # Crude method, nudge to grid preserving horizontal order across levels
    # pos[!is.na(pos)] <- rank(pos[!is.na(pos)])

    # Subtree local method, each subtree (individual, partner and all
    # descendants) occupies an interval exclusively
    #
    # NA indicates empty position; all non-empty positions are zero-ed out
    pos[!is.na(pos)] <- 0
    # Process levels one at a time, top-down
    for(i in 1:nlev) {
        # Process each level left-to-right in extended sibships
        # A extended sibship is a set of full siblings with their partners.
        # Non-founder partners are assumed to be depicted a second time on the
        # right with their ancentors and siblings.
        j <- 1
        while(j <= ncol(pos)) {
            # Empty cell, move on
            if(nid[i, j] == 0) {
                j <- j + 1
                next
            }
            # This case shouldn't occur
            # i.e. person is neither a descendant nor spouse
            if(fam[i, j] == 0 && spouse[i, j] == 0) {
                stop(sprintf("Unsupported pedigree: fam = 0 and spouse = 0 for (i, j) = (%i, %i)", i, j))
                next
            }
            # Pair of founders, insert at the right of the current level
            if(fam[i, j] == 0 && fam[i, j + 1] == 0) {
                pos <- insert.right(pos, i, c(j, j + 1))
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
                stop(sprintf("Unexpected condition: starting in the middle of a sibship at (i, j) = (%i, %i)", i, j))
            k <- max(sibs)
            while(spouse[i, k] != 0 && fam[i, k + 1] == 0)
                k <- k + 1
            dad <- which(nid[i - 1, ] == ped$findex[nid[i, min(sibs)]])
            mom <- which(nid[i - 1, ] == ped$mindex[nid[i, min(sibs)]])
            if(length(dad) != 1 || length(mom) != 1)
                stop(sprintf("Error determining parents: (i, j, k, sibs_min, sibs_max, n_dad, n_mom) = (%i, %i, %i, %i, %i, %i, %i)",
                        i, j, k, min(sibs), max(sibs), length(dad),
                        length(mom)))
            if(abs(mom - dad) != 1)
                stop(sprintf("Error determining parents: (i, j, k, sibs_min, sibs_max, dad, mom) = (%i, %i, %i, %i, %i, %i, %i)",
                        i, j, k, min(sibs), max(sibs), dad, mom))
            jmid <- sibs[floor((length(sibs) + 1)/2)]
            while(spouse[i, jmid] != 0 && fam[i, jmid + 1] == 0)
                jmid <- jmid + 1
            lxsibs <- j:jmid
            if(jmid < k)
                rxsibs <- ((jmid + 1):k)
            else
                rxsibs <- integer(0)
            pos <- insert.flanking(pos, i, lxsibs, rxsibs, dad, mom)
            j <- k + 1
            next
        }
    }
    if(any(pos[!is.na(pos)] == 0))
        warning("Some persons were not assigned positions, defaulting to zero")

    pos[is.na(pos)] <- 0
    align.result$pos <- pos

    return(align.result)
}

insert.right <- function(pos, i, js) {
    nj <- length(js)
    m <- max(pos[i, ], na.rm = TRUE)
    pos[!is.na(pos) & pos > m] <- pos[!is.na(pos) & pos > m] + nj
    pos[i, js] <- m + (1:nj)
    return(pos)
}

insert.flanking <- function(pos, i, lxsibs, rxsibs, dad, mom) {
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
