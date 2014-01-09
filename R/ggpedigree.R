#' Generate a ggplot object to draw the given pedigree
#'
#' @param ped
#'     pedigree object, as returned by kinship2::pedigree
#' @param pedalign
#'     pedigree alignment, as returned by align.pedigree.1d. Optional, defaults
#'     to value returned by kinship2::align.pedigree.
#' @param symbolfill
#' @param symbolcolor
#' @param symbolsize
#' @param symbolalpha
#'     (optional) vectors of same length as ped$id, indicating the fill, color,
#'     size and transparency of the corresponding symbols. Vectors can be of
#'     any type that can be coerced to a factor i.e. logical, integer, factor
#'     or character. Each vector must have two unique values (e.g. "affected"
#'     and "unaffected"); NA is not allowed. If a factor is provided, the first
#'     and second levels are mapped to the default and alternative symbol
#'     renderings i.e. to empty/grayscale/regular-sized/opaque symbols and to
#'     solid/color/large/transparent symbols.
#' @param labelsize
#'     size multiplier, default 1
#' @return
#'     ggplot object with 3 layers: symbols, segments and labels
#' @import kinship2 ggplot2
#' @export
ggpedigree <- function(ped, pedalign, symbolfill, symbolcolor, symbolalpha,
    symbolsize = 1, labelsize = 1) {

    # Unpack and check arguments
    if(missing(pedalign))
        pedalign <- align.pedigree(ped)
    id <- ped$id
    sex <- ped$sex
    affected <- ped$affected
    status <- ped$status
    nid <- pedalign$nid
    pos <- pedalign$pos
    if(is.matrix(affected))
        if(ncol(affected) > 1)
            stop(paste("Multiple traits i.e. pedigree with multi-column",
                    "'affected' element is not supported"))
        else
            affected <- as.vector(affected)

    # Plot symbols and legends

    # Generate plotting data from ped and pedalign
    # Values in pos give horizontal position, row number of pos gives vertical
    xpos <- pos
    ypos <- matrix(rep(1:nrow(nid), ncol(nid)), nrow = nrow(nid), byrow = FALSE)
    # xpos and ypos are matrices corresponding to the 2d layout of the pedigree
    # Extract vector, subset and order so xpos, ypos correspond to id
    posSubset <- as.vector(nid) != 0
    posOrder <- order(nid[posSubset])
    xpos <- pos[posSubset][posOrder]
    ypos <- ypos[posSubset][posOrder]
    # Sanity check that data elements have same length
    if(length(id) != length(xpos) || length(id) != length(ypos))
        stop("Internal error: xpos and/or ypos do not match id in length")
    # Build ggplot data
    symbolData <- data.frame(id = id, x = xpos, y = ypos, sex = sex,
        affected = affected, status

    # Plot lines showing relationships

}
