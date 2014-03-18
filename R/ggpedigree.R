#' Generate A ggplot Object To Draw The Given Pedigree
#'
#' This function is used to render a pedigree using ggplot graphics (based on
#' grid graphics), unlike its cognate plot.pedigree from kinship2 which uses
#' base graphics. The results differ primarily in how symbols are annotated to
#' indicate trait values, sample availability and/or other variables. First,
#' this function supports quantitative traits in addition to qualitative ones.
#' Second, when multiple variables are shown, instead of subdividing each
#' symbol into multiple sectors as in the case of plot.pedigree, ggpedigree
#' uses symbol attributes (aka ggplot aesthetics) -- namely fill color, border
#' color, size and transparency. This can make the plot easier to read.
#'
#' @param ped
#'      Pedigree object, as returned by kinship2::pedigree
#' @param pedalign
#'      (Optional) Pedigree alignment, as returned by align.pedigree.1d.
#'      Defaults to value returned by kinship2::align.pedigree.
#' @param symbolfill
#' @param symbolborder
#' @param symbolsize
#' @param symbolalpha
#'      (Optional) Factors or numeric vectors specifying the values of
#'      qualitative or quantitative traits (or other variables) respectively.
#'      Must correspond to ped$id in length and order of values. Used to
#'      annotate the corresponding pedigree symbols using the fill color,
#'      border color, size and transparency attributes. Logical, integer or
#'      character vectors can also be used, and will be coerced into factors
#'      (not recommended, see Details). If symbolfill is not provided and
#'      ped$affected is available, it is used; use NULL to disable this
#'      behavior.
#' @param fillvalues
#' @param bordervalues
#' @param sizevalues
#' @param alphavalues
#'      (Optional) Named vectors used to map trait or variable values to symbol
#'      attribute values. For factor variables, provide attribute values for
#'      the corresponding factor levels. Alternatively, for factor and numeric
#'      variables, specify the high and low attribute values. Default values
#'      are provided, except for factors with more than 3 levels.
#' @param labelsize
#'      (Optional) Label size multiplier, default 1.  Use 0 to avoid labeling symbols.
#' @return
#'      ggplot object, including legend.
#' @section Details:
#'      Following convention, symbol shape is used to indicate sex. The other
#'      symbol attributes - fill color, border color, size and transparency are
#'      available to the caller. Note that fill color is typically used to
#'      indicate affection status (for the primary trait, if there are
#'      multiple). For convenience, if symbolfill is not specified, it is
#'      automatically set to ped$affected if available. If the latter is a
#'      matrix, the first column is used. The values 0 and 1 are mapped to
#'      unaffected and affected, per the kinship2 documentation.
#'
#'      To maintain readability,
#'      * minimize the number of distinct attribute values in each variable
#'      * avoid using a hue-based scale for more than one attribute
#'      * use a large size range when specifying sizevalues
#'      * avoid using transparency to indicate more than two values
#'
#'      For qualitative/categorical variables as symbol{fill,border,size,alpha}
#'      arguments, use factors instead of logical, integer and character
#'      vectors for two reasons. First, you have control over the order of the
#'      levels which determines how levels are mapped to attribute values. And
#'      second, the levels are listed in the legend and descriptive string
#'      levels are easy to interpret compared to integer codes such as 0 and 1.
#'
#'      Use NA values in the symbol{fill,border,size,alpha} arguments at your
#'      own risk; they are not supported. If needed, they can be assigned an
#'      attribute value by adding a ggplot scale to the return value of
#'      ggpedigree. For instance, \code{ggpedigree(...) +
#'      scale_fill_gradient(..., na.value = "gray50")}.
#'
#'      To fine tune the result, add the appropriate ggplot theme or legend
#'      option to the return value.
#' @examples
#'      data(kinshipExtra)
#'      ggpedigree(ped1)
#'      ggpedigree(ped1, align.pedigree.1d(ped1))
#'      ggpedigree(ped1, symbolfill = ped1$affected)
#'      ggpedigree(ped1,
#'          symbolfill = factor(c("unaffected", "affected")[ped1$affected + 1],
#'              levels = c("unaffected", "affected"))
#'      ggpedigree(ped1, symbolfill = ped1$affected) +
#'          scale_fill_manual(values = c("white", "blue"))
#' @import kinship2 ggplot2 scales
#' @export
ggpedigree <- function(ped,
        pedalign     = align.pedigree(ped),
        symbolfill,
        symbolborder = NULL,
        symbolsize   = NULL,
        symbolalpha  = NULL,
        fillvalues   = NULL,
        bordervalues = NULL,
        sizevalues   = NULL,
        alphavalues  = NULL,
        labelsize    = 1) {

    # CHECK AND CONVERT ARGUMENTS #############################################

    # Pedigree
    id       <- ped$id
    sex      <- ped$sex
    status   <- ped$status
    affected <- ped$affected

    # Pedigree alignment should match pedigree
    nid    <- pedalign$nid
    if(!setequal(seq_along(id), nid[nid != 0]))
        stop("'pedalign' does not match 'ped'")
    n      <- pedalign$n
    fam    <- pedalign$fam
    spouse <- pedalign$spouse
    pos    <- pedalign$pos

    # Symbol attributes
    symbolshape <- sex
    # If appropriate, assign a default value to symbolfill using ped$affected
    if(missing(symbolfill)) {
        if(is.null(affected))
            symbolfill <- NULL
        else {
            if(is.matrix(affected))
                symbolfill <- affected[, 1]
            else
                symbolfill <- affected
            traitValues <- c("unaffected", "affected")
            symbolfill <- factor(traitValues[symbolfill + 1],
                    levels = traitValues)
        }
    }
    # Gather symbol attributes
    symbolAttrs <- list(
        "shape"  = symbolshape,
        "fill"   = symbolfill,
        "border" = symbolborder,
        "size"   = symbolsize,
        "alpha"  = symbolalpha)
    # Hide legend key for unused (NULL-valued) symbol attributes
    showAttrKeys <- lapply(symbolAttrs, function(x)
        ifelse(is.null(x), FALSE, TRUE))
    # Assign "mock" default values to unused (NULL-valued) symbol attributes
    symbolAttrs <- lapply(names(symbolAttrs), function(x) {
        ifelse(is.null(symbolAttrs[[x]]),
            defaultSymbolAttr(name = x, n = length(id)),
            symbolAttrs[[x]])
        })
    # Validate symbol{fill,border,size,alpha}
    symbolAttrs <- lapply(names(symbolAttrs), function(x) {
        validSymbolAttr(symbolAttr = symbolAttrs[[x]], n = length(id), name = x)
        })

    # Attribute value arguments
    shapevalues <- as.integer(c("male" = 1, "female" = 2, "unknown" = 3,
            "terminated" = 4))
    attrValues <- list(
        "shape"  = shapevalues,
        "fill"   = fillvalues,
        "border" = bordervalues,
        "size"   = sizevalues,
        "alpha"  = alphavalues)

    # Transform symbol attributes
    symbolAttrs <- within(symbolAttrs, {
        shape <- attrValues$shape[as.character(shape)]
        fill  <- fill
        size  <- NULL
    })

    if(!is.numeric(labelsize) || length(labelsize) != 1 || !is.finite(labelsize)
        || labelsize < 0)
        stop("'labelsize' must be a single positive finite numeric value")

    # GENERATE PLOTTING DATA ##################################################

    # Generate plotting data from ped and pedalign
    # Values in pos give horizontal position, row number of pos gives vertical
    xpos <- pos
    ypos <- matrix(rep(1:nrow(nid), ncol(nid)), nrow = nrow(nid), byrow = FALSE)
    # xpos and ypos are matrices corresponding to the 2d layout of the pedigree
    # Extract a vector; subset and order it so xpos, ypos correspond to id
    posSubset <- as.vector(nid) != 0
    posOrder <- order(nid[posSubset])
    xpos <- pos[posSubset][posOrder]
    ypos <- ypos[posSubset][posOrder]
    # Sanity check that data elements have same length
    if(length(id) != length(xpos) || length(id) != length(ypos))
        stop("Internal error: xpos and/or ypos do not match id in length")
    # Build ggplot data
    symbolData <- do.call(data.frame, c(list(id = id, x = xpos, y = ypos,
                shape = sex), symbolAttrs))
    symbolAttrMapping <- aes_all(c("x", "y", "shape", names(symbolAttrs)))

    # INITIALIZE PLOT #########################################################
    plt <- ggplot(mapping = aes(x = x, y = y)) +
        theme(panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())

    # Build ggplot symbol layer
    plt <- plt +
        geom_point(
            data = symbolData,
            mapping = aes(
                shape = as.integer(c(15, 16, 18, 17)[shape]))) +
        geom_point(
            data = symbolData,
            mapping = aes(
                shape = as.integer(c(15, 16, 18, 17)[shape]))) +
        geom_point(
            data = symbolData,
            mapping = aes(
                shape = as.integer(c(22, 21, 23, 24)[shape]))) +
        scale_x_continuous(expand = c(0.1, 0.1)) +
        scale_y_continuous(expand = c(0.1, 0.1)) +
        scale_shape_identity() +
        guides(shape = "none", fill = "none", color = "none", size = "none",
            alpha = "none")

    if("size" %in% symbolAttrMapping)
        plt <- plt +
            geom_point(data = symbolData, mapping = symbolAttrMapping) +
            scale_size_manual(values = c(7.5, 12.5))
    else
        plt <- plt +
            geom_point(data = symbolData, mapping = symbolAttrMapping,
                size = 7.5)
    plt <- plt +
        scale_y_continuous(trans = reverse_trans(), expand = c(0.1, 0)) +
        scale_shape_manual(values = c("male" = 22, "female" = 21,
                "unknown" = 23, "terminated" = 24)) +
        scale_fill_manual(values = c("white", "gray25")) +
        scale_color_manual(values = c("black", "red")) +
        scale_alpha_manual(values = c(1, 0.5)) +
        guides(shape = guide_legend(title = NULL),
            fill  = guide_legend(title = NULL, override.aes = list(shape = 22)),
            color = guide_legend(title = NULL, override.aes = list(shape = 22)),
            size  = guide_legend(title = NULL, override.aes = list(shape = 22)),
            alpha = guide_legend(title = NULL, override.aes = list(shape = 22)))
    # if("color" %in% symbolAttrMapping) {
    #     if("size" %in% symbolAttrMapping)
    #         plt <- plt +
    #             geom_point(data = symbolData,
    #                 mapping = aes_all(c("x", "y", "shape", setdiff(names(symbolAttrs), "fill")))) +
    #             scale_size_manual(values = c(7.5, 12.5) * 1.25)
    #     else
    #         plt <- plt +
    #             geom_point(data = symbolData,
    #                 mapping = aes_all(c("x", "y", "shape", setdiff(names(symbolAttrs), "fill"))),
    #                 size = 7.5 * 1.25)
    # }

    # PLOT LINE SEGMENTS SHOWING STATUS (ALIVE or DECEASED)

    # PLOT LINE SEGMENTS SHOWING RELATIONSHIPS

    # PLOT TEXT LABELS UNDER SYMBOLS

    return(plt)
}

#' Default symbol attributes
#'
#' @param name
#'      String. One of: shape, fill, border, size, and alpha
#' @param n
#'      Integer. Length of the attribute vector
#' @return
#'      Factor or numeric vector. Default value of symbol attributes suitable
#'      for symbol{shape,fill,border,size,alpha}
defaultSymbolAttr <- function(name, n) {
    defaultAttr <- switch(name,
        "shape"  = factor(rep("male", n)),
        "fill"   = factor(rep("unaffected", n)),
        "border" = factor(rep("default", n)),
        "size"   = as.numeric(rep(1, n)),
        "alpha"  = factor(rep("default", n)),
        stop(sprintf(paste("Internal error: '%s' is an invalid 'name'",
                    "argument; must be one of:",
                    "shape, fill, border, size, and alpha"),
                name)))
    return(defaultAttr)
}

#' Check validity of symbol attribute
#'
#' @param symbolAttr
#'      Symbol attribute to be checked.
#' @param n
#'      Integer. Required length of the attribute vector
#' @param name
#'      String. Name of the attribute being checked. One of: shape, fill,
#'      border, size, and alpha
#' @return
#'      Factor or numeric vector of length n
validSymbolAttr <- function(symbolAttr, n, name) {
    symbolAttr <- switch(class(symbolAttr),
        "factor" = symbolAttr,
        "numeric" = symbolAttr,
        "character" = factor(symbolAttr),
        "logical" = factor(symbolAttr),
        "integer" = factor(symbolAttr),
        stop(sprintf(paste("symbol%s must be a vector of type factor or,",
                    "numeric (or character, logical or integer)"),
                name)))
    if(length(symbolAttr) != n)
        stop(sprintf("symbol%s has length that does not match 'ped'",
                name))
    # if(any(is.na(symbolAttr)))
    #     stop(sprintf("%s: %s has NA(s), not allowed",
    #             msgPrefix, symbolAttr))
    return(symbolAttr)
}

#' Default attribute values
#'
#' @param name
#'      String. One of: fill, border, size, and alpha
#' @param symbolAttr
#'      Factor or numeric vector. Symbol attributes.
#' @return
#'      Character vector. Default value of attribute values suitable for
#'      {fill,border,size,alpha}values
#' @import RColorBrewer
defaultAttrValue <- function(name, symbolAttr) {
    if(is.factor(symbolAttr)) {
        if(nlevels(symbolAttr) > 3)
            stop(sprintf(paste("'symbol%s' argument has more than 3 levels.",
                        "Specify '%svalues'"), name, name))
        attrValues <- switch(name,
            "fill" = c("white", "gray25", "gray65"),
            "color" = brewer.pal(3, "Dark2"),
            "size" = c(3, 15),
            "alpha" = c(1, 0.2),
            stop(sprintf(paste("Internal error: '%s' is an invalid 'name'",
                        "argument; must be one of:",
                        "shape, fill, border, size, and alpha"),
                    name)))
    } else if(is.numeric(symbolAttr)) {
        attrValues <- switch(name,
            "fill" = c("white", "gray25", "gray65"),
            "color" = brewer.pal(3, "Dark2"),
            "size" = c(3, 15),
            "alpha" = c(1, 0.2),
            stop(sprintf(paste("Internal error: '%s' is an invalid 'name'",
                        "argument; must be one of:",
                        "shape, fill, border, size, and alpha"),
                    name)))
    } else {
        stop(sprintf(paste("Internal error: '%s' is an invalid 'name'",
                    "argument; must be one of:",
                    "shape, fill, border, size, and alpha"),
                name)))
    }
        # bordervalues = c("black", "red", "cyan"),
        # sizevalues   = c(7.5, 10, 5),
        # alphavalues  = c(1, 0.3, 0.6),
}
