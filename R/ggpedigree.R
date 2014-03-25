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
#'      the corresponding factor levels. For numeric variables, specify the
#'      high and low attribute values. Default values are provided, except for
#'      factors with more than 3 levels.
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

    # CHECK ARGUMENTS #########################################################

    # Pedigree alignment should match pedigree
    if(!setequal(seq_along(ped$id), pedalign$nid[pedalign$nid != 0]))
        stop("'pedalign' does not match 'ped'")

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
    symbolAttrs <- list(
        "shape"  = ped$sex,
        "fill"   = symbolfill,
        "border" = symbolborder,
        "size"   = symbolsize,
        "alpha"  = symbolalpha)
    # Hide legend key for unused (NULL-valued) symbol attributes
    symbolAttrWasNull <- lapply(symbolAttrs, is.null)
    symbolAttrs <- validSymbolAttrs(symbolAttrs, length(ped$id))

    attrValues <- list(
        "shape"  = c("male" = 1, "female" = 2, "unknown" = 3, "terminated" = 4),
        "fill"   = fillvalues,
        "border" = bordervalues,
        "size"   = sizevalues,
        "alpha"  = alphavalues)
    attrValueWasNull <- lapply(attrValues, is.null)
    attrValues <- validAttrValues(attrValues, symbolAttrs)

    if(!is.numeric(labelsize) || length(labelsize) != 1 || !is.finite(labelsize)
        || labelsize < 0)
        stop("'labelsize' must be a single positive finite numeric value")

    # GENERATE PLOTTING DATA ##################################################

    # Special handling of shape and size is needed due to the non-uniform
    # apparent size of the standard plotting symbol shapes
    symbolAttrs <- within(symbolAttrs, {
        shape <- attrValues$shape[as.character(shape)]
        size  <- ifelse(is.numeric(size),
                scaleMinMax(size,
                    attrValues$size["low"],
                    attrValues$size["high"]),
                attrValues$size[as.character(size)])
    })
    symbolData   <- getSymbolData(ped, pedalign, symbolAttrs)
    statusData   <- getStatusData(symbolData, attrValues$size)
    relationData <- getRelationData(ped, pedalign)

    # INITIALIZE PLOT #########################################################
    plt <- ggplot() +
        theme(
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text  = element_blank(),
            axis.ticks = element_blank(),
            axis.line  = element_blank())

    # PLOT LINE SEGMENTS SHOWING RELATIONSHIPS

    # PLOT SYMBOLS SHOWING INDIVIDUALS
    plt <- plt +
        geom_point(                # cover up relationship lines under symbols
            data    = symbolData,
            mapping = aes(
                x     = x,
                y     = y,
                shape = as.integer(c(15, 16, 18, 17)[shape]),
                size  = as.numeric(size * c(1.25, 1.45, 1.51, 1.05)[shape])),
            color   = "gray90",
            alpha   = 1) +
        geom_point(                # symbol border
            data    = symbolData,
            mapping = aes(
                x     = x,
                y     = y,
                shape = as.integer(c(15, 16, 18, 17)[shape]),
                color = border,
                size  = as.numeric(size * c(1.25, 1.45, 1.51, 1.05)[shape]),
                alpha = alpha)) +
        geom_point(                # spacer
            data    = symbolData,
            mapping = aes(
                x     = x,
                y     = y,
                shape = as.integer(c(15, 16, 18, 17)[shape]),
                size  = as.numeric(size * c(1.05, 1.24, 1.26, 0.82)[shape]),
                alpha = alpha),
            color   = "gray90")) +
        geom_point(                # symbol fill
            data    = symbolData,
            mapping = aes(
                x     = x,
                y     = y,
                shape = as.integer(c(22, 21, 23, 24)[shape]),
                fill  = fill,
                size  = as.numeric(size * c(1, 1.05, 0.8, 0.6)[shape]),
                alpha = alpha),
            color = "gray20") +
        scale_x_continuous(expand = c(0.1, 0.1)) +
        scale_y_continuous(expand = c(0.1, 0.1), trans = reverse_trans()) +
        scale_shape_identity() +
        ifelse(is.factor(symbolAttrs$fill),
            scale_fill_manual(values = attrValues$fill),
            scale_fill_gradient(
                low  = attrValues$fill["low"],
                high = attrValues$fill["high"])) +
        ifelse(is.factor(symbolAttrs$border),
            scale_color_manual(values = attrValues$border),
            scale_color_gradient(
                low  = attrValues$border["low"],
                high = attrValues$border["high"])) +
        scale_size_identity() +
        ifelse(is.factor(symbolAttrs$alpha),
            scale_alpha_manual(values = attrValues$alpha),
            scale_alpha_continuous(range = c(
                    low  = attrValues$alpha["low"],
                    high = attrValues$alpha["high"]))) +
        guides(
            shape = "none",
            fill  = "none",
            color = "none",
            size  = "none",
            alpha = "none")

    # PLOT LINE SEGMENTS SHOWING STATUS (ALIVE or DECEASED)

    # PLOT TEXT LABELS UNDER SYMBOLS

    # PLOT LEGEND
        # guides(shape = guide_legend(title = NULL),
        #     fill  = guide_legend(title = NULL, override.aes = list(shape = 22)),
        #     color = guide_legend(title = NULL, override.aes = list(shape = 22)),
        #     size  = guide_legend(title = NULL, override.aes = list(shape = 22)),
        #     alpha = guide_legend(title = NULL, override.aes = list(shape = 22)))

    return(plt)
}

#' Return valid symbol attributes
#'
#' @param symbolAttrs
#'      Named list of symbol attributes, with elements shape, fill, border,
#'      size, and alpha
#' @param n
#'      Integer giving the required length of each symbol attribute
#' @return
#'      Named list of symbol attributes, with default attributes in place of
#'      NULLs and checking all attributes for validity
validSymbolAttrs <- function(symbolAttrs, n) {
    # Assign mock default values to unused (NULL-valued) symbol attributes
    symbolAttrs <- lapply(names(symbolAttrs), function(name) {
        ifelse(is.null(symbolAttrs[[name]]),
            defaultSymbolAttr(name, n),
            symbolAttrs[[name]])
        })
    # Check that attributes are valid
    symbolAttrs <- lapply(names(symbolAttrs), function(name) {
        checkSymbolAttr(symbolAttrs[[name]], n, name)
        })
    return(symbolAttrs)
}

#' Default symbol attributes
#'
#' @param name
#'      String. One of: shape, fill, border, size, and alpha
#' @inheritParams validSymbolAttrs
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
#' @inheritParams validSymbolAttrs
#' @inheritParams defaultSymbolAttr
#' @return
#'      Factor or numeric vector of length n
checkSymbolAttr <- function(symbolAttr, n, name) {
    symbolAttr <- switch(class(symbolAttr),
        "factor"    = symbolAttr,
        "numeric"   = symbolAttr,
        "character" = factor(symbolAttr),
        "logical"   = factor(symbolAttr),
        "integer"   = factor(symbolAttr),
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

#' Return valid attribute values
#'
#' @param attrValues
#'      Named list of attribute values, with elements shape, fill, border,
#'      size and alpha
#' @inheritParams validSymbolAttrs
#' @return
#'      Named list of valid attribute values, with NULL values replaced with
#'      suitable defaults
validAttrValues <- function(attrValues, symbolAttrs) {
    attrValues <- lapply(names(attrValues), function(name) {
        ifelse(is.null(attrValues[[name]]),
            defaultAttrValue(name, symbolAttrs[[name]]),
            attrValues[[name]])
        })
    attrValues <- lapply(names(attrValues), function(name) {
        checkAttrValue(attrValues[[name]], symbolAttrs[[name]], name)
    })
    return(attrValues)
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
#' @import RColorBrewer, scales
defaultAttrValue <- function(name, symbolAttr) {
    if(is.factor(symbolAttr)) {
        if(nlevels(symbolAttr) > 3)
            stop(sprintf(paste("'symbol%s' argument has more than 3 levels.",
                        "Specify '%svalues'"), name, name))
        attrValue <- switch(name,
            "fill"   = switch(nlevels(symbolAttr),
                        c("white"),
                        c("white", "gray25"),
                        c("white", "gray65", "gray25")),
            "border" = switch(nlevels(symbolAttr),
                        c("gray90"),
                        c("gray90", "#C66E00"),
                        c("gray90", "#009DCF", "#C66E00")),
            "size"   = c(1, 2, 3)[1:nlevels(symbolAttr)],
            "alpha"  = switch(nlevels(symbolAttr),
                        c(1),
                        c(1, 0.2),
                        c(1, 0.6, 0.2)),
            stop(sprintf(paste("Internal error: '%s' is an invalid 'name'",
                        "argument; must be one of:",
                        "shape, fill, border, size, and alpha"),
                    name)))
        names(attrValue) <- levels(symbolAttr)
    } else if(is.numeric(symbolAttr)) {
        attrValue <- switch(name,
            "fill"   = c(low = "white",   high = "black"),
            "border" = c(low = "#fee5d9", high = "#a50f15"),
            "size"   = c(low = 3,         high = 15),
            "alpha"  = c(low = 0.2,       high = 1),
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
    return(attrValue)
}

#' Check validity of attribute values
#'
#' @param attrValue
#'      Named vector used to map symbol attributes to plot values
#' @inheritParams defaultAttrValue
#' @return
#'      Named vector of valid attribute values
checkAttrValue <- function(attrValue, symbolAttr, name) {
    if(is.factor(symbolAttr)) {
        if(!setequal(levels(symbolAttr), names(attrValue))) {
            if(all(levels(symbolAttr) %in% names(attrValue))) {
                attrValue <- attrValue[levels(symbolAttr)]
            } else if(is.null(names(attrValue)) &&
                length(attrValue) == nlevels(symbolAttr)) {
                names(attrValue) <- levels(symbolAttr)
            } else {
                stop(sprintf(paste("%svalues must be a named vector,",
                        "with values for each level of the factor symbol%s"),
                    name, name))
            }
        }
    } else if (is.numeric(symbolAttr)) {
        if(!setequal(c("low", "high"), names(attrValue))) {
            if(all(c("low", "high") %in% names(attrValue))) {
                attrValue <- attrValue[c("low", "high")]
            } else if(is.null(names(attrValue)) && length(attrValue) == 2) {
                names(attrValue) <- c("low", "high")
            } else {
                stop(sprintf(paste("%svalues must be a named vector with",
                            "values named 'low' and 'high' to scale symbol%s"),
                        name, name))
            }
        }
    } else {
        stop(sprintf(paste("Internal error: symbol%s must be a factor or",
                    "numeric vector"), name))
    }
    # fill and border values take multiple valid formats e.g. "blue", "#aabbcc"
    # size and alpha values must be numeric
    if(name %in% c("size", "alpha")) {
        if(!all(is.numeric(attrValue))) {
            stop(sprintf("%svalues must be numeric", name))
        }
    }
    return(attrValue)
}

#' Scale linearly between two bounds
#'
#' @param x
#'      Numeric values to be scaled
#' @param newMin
#' @param newMax
#'      Lower and upper bounds of the scaled values
#' @return
#'      Numeric values linearly scaled between min and max
scaleMinMax <- function(x, newMin, newMax) {
    return( newMin + (newMax - newMin) * (x - min(x)) / (max(x) - min(x)) )
}

#' Construct symbol plotting data
#'
#' @inheritParams ggpedigree
#' @inheritParams validSymbolAttrs
#' @return
#'      dataframe to be use as ggplot data for symbol layers
getSymbolData <- function(ped, pedalign, symbolAttrs) {
    # Values in pos give horizontal position, row number of pos gives vertical
    xpos <- pedalign$pos
    ypos <- matrix(rep(1:nrow(pedalign$nid), ncol(pedalign$nid)),
            nrow = nrow(pedalign$nid), byrow = FALSE)
    # xpos and ypos are matrices corresponding to the 2d layout of the pedigree
    # Extract a vector; subset and order it so xpos, ypos correspond to ped$id
    posSubset <- as.vector(pedalign$nid) != 0
    posOrder <- order(pedalign$nid[posSubset])
    xpos <- xpos[posSubset][posOrder]
    ypos <- ypos[posSubset][posOrder]
    # Sanity check that data elements have same length
    if(length(ped$id) != length(xpos) || length(ped$id) != length(ypos))
        stop("Internal error: xpos and/or ypos do not match ped$id in length")
    symbolData <- do.call(data.frame, c(list(id = ped$id, x = xpos, y = ypos),
            symbolAttrs))
    return(symbolData)
}

#' Construct symbol status plotting data to indicate deceased individuals
#'
#' @param symbolData
#'      dataframe with symbol plotting data
#' @param sizeValues
#'      named vector mapping symbolData$size
getStatusData <- function(symbolData, attrValues$size) {
    dead <- ped$id[ped$status != 0]
    symbolData <- symbolData[symbolData$id %in% dead, c("id", "x", "y", "size")]
    statusData
}
