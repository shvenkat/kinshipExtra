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
#'      (Optional) Pedigree alignment list, defaults to alignped(ped)
#' @param symbolfill,symbolborder,symbolsize,symbolalpha
#'      (Optional) Factors or numeric vectors specifying the values of
#'      qualitative or quantitative traits (or other variables) respectively.
#'      Must correspond to ped$id in length and order of values. Used to
#'      annotate the corresponding pedigree symbols using the fill color,
#'      border color, size and transparency attributes. Logical, integer or
#'      character vectors can also be used, and will be coerced into factors
#'      (not recommended, see Details). If symbolfill is not provided and
#'      ped$affected is available, it is used; use NULL to disable this
#'      behavior.
#' @param fillvalues,bordervalues,sizevalues,alphavalues
#'      (Optional) Named vectors used to map trait or variable values to symbol
#'      attribute values. For factor variables, provide attribute values for
#'      the corresponding factor levels. For numeric variables, specify the
#'      high and low attribute values. Default values are provided, except for
#'      factors with more than 3 levels.
#' @param fillkey,borderkey,sizekey,alphakey
#'      (Optional) String value to modify the appearance of a symbol attribute
#'      key. NULL (the default) returns no key, "" returns an untitled key and
#'      any other string is used as the key title.
#' @param labeltext
#'      (Optional) Character vector of labels to be printed below the symbols
#'      in addition to ped$id.
#' @param labelsize
#'      (Optional) Label size multiplier, default 1. Use 0 to avoid labeling
#'      symbols.
#' @param labeloffset
#'      (Optional) Label offset value, default 0. Used to adjust vertical space
#'      between symbols and labels.
#' @param bgcolor
#'      (Optional) Plot background color, default "gray90" i.e. light gray
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
#'      * minimize the number of distinct attribute values in each variable,
#'      * avoid using a hue-based scale for more than one attribute,
#'      * use a large size range when specifying sizevalues, and
#'      * avoid using transparency to indicate more than two values.
#'
#'      For qualitative/categorical variables as symbol.. arguments, use
#'      factors instead of logical, integer and character vectors for two
#'      reasons. First, you have control over the order of the levels which
#'      determines how levels are mapped to attribute values. And second, the
#'      levels are listed in the legend and descriptive string levels are easy
#'      to interpret compared to integer codes such as 0 and 1.
#'
#'      Use NA values in the symbol.. arguments at your own risk; they are not
#'      supported. If needed, they can be assigned an attribute value by adding
#'      a ggplot scale to the return value of ggpedigree. For instance,
#'      \code{ggpedigree(...) + scale_fill_gradient(..., na.value = "gray50")}.
#'
#'      To fine tune the result, try adding the appropriate ggplot theme or
#'      legend option to the return value. Such use has not been tested.
#' @examples
#'      library(kinship2)
#'      library(ggplot2)
#'      data(simpleped)
#'      ggpedigree(ped1)
#'      ggpedigree(ped1, alignped(ped1, method = "1d"))
#'      ggpedigree(ped1, symbolfill = ped1$affected)
#'      ggpedigree(ped1,
#'          symbolfill = factor(c("unaffected", "affected")[ped1$affected + 1],
#'                  levels = c("unaffected", "affected")))
#'      ggpedigree(ped1,
#'          symbolfill = factor(ped1$affected, levels = c("0", "1"))) +
#'          scale_fill_manual(values = c("white", "red"))
#' @import kinship2 ggplot2 scales
#' @export
ggpedigree <- function(ped,
                       pedalign = alignped(ped),
                       symbolfill,        symbolborder = NULL,
                       symbolsize = NULL, symbolalpha  = NULL,
                       fillvalues = NULL, bordervalues = NULL,
                       sizevalues = NULL, alphavalues  = NULL,
                       fillkey    = NULL, borderkey    = NULL,
                       sizekey    = NULL, alphakey     = NULL,
                       labeltext  = NULL, labelsize    = 1,
                       labeloffset = 0,   bgcolor      = "gray90") {

    # CHECK ARGUMENTS -------------------------------------------------------

    # Pedigree alignment should match pedigree
    if(!setequal(seq_along(ped$id), pedalign$nid[pedalign$nid != 0]))
        stop("'pedalign' does not match 'ped'")

    # If appropriate, assign a default value to symbolfill using ped$affected
    if(missing(symbolfill)) {
        if(is.null(ped$affected))
            symbolfill <- NULL
        else {
            if(is.matrix(ped$affected))
                symbolfill <- ped$affected[, 1]
            else
                symbolfill <- ped$affected
            traitValues <- c("unaffected", "affected")
            symbolfill <- factor(traitValues[symbolfill + 1],
                                 levels = traitValues)
        }
    }

    symbolAttrs <- list("shape"  = ped$sex,
                        "fill"   = symbolfill, "border" = symbolborder,
                        "size"   = symbolsize, "alpha"  = symbolalpha)
    symbolAttrWasNull <- lapply(symbolAttrs, is.null)
    symbolAttrs <- validSymbolAttrs(symbolAttrs, length(ped$id))

    attrValues <- list("shape"  = c("male" = 1, "female" = 2,
                                    "unknown" = 3, "terminated" = 4),
                       "fill"   = fillvalues, "border" = bordervalues,
                       "size"   = sizevalues, "alpha"  = alphavalues)
    # attrValueWasNull <- lapply(attrValues, is.null)
    attrValues <- validAttrValues(attrValues, symbolAttrs, bgcolor)

    attrKeys <- list("shape"  = NULL,
                     "fill"   = fillkey,  "border" = borderkey,
                     "size"   = sizekey,  "alpha"  = alphakey)
    attrKeys <- validAttrKeys(attrKeys, symbolAttrWasNull)

    if(!is.null(labeltext)) {
        if(!is.character(labeltext) || length(labeltext) != length(ped$id))
            stop(paste("labeltext must be a character vector of the same",
                       "length as ped$id"))
    }
    if(!is.numeric(labelsize) || length(labelsize) != 1 ||
       !is.finite(labelsize) || labelsize < 0)
        stop("'labelsize' must be a single positive finite numeric value")

    if(!is.numeric(labeloffset) || length(labeloffset) != 1 ||
       !is.finite(labeloffset))
        stop("'labeloffset' must be a single finite numeric value")

    # GENERATE PLOTTING DATA ------------------------------------------------

    # Prepare to tweak symbol size when building the plot so that different
    # shapes appear similar in size
    symbolAttrs <- within(symbolAttrs, {
        shape <- attrValues$shape[as.character(shape)]
        size  <- switch(ifelse(is.numeric(size), 1, 2),
                        scaleMinMax(size,
                                    attrValues$size["low"],
                                    attrValues$size["high"]),
                        attrValues$size[as.character(size)])
    })
    relationData <- getRelationData(ped, pedalign)
    symbolData   <- getSymbolData(ped, pedalign, symbolAttrs)
    statusData   <- getStatusData(ped, symbolData)
    labelData    <- getLabelData(symbolData, labeltext, labelsize, labeloffset)

    # BUILD PLOT ------------------------------------------------------------

    # Many ggplot layers are required to generate the desired plot. Symbols,
    # line segments and text labels are plotted on separate layers. In
    # addition, symbol fill and border are on separate layers. Only appropriate
    # layers contribute to the plot legend/key. Mock layers are used either for
    # their contribution to the legend or to mask certain elements from layers
    # below.

    # To minimize differences in the apparent size of various symbol shapes,
    # the shape and size are transformed prior to being passed as ggplot
    # aesthetics. The integer shape values 1-4 are mapped to integer codes for
    # plot symbols. Size values are adjusted by a shape-dependent value.

    # Base plot and theme
    plt <- ggplot() +
           theme(panel.background = element_rect(fill = bgcolor),
                 panel.grid = element_blank(),
                 axis.title = element_blank(),
                 axis.text  = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line  = element_blank())

    # Add a fully transparent mock layer to get a border key with lines
    plt <- plt +
           geom_line(data = symbolData,
                     mapping = aes(x = x, y = y, color = border),
                     size = 3, alpha = 0, linetype = 1)

    # Plot symbol borders first so they appear UNDER relationship lines
    plt <- plt +
           geom_point(data = symbolData,       # symbol border
                      mapping = aes(x     = x,
                                    y     = y,
                                    shape = as.integer(c(15, 16, 18, 17)[shape]),
                                    color = border,
                                    size  = as.numeric(size * c(1.25, 1.4, 1.46, 1.02)[shape]),
                                    alpha = alpha),
                      show_guide = FALSE) +
           geom_point(data = symbolData,       # spacer between border and fill
                      mapping = aes(x     = x,
                                    y     = y,
                                    shape = as.integer(c(15, 16, 18, 17)[shape]),
                                    size  = as.numeric(size * c(1.05, 1.2, 1.26, 0.82)[shape]),
                                    alpha = alpha),
                      color = bgcolor,
                      show_guide = FALSE)

    # Plot line segments showing relationships
    plt <- plt +
           geom_segment(data = relationData,
                        mapping = aes(x     = x,
                                      xend  = xend,
                                      y     = y,
                                      yend  = yend),
                        alpha = 1,
                        color = "gray20",
                        show_guide = FALSE)

    # Plot symbols showing individuals
    plt <- plt +
           geom_point(data = symbolData,       # mask lines under symbols
                      mapping = aes(x     = x,
                                    y     = y,
                                    shape = as.integer(c(22, 21, 23, 24)[shape]),
                                    size  = as.numeric(size * c(1, 1, 0.8, 0.6)[shape])),
                      fill = bgcolor,
                      color = bgcolor,
                      alpha = 1,
                      show_guide = FALSE) +
           geom_point(data = symbolData,       # symbol fill and outline
                      mapping = aes(x     = x,
                                    y     = y,
                                    shape = as.integer(c(22, 21, 23, 24)[shape]),
                                    fill  = fill,
                                    size  = as.numeric(size * c(1, 1, 0.8, 0.6)[shape]),
                                    alpha = alpha),
                      color = "gray20")

    # Plot line segments showing status (alive/censored or deceased)
    if(nrow(statusData) > 0) {
        plt <- plt +
               geom_segment(data = statusData,
                            mapping = aes(x     = x,
                                          xend  = xend,
                                          y     = y,
                                          yend  = yend,
                                          alpha = alpha),
                            color = "gray20",
                            show_guide = FALSE)
    }

    # Plot text labels under symbols
    plt <- plt +
           geom_text(data    = labelData,
                     mapping = aes(x     = x,
                                   y     = y,
                                   vjust = vjust,
                                   label = label,
                                   size  = size,
                                   alpha = alpha),
                     color   = "black",
                     hjust   = 0.5,
                     show_guide = FALSE)

    # Set the scale for ggplot aesthetics
    plt <- plt +
           scale_x_continuous(limits = range(symbolData$x)
                                       + c(-0.2, 0.2)
                                       * 0.1 * mean(symbolData$size)) +
           scale_y_continuous(trans = reverse_trans(),
                              limits = rev(range(symbolData$y)
                                           + c(-0.4, 0.2 + 0.4 * labelsize)  # Accomodate text labels
                                           * 0.1 * mean(symbolData$size))) + # Scale to symbol size
           scale_shape_identity() +
           switch(ifelse(is.factor(symbolAttrs$fill), 1, 2),
                  scale_fill_manual(values = attrValues$fill, drop = FALSE),
                  scale_fill_gradient(low  = attrValues$fill["low"],
                                      high = attrValues$fill["high"])) +
           switch(ifelse(is.factor(symbolAttrs$border), 1, 2),
                  scale_color_manual(values = attrValues$border, drop = FALSE),
                  scale_color_gradient(low  = attrValues$border["low"],
                                       high = attrValues$border["high"])) +
           scale_size_identity() +
           switch(ifelse(is.factor(symbolAttrs$alpha), 1, 2),
                  scale_alpha_manual(values = attrValues$alpha, drop = FALSE),
                  scale_alpha_continuous(range = c(low  = attrValues$alpha["low"],
                                                   high = attrValues$alpha["high"])))

    # Add legend matching symbol attributes
    gds <- list(shape = guide_legend(title = attrKeys$shape,
                                     order = 1),
                fill  = guide_legend(title = attrKeys$fill,
                                     override.aes = list(shape = 21,
                                                         size = mean(symbolAttrs$size)/2),
                                     order = 2),
                color = guide_legend(title = attrKeys$border,
                                     override.aes = list(size = 1, alpha = 1),
                                     order = 3),
                size  = guide_legend(title = attrKeys$size,
                                     override.aes = list(shape = 21),
                                     order = 4),
                alpha = guide_legend(title = attrKeys$alpha,
                                     override.aes = list(shape = 21,
                                                         size = mean(symbolAttrs$size)/2),
                                     order = 5))
    if(!is.null(attrKeys$shape)  && is.na(attrKeys$shape))  gds$shape <- FALSE
    if(!is.null(attrKeys$fill)   && is.na(attrKeys$fill))   gds$fill  <- FALSE
    if(!is.null(attrKeys$border) && is.na(attrKeys$border)) gds$color <- FALSE
    if(!is.null(attrKeys$size)   && is.na(attrKeys$size))   gds$size  <- FALSE
    if(!is.null(attrKeys$alpha)  && is.na(attrKeys$alpha))  gds$alpha <- FALSE
    plt <- plt +
           do.call(guides, gds) +
           theme(legend.key = element_rect(fill = bgcolor),
                 legend.position = "right")

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
    attrNames <- names(symbolAttrs)
    names(attrNames) <- attrNames
    # Assign mock default values to unused (NULL-valued) symbol attributes
    symbolAttrs <- lapply(attrNames,
                          function(name) {
                              if(is.null(symbolAttrs[[name]]))
                                  return(defaultSymbolAttr(name, n))
                              else
                                  return(symbolAttrs[[name]])
                          })
    # Check that attributes are valid
    symbolAttrs <- lapply(attrNames,
                          function(name) {
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
                          stop(sprintf(paste("Internal error: '%s' is an",
                                             "invalid 'name' argument; must",
                                             "be one of: shape, fill, border",
                                             "size, and alpha"),
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
                         stop(sprintf(paste("symbol%s must be a vector of",
                                            "type factor or numeric (or",
                                            "character, logical or integer)"),
                                      name)))
    if(length(symbolAttr) != n)
        stop(sprintf("symbol%s has length that does not match 'ped'", name))
    # if(any(is.na(symbolAttr)))
    #     stop(sprintf("%s: %s has NA(s), not allowed",
    #                  msgPrefix, symbolAttr))
    return(symbolAttr)
}

#' Return valid attribute values
#'
#' @param attrValues
#'      Named list of attribute values, with elements shape, fill, border,
#'      size and alpha
#' @inheritParams validSymbolAttrs
#' @param bgColor
#'      Plot background color
#' @return
#'      Named list of valid attribute values, with NULL values replaced with
#'      suitable defaults
validAttrValues <- function(attrValues, symbolAttrs, bgColor) {
    attrNames <- names(symbolAttrs)
    names(attrNames) <- attrNames
    attrValues <- lapply(attrNames,
                         function(name) {
                             if(is.null(attrValues[[name]]))
                                 return(defaultAttrValue(name,
                                                         symbolAttrs[[name]],
                                                         bgColor))
                             else
                                 return(attrValues[[name]])
                         })
    attrValues <- lapply(attrNames,
                         function(name) {
                             checkAttrValue(attrValues[[name]],
                                            symbolAttrs[[name]], name)
                         })
    return(attrValues)
}

#' Default attribute values
#'
#' @param name
#'      String. One of: fill, border, size, and alpha
#' @param symbolAttr
#'      Factor or numeric vector. Symbol attributes.
#' @inheritParams validAttrValues
#' @return
#'      Character vector. Default value of attribute values suitable for
#'      {fill,border,size,alpha}values
defaultAttrValue <- function(name, symbolAttr, bgColor) {
    if(is.factor(symbolAttr)) {
        if(nlevels(symbolAttr) > 3)
            stop(sprintf(paste("'symbol%s' argument has more than 3 levels.",
                               "Specify '%svalues'"),
                         name, name))
        attrValue <- switch(name,
                            "fill"   = switch(nlevels(symbolAttr),
                                              c("white"),
                                              c("white", "gray25"),
                                              c("white", "gray65", "gray25")),
                            "border" = switch(nlevels(symbolAttr),
                                              c(bgColor),
                                              c(bgColor, "#C66E00"),
                                              c(bgColor, "#009DCF", "#C66E00")),
                            "size"   = c(1, 2, 3)[1:nlevels(symbolAttr)],
                            "alpha"  = switch(nlevels(symbolAttr),
                                              c(1),
                                              c(1, 0.2),
                                              c(1, 0.6, 0.2)),
                            stop(sprintf(paste("Internal error: '%s' is an",
                                               "invalid 'name' argument; must",
                                               "be one of: shape, fill,",
                                               "border, size, and alpha"),
                                         name)))
        names(attrValue) <- levels(symbolAttr)
    } else if(is.numeric(symbolAttr)) {
        attrValue <- switch(name,
                            "fill"   = c(low = "white",   high = "black"),
                            "border" = c(low = "#fee5d9", high = "#a50f15"),
                            "size"   = c(low = 3,         high = 15),
                            "alpha"  = c(low = 0.2,       high = 1),
                            stop(sprintf(paste("Internal error: '%s' is an",
                                               "invalid 'name' argument; must",
                                               "be one of: shape, fill,",
                                               "border, size, and alpha"),
                                         name)))
    } else {
        stop(sprintf(paste("Internal error: symbolAttr is neither a factor",
                           "nor a numeric vector, but a", class(symbolAttr)),
                     name))
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
                stop(sprintf(paste("%svalues must be a named vector, with",
                                   "values for each level of the factor",
                                   "symbol%s"),
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
                                   "values named 'low' and 'high' to scale",
                                   "symbol%s"),
                             name, name))
            }
        }
    } else {
        stop(sprintf(paste("Internal error: symbol%s must be a factor or",
                           "numeric vector"),
                     name))
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

#' Return valid attribute key modifiers
#'
#' @param attrKeys
#'      Named list of attribute keys, with elements shape, fill, border, size
#'      and alpha
#' @inheritParams validSymbolAttrs
#' @return
#'      Named list of attribute keys
validAttrKeys <- function(attrKeys, symbolAttrs) {
    if(!setequal(names(attrKeys), names(symbolAttrs)))
        stop("Internal error: invalid attrKeys")
    attrNames <- names(attrKeys)
    names(attrNames) <- attrNames
    attrKeys <- lapply(attrNames,
                       function(x) {
                           y <- attrKeys[[x]]
                           if(is.null(symbolAttrs[[x]]) || is.null(y))
                               return(NA_character_)
                           else if((is.character(y) && length(y) == 1))
                               if(nchar(y) == 0)
                                   return(NULL)
                               else
                                   return(y)
                           else
                               stop(sprintf(paste("%skey must be a single",
                                                  "string value or NULL"),
                                            x))
                        })
    return(attrKeys)
}

#' Scale linearly between two bounds
#'
#' @param x
#'      Numeric values to be scaled
#' @param newMin,newMax
#'      Lower and upper bounds of the scaled values
#' @return
#'      Numeric values linearly scaled between min and max
scaleMinMax <- function(x, newMin, newMax) {
    if(max(x) > min(x))
        return(newMin + (newMax - newMin) * (x - min(x)) / (max(x) - min(x)))
    else
        return(rep((newMin + newMax) / 2, length(x)))
}

#' Construct line segment plotting data to show relationships
#'
#' @inheritParams ggpedigree
#' @return
#'      dataframe to be used as ggplot data for relation layer
getRelationData <- function(ped, pedalign) {
    n      <- pedalign$n
    nid    <- pedalign$nid
    spouse <- pedalign$spouse
    fam    <- pedalign$fam
    pos    <- pedalign$pos
    segmts <- list()
    for(i in seq_along(n)) {
        for(j in 1:n[i]) {
            if(spouse[i, j] == 1) {
                x      <- pos[i, j]
                xend   <- pos[i, j+1]
                y      <- i
                yend   <- i
                segmts <- append(segmts, list(c(x, xend, y, yend)))
            }
            if(fam[i, j] != 0) {
                js <- which(fam[i, ] == fam[i, j])
                if(j == js[1]) {
                    xmin <- pos[i, min(js)]
                    xmax <- pos[i, max(js)]
                    jdad <- which(nid[i-1, ] == ped$findex[nid[i, j]])
                    jmom <- which(nid[i-1, ] == ped$mindex[nid[i, j]])
                    xlparent <- min(pos[i-1, jdad], pos[i-1, jmom])
                    xrparent <- max(pos[i-1, jdad], pos[i-1, jmom])
                    xmparent <- (xlparent + xrparent) / 2

                    segmts <- append(segmts,
                        list(c(xmparent, xmparent, i-1, i-0.5)))
                    for(k in js) {
                        segmts <- append(segmts,
                                         list(c(pos[i, k], pos[i, k],
                                                i, i-0.5)))
                    }
                    segmts <- append(segmts,
                        list(c(min(xmin, xmparent), max(xmax, xmparent),
                               i-0.5, i-0.5)))
                }
            }
        }
    }
    relationData <- as.data.frame(t(do.call(data.frame, segmts)))
    names(relationData) <- c("x", "xend", "y", "yend")
    return(relationData)
}

#' Construct symbol plotting data
#'
#' @inheritParams ggpedigree
#' @inheritParams validSymbolAttrs
#' @return
#'      dataframe to be used as ggplot data for symbol layers
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
    symbolData <- do.call(data.frame,
                          c(list(id = ped$id, x = xpos, y = ypos),
                            symbolAttrs))
    return(symbolData)
}

#' Construct symbol status plotting data to indicate deceased individuals
#'
#' @inheritParams ggpedigree
#' @param symbolData
#'      dataframe with symbol plotting data
#' @return
#'      dataframe to be used as ggplot data for status layer
getStatusData <- function(ped, symbolData) {
    dead <- ped$id[ped$status == 1]
    symbolData <- symbolData[symbolData$id %in% dead,
                             c("id", "x", "y", "size", "alpha")]
    offset <- 0.06 * symbolData$size/2 * cos(pi/4)
    statusData <- data.frame(id    = symbolData$id,
                             x     = symbolData$x - offset,
                             xend  = symbolData$x + offset,
                             y     = symbolData$y + offset,
                             yend  = symbolData$y - offset,
                             alpha = symbolData$alpha)
    return(statusData)
}

#' Construct label plotting data to annotate symbols
#'
#' @inheritParams getStatusData
#' @param labelText
#'      character vector of labels used in addition to ped$id for annotation
#' @param labelSize
#'      single numeric value used as size multiplier
#' @param labelOffset
#'      single numeric value used to adjust vertical offset of labels
#' @return
#'      dataframe to be used as ggplot data for label layer
getLabelData <- function(symbolData, labelText, labelSize, labelOffset) {
    labelData <- data.frame(label = as.character(symbolData$id),
                            x     = symbolData$x,
                            y     = symbolData$y,
                            vjust = 3.5 - 0.1 * symbolData$size + labelOffset,
                            size  = symbolData$size / 2 * labelSize,
                            alpha = symbolData$alpha)
    if(!is.null(labelText))
        labelData$label <- paste(labelData$label, labelText, sep = "\n")
    return(labelData)
}

# TODO: support the case when an individual is shown more than once
