# Function pca_eigencorplot
# Built upon PCAtools::eigencorplot. If returnPlot==FALSE, it returns a list witht the correlation and P/FDR values
pca_eigencorplot <- function (pcaobj, components = getComponents(pcaobj, seq_len(10)),
                metavars, titleX = "", cexTitleX = 1, rotTitleX = 0, colTitleX = "black",
                fontTitleX = 2, titleY = "", cexTitleY = 1, rotTitleY = 0,
                colTitleY = "black", fontTitleY = 2, cexLabX = 1, rotLabX = 0,
                colLabX = "black", fontLabX = 2, cexLabY = 1, rotLabY = 0,
                colLabY = "black", fontLabY = 2, posLab = "bottomleft", col = c("blue4",
                                                                                "blue3", "blue2", "blue1", "white", "red1", "red2", "red3",
                                                                                "red4"), posColKey = "right", cexLabColKey = 1, cexCorval = 1,
                colCorval = "black", fontCorval = 1, scale = TRUE, main = "",
                cexMain = 2, rotMain = 0, colMain = "black", fontMain = 2,
                corFUN = "pearson", corUSE = "pairwise.complete.obs", corMultipleTestCorrection = "none",
                signifSymbols = c("***", "**", "*", ""), signifCutpoints = c(0,
                                                                             0.001, 0.01, 0.05, 1), colFrame = "white", plotRsquared = FALSE,
                returnPlot = TRUE)
{
  require(PCAtools)
  require(lattice)
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata
  for (i in seq_len(length(components))) {
    if (!is.numeric(data[, components[i]])) {
      warning(components[i], " is not numeric - please check the source data",
              " as everything will be converted to a matrix")
    }
  }
  for (i in seq_len(length(metavars))) {
    if (!is.numeric(metadata[, metavars[i]])) {
      warning(metavars[i], " is not numeric - please check the source data",
              " as non-numeric variables will be coerced to numeric")
    }
  }
  xvals <- data.matrix(data[, which(colnames(data) %in% components)])
  yvals <- metadata[, which(colnames(metadata) %in% metavars)]
  chararcter_columns = unlist(lapply(yvals, is.numeric))
  chararcter_columns = !chararcter_columns
  chararcter_columns = names(which(chararcter_columns))
  for (c in chararcter_columns) {
    yvals[, eval(quote(c))] = as.numeric(as.factor(yvals[,
                                                         eval(quote(c))]))
  }
  yvals <- data.matrix(yvals)
  corvals <- cor(xvals, yvals, use = corUSE, method = corFUN)
  N <- ncol(xvals) * ncol(yvals)
  pvals <- data.frame(pval = numeric(N), i = numeric(N), j = numeric(N))
  k <- 0
  for (i in seq_len(ncol(xvals))) {
    for (j in seq_len(ncol(yvals))) {
      k <- k + 1
      pvals[k, "pval"] <- cor.test(xvals[, i], yvals[,
                                                     j], use = corUSE, method = corFUN)$p.value
      pvals[k, "i"] <- colnames(xvals)[i]
      pvals[k, "j"] <- colnames(yvals)[j]
    }
  }
  if (corMultipleTestCorrection != "none") {
    pvals$pval <- p.adjust(pvals$pval, method = corMultipleTestCorrection)
  }
  pvals <- reshape2::dcast(pvals, i ~ j, value.var = "pval")
  rownames(pvals) <- pvals$i
  pvals$i <- NULL
  pvals <- pvals[match(rownames(corvals), rownames(pvals)),
  ]
  pvals <- pvals[colnames(corvals)]

  if (plotRsquared == TRUE) {
    corvals <- corvals^2
  }
  if (scale == FALSE && plotRsquared == TRUE) {
    iUpperRange <- 1
    iLowerRange <- 0
  }
  else if (scale == FALSE && plotRsquared == FALSE) {
    iUpperRange <- 1
    iLowerRange <- -1
  }
  else if (scale == TRUE) {
    max <- max(corvals)
    min <- min(corvals)
    if (abs(max) > abs(min)) {
      iUpperRange <- max + 0.01
      iLowerRange <- (max * (-1)) - 0.01
    }
    else {
      iUpperRange <- abs(min) + 0.01
      iLowerRange <- min - 0.01
    }
    if (plotRsquared == TRUE) {
      iUpperRange <- max + 0.1
      iLowerRange <- 0
    }
  }
  cols <- colorRampPalette(col)
  signif <- corvals
  for (i in seq_len(ncol(pvals))) {
    signif[, i] <- c(symnum(pvals[, i], corr = FALSE, na = FALSE,
                            cutpoints = signifCutpoints, symbols = signifSymbols))
  }
  plotLabels <- corvals
  for (i in seq_len(nrow(corvals))) {
    for (j in seq_len(ncol(corvals))) {
      plotLabels[i, j] <- paste(round(corvals[i, j], 2),
                                signif[i, j], sep = "")
      colnames(plotLabels)[j] <- colnames(corvals)[j]
    }
    rownames(plotLabels)[i] <- rownames(corvals)[i]
  }
  if (posLab == "bottomleft") {
    posLab = 1
    axisTicks = c(1, 0)
  }
  else if (posLab == "topright") {
    posLab = 2
    axisTicks = c(0, 1)
  }
  else if (posLab == "all") {
    posLab = 3
    axisTicks = c(1, 1)
  }
  else if (posLab == "none") {
    posLab = 0
    axisTicks = c(0, 0)
  }
  labels <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    ltext(x, y, labels = plotLabels, cex = cexCorval, col = colCorval,
          font = fontCorval)
  }
  l <- levelplot(data.matrix(corvals), xlab = list(label = titleX,
                                                   cex = cexTitleX, rot = rotTitleX, col = colTitleX, font = fontTitleX),
                 ylab = list(label = titleY, cex = cexTitleY, rot = rotTitleY,
                             col = colTitleY, font = fontTitleY), panel = labels,
                 pretty = TRUE, par.settings = list(panel.background = list(col = colFrame)),
                 scales = list(x = list(cex = cexLabX, rot = rotLabX,
                                        col = colLabX, font = fontLabX), y = list(cex = cexLabY,
                                                                                  rot = rotLabY, col = colLabY, font = fontLabY), tck = axisTicks,
                               alternating = posLab), aspect = "fill", col.regions = cols,
                 cuts = 100, at = seq(iLowerRange, iUpperRange, 0.01),
                 main = list(label = main, cex = cexMain, rot = rotMain,
                             col = colMain, font = fontMain), colorkey = list(space = posColKey,
                                                                              labels = list(cex = cexLabColKey)))
  if (returnPlot == TRUE) {
    return(l)
  }
  else if (returnPlot == FALSE) {
    list(corvals = corvals, pvals =pvals)
  }
}
