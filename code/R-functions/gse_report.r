gse.dotplot <- function(object, x = "p.adjust", y = "Description", color = "p.adjust", size = "setSize",
                        x2log10 = FALSE, color2log10 = FALSE,
                        showCategory = 10, split = NULL,
                        title = "",
                        orderBy = "x", label_format = 30, decreasing = FALSE,
                        xlab = x,
                        scale_size_range = c(3, 8),
                        font_size_y = 8,
                        dot_alpha = 1,
                        vline = NULL) {
  df <- data.frame(
    x = as.numeric(object[, x]),
    y = object[, y],
    color = object[, color],
    size = object[, size],
    stringsAsFactors = FALSE
  ) %>%
    head(showCategory)

  if (!orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }

  label_func <- default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }

  # Set color
  use_color_scale <- scale_color_viridis_c()
  color_title <- color

  size_title <- size


  if (x2log10) {
    df$x <- -log10(df$x)
    xlab <- bquote("-log"[10] ~ .(paste0("(", x, ")")))
  }
  if (color2log10) {
    df$color <- -log10(df$color)
    use_color_scale <- scale_color_viridis_c(direction = -1)
    color_title <- bquote("-log"[10] ~ .(paste0("(", color, ")")))
  }


  # Set order Y-axis
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$y <- factor(df$y, levels = rev(unique(df$y[idx])))

  res <- ggplot(df, aes(x = x, y = y, size = size, color = color)) +
    geom_point(alpha = dot_alpha) +
    use_color_scale +
    scale_y_discrete(labels = label_func) +
    labs(
      x = xlab,
      y = NULL,
      title = title,
      color = color_title,
      size = size_title
    ) +
    scale_size(range = scale_size_range) +
    theme(axis.text.y = element_text(size = font_size_y))

  if (!is.null(vline)) {
    res <- res + geom_vline(xintercept = vline, lty = 2)
  }

  return(res)
}




##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.

ep_str_wrap <- function(string, width) {
  x <- gregexpr(" ", string)
  vapply(seq_along(x),
    FUN = function(i) {
      y <- x[[i]]
      n <- nchar(string[i])
      len <- (c(y, n) - c(0, y)) ## length + 1
      idx <- len > width
      j <- which(!idx)
      if (length(j) && max(j) == length(len)) {
        j <- j[-length(j)]
      }
      if (length(j)) {
        idx[j] <- len[j] + len[j + 1] > width
      }
      idx <- idx[-length(idx)] ## length - 1
      start <- c(1, y[idx] + 1)
      end <- c(y[idx] - 1, n)
      words <- substring(string[i], start, end)
      paste0(words, collapse = "\n")
    },
    FUN.VALUE = character(1)
  )
}

default_labeller <- function(n) {
  function(str) {
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

# Object gostplot_pal : color palette for gprofiler2::gostplot
gostplot_pal <-  c(
  `GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618",
  MF = "#dc3912", BP = "#ff9900", CC = "#109618",
  KEGG = "#dd4477", REAC = "#3366cc", WP = "#0099c6",
  TF = "#5574a6", MIRNA = "#22aa99", HPA = "#6633cc",
  CORUM = "#66aa00", HP = "#990099",
  msigdb.c2.cp = "#3355cc", msigdb.c2.cgp = "#6633cc", msigdb.h = "#66aa00", msigdb.c3 = "#990099"
)

# Function gse.gostplot
# Adapted from gprofiler2::gostplot
gse.gostplot <- function(object, capped = TRUE, interactive = TRUE, ylab = bquote("-log"[10] ~ .("(p-adj)")),
                         pal = c(
                           `GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618",
                           MF = "#dc3912", BP = "#ff9900", CC = "#109618", KEGG = "#dd4477",
                           REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99",
                           HPA = "#6633cc", CORUM = "#66aa00", HP = "#990099",
                           msigdb.c2.cp = "#6633cc", msigdb.h = "#66aa00", msigdb.c3 = "#990099"
                         )) {
  if (is.null(pal)) {
    pal <- c(
      `GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618",
      MF = "#dc3912", BP = "#ff9900", CC = "#109618",
      KEGG = "#dd4477", REAC = "#3366cc", WP = "#0099c6",
      TF = "#5574a6", MIRNA = "#22aa99", HPA = "#6633cc",
      CORUM = "#66aa00", HP = "#990099",
      msigdb.c2.cp = "#6633cc", msigdb.h = "#66aa00", msigdb.c3 = "#990099"
    )
  }

  df <- object$result %>% filter(significant)
  meta <- object$meta
  source_order <- logpval <- term_id <- opacity <- NULL
  term_size <- term_name <- p_value <- term_size_scaled <- NULL

  essential_names <- c(
    "source_order", "term_size", "term_name",
    "term_id", "source", "significant"
  )
  if (!(all(essential_names %in% colnames(df)))) {
    stop(paste(
      "The following columns are missing from the result:",
      paste0(setdiff(essential_names, colnames(df)), collapse = ", ")
    ))
  }
  if (!any(grepl("p_value", colnames(df)))) {
    stop("Column 'p_value(s)' is missing from the result")
  }
  widthscale <- unlist(lapply(
    meta$query_metadata$sources,
    function(x) meta$result_metadata[[x]][["number_of_terms"]]
  ))
  names(widthscale) <- meta$query_metadata$sources
  space <- 1000
  starts <- c()
  start <- 1
  starts[1] <- start
  if (!length(widthscale) < 2) {
    for (idx in 2:length(widthscale)) {
      starts[idx] <- starts[idx - 1] + space + widthscale[idx -
        1]
    }
  }
  names(starts) <- names(widthscale)
  if (is.null(names(pal))) {
    names(pal) <- meta$query_metadata$sources[1:length(pal)]
  }
  sourcediff <- setdiff(meta$query_metadata$sources, names(pal))
  colors <- grDevices::colors(distinct = TRUE)[grep("gr(a|e)y|white|snow|khaki|lightyellow",
    grDevices::colors(distinct = TRUE),
    invert = TRUE
  )]
  if (length(sourcediff) > 0) {
    use_cols <- sample(colors, length(sourcediff), replace = FALSE)
    pal[sourcediff] <- use_cols
  }
  if ("p_values" %in% colnames(df)) {
    p_values <- query <- significant <- NULL
    # df$query <- list(names(meta$query_metadata$queries))
    df <- tidyr::unnest(data = df, cols = c(
      p_values, query,
      significant
    ))
    df <- dplyr::rename(df, p_value = p_values)
  }
  logScale <- function(input, input_start = 1, input_end = 50000,
                       output_start = 2, output_end = 10) {
    m <- (output_end - output_start) / (log(input_end) - log(input_start))
    b <- -m * log(input_start) + output_start
    output <- m * log(input) + b
    return(output)
  }
  xScale <- function(input, input_start = 1, input_end = sum(widthscale) +
                       (length(widthscale) - 1) * space, output_start = 2, output_end = 200) {
    m <- (output_end - output_start) / (input_end - input_start)
    b <- -m * input_start + output_start
    output <- m * input + b
    return(output)
  }
  df$logpval <- -log10(df$p_value)
  df$opacity <- ifelse(df$significant, 0.8, ifelse(df$p_value ==
    1, 0, 0.3))
  df$term_size_scaled <- logScale(df$term_size)
  df <- df %>%
    dplyr::group_by(source) %>%
    dplyr::mutate(order = xScale(source_order,
      input_start = 1, input_end = widthscale[source], output_start = starts[source],
      output_end = starts[source] + widthscale[source]
    ))
  df$order <- xScale(df$order)
  if (capped) {
    df$logpval[df$logpval > 16] <- 17
    ymin <- -1
    ymax <- 18.5
    ticklabels <- c(
      "0", "2", "4", "6", "8", "10", "12",
      "14", ">16"
    )
    tickvals <- c(0, 2, 4, 6, 8, 10, 12, 14, 16)
  }
  else {
    ymin <- -1
    ymax <- ceiling(max(df$logpval)) + 5
    ticklabels <- ggplot2::waiver()
    tickvals <- ggplot2::waiver()
  }
  if (interactive) {
    sd <- crosstalk::SharedData$new(df, key = ~term_id)
  }
  else {
    sd <- df
  }
  p <- ggplot2::ggplot(data = sd, ggplot2::aes(
    x = order, y = logpval,
    text = paste(
      term_id, paste0("(", term_size, ")"), "<br>",
      term_name, "<br>", formatC(p_value,
        format = "e",
        digits = 3
      )
    )
  )) +
    ggplot2::geom_point(ggplot2::aes(color = source, size = term_size_scaled, alpha = opacity), show.legend = FALSE) +
    ggplot2::facet_wrap(~query, ncol = 1, scales = "free_x", shrink = FALSE) +
  ggplot2::labs(y = ylab) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    panel.border = ggplot2::element_blank(), strip.text = ggplot2::element_text(
      size = 12,
      colour = "darkgrey"
    ), strip.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(
      size = 8,
      angle = 45, hjust = 1
    ), axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(
      color = "grey",
      size = 0.5
    ), axis.line.x = ggplot2::element_line(
      color = "grey",
      size = 0.1
    ), axis.line.y = ggplot2::element_line(
      size = 0.5,
      color = "grey"
    ), strip.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0, size = 10
    ), plot.margin = ggplot2::margin(
      t = 0,
      r = 5, b = 20, l = 20, unit = "pt"
    ), axis.title.y = ggplot2::element_text(
      size = 10,
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)
    )
  ) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::scale_alpha(range = c(
    0,
    0.8
  ), limits = c(0, 0.8)) +
  ggplot2::scale_y_continuous(expand = c(
    0,
    0
  ), limits = c(ymin, ymax), labels = ticklabels, breaks = tickvals) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(
    0,
    210
  ), breaks = (xScale(starts) + xScale(starts +
    widthscale)) / 2, labels = names(widthscale))
  for (s in names(widthscale)) {
    xstart = xScale(starts[s])
    xend = xScale(starts[s] + widthscale[s])
    p <- p + ggplot2::annotate("segment", x = xstart, xend = xend,
                               y = -1, yend = -1, size = 3, colour = pal[s])
  }
  if (capped) {
    p <- p + ggplot2::annotate(geom = "text", x = 180, y = 16.2,
                               label = "values above this threshold are capped",
                               size = 2, color = "grey") + ggplot2::geom_hline(yintercept = 16,
                                                                               linetype = "dashed", size = 0.2, color = "grey")
  }
  if (interactive) {
    p <- p %>% plotly::ggplotly(tooltip = "text")
    p <- p %>% plotly::highlight(on = "plotly_click", off = "plotly_doubleclick",
                                 dynamic = FALSE, persistent = FALSE)
  }
  return(p)
}

# Function gost.topGO
# Convert topGO returned by gse_omnibus into a list that can be used by gse.gostplot
gost.topGO <- function(x) {
  queries <- names(x)
  sources <- map(x, names) %>%
    unlist() %>%
    unique()

  qi <- queries[[1]]
  si <- sources[[1]]

  res_result <- foreach(qi = queries, .combine = rbind) %do% {
    foreach(si = sources, .combine = rbind) %do% {
      terms_sorted <- x[[qi]][[si]]$GOresult@score %>%
        names() %>%
        sort()
      x[[qi]][[si]]$GOtab %>%
        mutate(
          source = si,
          query = qi,
          p_value = elimFisher,
          term_size = Annotated,
          term_name = Term,
          term_id = GO.ID,
          significant = p_value <= 0.05,
          source_order = match(term_id, terms_sorted)
        )
    }
  }

  res_meta <- list()
  res_meta$query_metadata <- list(sources = sources)
  res_meta$result_metadata <- foreach(si = sources) %do% {
    list(
      number_of_terms = x[[1]][[si]]$GOresult@score %>% length()
    )
  }
  names(res_meta$result_metadata) <- sources

  list(
    result = res_result,
    meta = res_meta
  )
}


# Function gost.enricher
# Convert enricher returned by gse_omnibus into a list that can be used by gse.gostplot
gost.enricher <- function(x, terms = NULL, sources = NULL, simplify = FALSE) {
  queries <- names(x)
  if (is.null(sources)) {
    sources <- map(x, names) %>%
      unlist() %>%
      unique()
  }
  qi <- queries[[1]]
  si <- sources[[1]]

  res_result <- foreach(qi = queries, .combine = rbind) %do% {
    if (class(x[[qi]]) %in% c("enrichResult", "gseaResult")) {
      if (is.null(terms)) {
        terms_sorted <- x[[qi]]@result$ID %>% sort()
      } else {
        terms_sorted <- terms[[si]]
      }
      x[[qi]]@result %>%
        mutate(
          source = si,
          query = qi,
          p_value = p.adjust,
          term_size = setSize,
          term_name = Description,
          term_id = ID,
          significant = p.adjust <= 0.05,
          source_order = match(ID, terms_sorted)
        )
    } else {
      foreach(si = sources, .combine = rbind) %do% {
        if (is.null(terms)) {
          terms_sorted <- x[[qi]][[si]]@result$ID %>% sort()
        } else {
          terms_sorted <- terms[[si]]
        }
        x[[qi]][[si]]@result %>%
          mutate(
            source = si,
            query = qi,
            p_value = p.adjust,
            term_size = setSize,
            term_name = Description,
            term_id = ID,
            significant = p.adjust <= 0.05,
            source_order = match(ID, terms_sorted)
          )
      }
    }
  }

  if (simplify & ("simplify" %in% names(res_result))) {
    res_result <- res_result %>% filter(simplify)
  }
  if ("simplify" %in% names(res_result)) {
    res_result$simplify <- NULL
  }
  res_meta <- list()
  res_meta$query_metadata <- list(sources = sources)

  if (length(sources) > 1) {
    res_meta$result_metadata <- foreach(si = sources) %do% {
      list(
        number_of_terms = ifelse(is.null(terms),
          nrow(x[[1]][[si]]@result),
          length(terms[[si]])
        )
      )
    }
  } else {
    res_meta$result_metadata[[sources]] <- list(
      number_of_terms = ifelse(is.null(terms),
        nrow(x[[1]]@result),
        length(terms[[sources]])
      )
    )
  }
  names(res_meta$result_metadata) <- sources

  list(
    result = res_result,
    meta = res_meta
  )
}


# Function gost.enricher.combine
# Combine objects from gost.enricher. Input variable should be a list (named or unamed) of outputs from gost.enricher. Should work also for other gost.** functions
gost.enricher.combine <- function(x) {
  res_result <- foreach(i = x, .combine = rbind) %do% {
    i$result
  }
  res_meta <- list()
  res_meta$query_metadata$sources <- map(x, function(xi) xi$meta$query_metadata$sources) %>% unlist()
  res_meta$result_metadata <- list()
  for (xi in x) {
    res_meta$result_metadata <- append(res_meta$result_metadata, xi$meta$result_metadata)
  }

  list(
    result = res_result,
    meta = res_meta
  )
}

# Function gost.edgeR_gse
# Convert returned df from edgeR analysis (camera, fry...) returned by gse_edger_omnibus into a list that can be used by gse.gostplot
gost.edgeR_gse <- function(x, terms = NULL, sources = NULL) {
  queries <- names(x)
  if (is.null(sources)) {
    sources <- map(x, names) %>%
      unlist() %>%
      unique()
  }
  qi <- queries[[1]]
  si <- sources[[1]]

  res_result <- foreach(qi = queries, .combine = rbind) %do% {
    foreach(si = sources, .combine = rbind) %do% {
      if (is.null(terms)) {
        terms_sorted <- x[[qi]][[si]]$Term %>% sort()
      } else {
        terms_sorted <- terms[[si]]
      }
      x[[qi]][[si]] %>%
        mutate(
          source = si,
          query = qi,
          p_value = FDR,
          term_size = NGenes,
          term_name = Term,
          term_id = Term,
          significant = FDR <= 0.05,
          source_order = match(Term, terms_sorted)
        )
    }
  }

  res_meta <- list()
  res_meta$query_metadata <- list(sources = sources)

  if (length(sources) > 1) {
    res_meta$result_metadata <- foreach(si = sources) %do% {
      list(
        number_of_terms = ifelse(is.null(terms),
                                 nrow(x[[1]][[si]]),
                                 length(terms[[si]])
        )
      )
    }
  } else {
    res_meta$result_metadata[[sources]] <- list(
      number_of_terms = ifelse(is.null(terms),
                               nrow(x[[1]]@result),
                               length(terms[[sources]])
      )
    )
  }
  names(res_meta$result_metadata) <- sources

  list(
    result = res_result,
    meta = res_meta
  )
}


# Function gost.gsva
# Convert returned df from gsva analysis returned by into a list that can be used by gostplot
gost.gsva <- function(x, method = 'gsva', terms = NULL, sources = NULL) {
  queries <- names(x)
  if (is.null(sources)) {
    sources <- map(x, names) %>%
      unlist() %>%
      unique()
  }
  qi <- queries[[1]]
  si <- sources[[1]]

  res_result <- foreach(qi = queries, .combine = rbind) %do% {
    foreach(si = sources, .combine = rbind) %do% {
      if (is.null(terms)) {
        terms_sorted <- x[[qi]][[si]][[method]]$gse$term_name %>% sort()
      } else {
        terms_sorted <- terms[[si]]
      }
      x[[qi]][[si]][[method]]$gse %>%
        mutate(
          source = si,
          query = qi,
          p_value = adj.P.Val,
          term_id = term_name,
          significant = adj.P.Val <= 0.05,
          Direction = ifelse(logFC > 0, 'Up', 'Down'),
          source_order = match(term_name, terms_sorted)
        )
    }
  }

  res_meta <- list()
  res_meta$query_metadata <- list(sources = sources)

  if (length(sources) > 1) {
    res_meta$result_metadata <- foreach(si = sources) %do% {
      list(
        number_of_terms = ifelse(is.null(terms),
                                 nrow(x[[1]][[si]][[method]]$gse),
                                 length(terms[[si]])
        )
      )
    }
  } else {
    res_meta$result_metadata[[sources]] <- list(
      number_of_terms = ifelse(is.null(terms),
                               nrow(x[[1]][[sources]][[method]]$gse),
                               length(terms[[sources]])
      )
    )
  }
  names(res_meta$result_metadata) <- sources

  list(
    result = res_result,
    meta = res_meta
  )
}
