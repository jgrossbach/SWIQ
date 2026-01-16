add_line <- function(p,
                     slope,
                     offset,
                     color = "red",
                     linewidth = 1) {
  p +
    geom_abline(
      slope = slope,
      intercept = offset,
      color = color,
      linewidth = linewidth
    )
}


#####histogramm with ggplot#####################################################
gg_histogram <- function(x, breaks = 30, color = "steelblue",
                         xlab = "Value", ylab = "Count", title = "Histogram") {
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(x = x)
  
  # Build and return the plot
  ggplot(df, aes(x = x)) +
    geom_histogram(bins = breaks, fill = color, color = "black", alpha = 0.8) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 14)
}

#####barplot with ggplot########################################################
gg_barplot <- function(values, 
                       title = "Barplot", 
                       xlab = "Index", 
                       ylab = "Value", 
                       bar_color = "steelblue") {
  
  # Check if input is numeric
  if(!is.numeric(values)) {
    stop("Input 'values' must be a numeric vector.")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(Index = seq_along(values), Value = values)
  
  # Create barplot using ggplot2
  p <- ggplot(df, aes(x = factor(Index), y = Value)) +
    geom_bar(stat = "identity", fill = bar_color) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

gg_boxplot <- function(x, y, 
                       title = NULL, 
                       xlab = NULL, 
                       ylab = NULL, 
                       fill = "steelblue",
                       box_labels = NULL,
                       show_outliers = TRUE) {
  # Basic checks
  if (missing(x) | missing(y)) {
    stop("Please provide both x (grouping) and y (numeric vector).")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric.")
  }
  if (!is.factor(x) && !is.character(x)) {
    x <- factor(x)
  }
  
  # Build data frame internally
  df <- data.frame(group = x, value = y)
  
  # Determine y-limits if outliers are hidden
  ylim_vals <- NULL
  if (!show_outliers) {
    # Compute the lower and upper whisker values for each group
    library(dplyr)
    whiskers <- df %>%
      group_by(group) %>%
      summarise(
        lower = quantile(value, 0.25) - 1.5 * IQR(value),
        upper = quantile(value, 0.75) + 1.5 * IQR(value),
        .groups = "drop"
      )
    ylim_vals <- c(min(whiskers$lower), max(whiskers$upper))
  }
  
  # Create plot
  p <- ggplot(df, aes(x = group, y = value)) +
    geom_boxplot(
      fill = fill,
      outlier.shape = if (show_outliers) 19 else NA
    ) +
    theme_minimal()
  
  # Apply adjusted y-limits if needed
  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals)
  }
  
  # Optional labels
  if (!is.null(title)) p <- p + ggtitle(title)
  
  # Use provided box_labels for the x-axis scale breaks
  if (!is.null(box_labels)) {
    if (length(box_labels) != length(levels(x))) {
      stop("The length of 'box_labels' must match the number of groups (levels of x).")
    }
    p <- p + scale_x_discrete(labels = box_labels)
    if (is.null(xlab)) {
      xlab <- "Group"
    }
  }
  
  if (!is.null(xlab))  p <- p + xlab(xlab)
  if (!is.null(ylab))  p <- p + ylab(ylab)
  
  return(p)
}


named_list <- function(name_vec){
  l <- lapply(name_vec,get)
  names(l) <- name_vec
  return(l)
}

gg_dotplot <- function(x, y,
                       col = "black",
                       main = NULL,
                       pch = 19,
                       xlab = NULL,
                       ylab = NULL,
                       xlim = NULL,
                       ylim = NULL) {
  
  # Load required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  # Build data frame
  df <- data.frame(x = x, y = y)
  
  # Start plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = col, shape = pch)
  
  # Labels
  if (!is.null(main)) p <- p + ggplot2::ggtitle(main)
  if (!is.null(xlab)) p <- p + ggplot2::xlab(xlab)
  if (!is.null(ylab)) p <- p + ggplot2::ylab(ylab)
  
  # Axis limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  
  return(p)
}
