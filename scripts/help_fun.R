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
                         xlab = "Value", ylab = "Count", title = "Histogram",
                         title_size = 14,
                         xlab_size = 12,
                         ylab_size = 12,
                         axis_text_size = 10) {
  
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  
  # Create a data frame for ggplot
  df <- data.frame(x = x)
  
  # Build and return the plot
  ggplot(df, aes(x = x)) +
    geom_histogram(
      bins = breaks,
      fill = color,
      color = "black",
      alpha = 0.8
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size),
      axis.title.x = element_text(size = xlab_size),
      axis.title.y = element_text(size = ylab_size),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size)
    )
}

#####barplot with ggplot########################################################
# gg_barplot <- function(values, 
#                        title = "Barplot", 
#                        xlab = "Index", 
#                        ylab = "Value", 
#                        bar_color = "steelblue") {
#   
#   # Check if input is numeric
#   if(!is.numeric(values)) {
#     stop("Input 'values' must be a numeric vector.")
#   }
#   
#   # Create a data frame for ggplot
#   df <- data.frame(Index = seq_along(values), Value = values)
#   
#   # Create barplot using ggplot2
#   p <- ggplot(df, aes(x = factor(Index), y = Value)) +
#     geom_bar(stat = "identity", fill = bar_color) +
#     labs(title = title, x = xlab, y = ylab) +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   return(p)
# }

gg_barplot <- function(values, 
                       bar_labels = NULL,
                       title = "Barplot", 
                       xlab = "Index", 
                       ylab = "Value", 
                       bar_color = "steelblue",
                       axis_text_size = 12,
                       title_size = 16,
                       xlab_size = 14,
                       ylab_size = 14) {
  
  # Check if input is numeric
  if (!is.numeric(values)) {
    stop("Input 'values' must be a numeric vector.")
  }
  
  # Use default labels if none are provided
  if (is.null(bar_labels)) {
    bar_labels <- seq_along(values)
  }
  
  # Check labels length
  if (length(bar_labels) != length(values)) {
    stop("Length of 'bar_labels' must match length of 'values'.")
  }
  
  # Create data frame for ggplot
  df <- data.frame(
    Label = factor(bar_labels, levels = bar_labels),
    Value = values
  )
  
  # Create barplot using ggplot2
  p <- ggplot(df, aes(x = Label, y = Value)) +
    geom_bar(stat = "identity", fill = bar_color) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = title_size
      ),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title.x = element_text(size = xlab_size),
      axis.title.y = element_text(size = ylab_size)
    )
  
  return(p)
}

gg_barplot_matrix <- function(x,
                              pair_labels = NULL,
                              title = "Barplot",
                              xlab = "Index",
                              ylab = "Value",
                              bar_color = c("steelblue", "tomato"),
                              legend_labels = NULL,
                              show_legend = TRUE,
                              axis_text_size = 12,
                              title_size = 16,
                              xlab_size = 14,
                              ylab_size = 14,
                              legend_text_size = 12,
                              legend_title = "Group",
                              within_pair_gap = 0.6,
                              between_pair_gap = 1.4,
                              bar_width = 0.5) {
  
  # Check input
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Input 'x' must be a numeric matrix.")
  }
  
  values <- as.vector(x)   # column-wise order
  
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  n_bars <- length(values)
  
  # Default pair labels
  if (is.null(pair_labels)) {
    pair_labels <- paste0("Group ", 1:n_cols)
  }
  
  if (length(pair_labels) != n_cols) {
    stop("Length of 'pair_labels' must match number of matrix columns.")
  }
  
  # Row-based color recycling
  row_colors <- rep(bar_color, length.out = n_rows)
  fill_colors <- rep(row_colors, times = n_cols)
  
  # Default legend labels
  if (is.null(legend_labels)) {
    legend_labels <- paste0("Row ", 1:n_rows)
  }
  
  if (length(legend_labels) != n_rows) {
    stop("Length of 'legend_labels' must match number of matrix rows.")
  }
  
  # Create custom x positions
  positions <- c()
  pair_centers <- c()
  row_group <- c()
  current_pos <- 1
  
  for (j in 1:n_cols) {
    pair_pos <- c()
    
    for (i in 1:n_rows) {
      positions <- c(positions, current_pos)
      pair_pos <- c(pair_pos, current_pos)
      row_group <- c(row_group, legend_labels[i])
      
      if (i < n_rows) {
        current_pos <- current_pos + within_pair_gap
      } else {
        current_pos <- current_pos + between_pair_gap
      }
    }
    
    pair_centers <- c(pair_centers, mean(pair_pos))
  }
  
  df <- data.frame(
    xpos = positions,
    Value = values,
    Row = factor(row_group, levels = legend_labels)
  )
  
  library(ggplot2)
  
  p <- ggplot(df, aes(x = xpos, y = Value, fill = Row)) +
    geom_col(width = bar_width) +
    scale_fill_manual(
      values = setNames(row_colors, legend_labels),
      name = legend_title
    ) +
    scale_x_continuous(
      breaks = pair_centers,
      labels = pair_labels
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = title_size),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      axis.title.x = element_text(size = xlab_size),
      axis.title.y = element_text(size = ylab_size),
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_text_size)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

gg_boxplot <- function(x, y, 
                       title = NULL, 
                       xlab = NULL, 
                       ylab = NULL, 
                       fill = "darkgreen",
                       box_labels = NULL,
                       show_outliers = TRUE,
                       title_size = 14,
                       xlab_size = 12,
                       ylab_size = 12,
                       axis_text_size = 10) {
  
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
  
  if (!is.null(xlab)) p <- p + xlab(xlab)
  if (!is.null(ylab)) p <- p + ylab(ylab)
  
  # Text size controls
  p <- p + theme(
    plot.title = element_text(size = title_size),
    axis.title.x = element_text(size = xlab_size),
    axis.title.y = element_text(size = ylab_size),
    axis.text.x = element_text(size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size)
  )
  
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
                       ylim = NULL,
                       main_size = 14,
                       xlab_size = 12,
                       ylab_size = 12,
                       axis_text_x_size = 10,
                       axis_text_y_size = 10) {
  
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
  
  # Text size controls
  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(size = main_size),
    axis.title.x = ggplot2::element_text(size = xlab_size),
    axis.title.y = ggplot2::element_text(size = ylab_size),
    axis.text.x = ggplot2::element_text(size = axis_text_x_size),
    axis.text.y = ggplot2::element_text(size = axis_text_y_size)
  )
  
  return(p)
}

###QTL_profile##################################################################
#plot the p-values for a single trait
#Arguments:
#           QTLlist: output of QTLgrouper
#           predictorLocs: location of predictors in the genome in absolute 
#                          cummulative positions
#           chrLen: vector containing the length of each chromosome with names
#           xlab: desired label for the x-axis
#           ylab: desired label for the y-axis
#           main: desired main-title
QTL_profile <- function(p_vec,
                        predictorLocs,
                        chrLen,
                        xlab,
                        ylab,
                        main,
                        col="black",
                        labcex=0.6,
                        chrLty=2,
                        chrCol="red"){
  plot.default(y = -log10(p_vec),
               x = predictorLocs,
               main = main,
               xlab = xlab,
               ylab = ylab,
               type = "h",
               col = col,
               xaxt = "n",
               cex.lab=labcex)
  chrLenCum <- sapply(1:length(chrLen),FUN=function(x){
    sum(chrLen[1:x])
  })
  for(i in 1:(length(chrLenCum)-1)){
    abline(v=chrLenCum[i],lty=chrLty,col=chrCol)
  }
  borders <- c(0,chrLenCum)
  for(i in 1:length(chrLen)){
    mtext(names(chrLen[i]),
          side=1,
          at=mean(borders[i:(i+1)]),
          cex=par()$cex*labcex*0.75)
  }
}
