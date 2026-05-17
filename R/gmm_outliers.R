#' Gaussian Mixture Model with Outlier Detection
#'
#' @description
#' Performs clustering and outlier detection using Gaussian Mixture Models
#' with uniform component for outliers. Model selection is performed using BIC.
#'
#' @param data Numeric matrix or data frame of observations (n x d)
#' @param K Integer, number of Gaussian components (default: 2)
#' @param n_init Integer, number of random initializations (default: 20)
#' @param max_iter Integer, maximum iterations for EM algorithm (default: 100)
#' @param tol Numeric, convergence tolerance (default: 1e-6)
#'
#' @return A list containing:
#' \item{has_uniform}{Boolean indicating if uniform component was selected}
#' \item{responsibilities}{Matrix of component responsibilities (n x (K+1))}
#' \item{map_assignments}{Cluster assignments by MAP}
#' \item{outliers}{Boolean vector indicating outliers (if has_uniform is TRUE)}
#' \item{log_likelihood_trace}{Log-likelihood values during EM iterations}
#' \item{BIC}{BIC value of selected model}
#' \item{parameters}{Estimated model parameters}
#' \item{K}{Number of Gaussian components used}
#' \item{n_iter}{Number of EM iterations performed}
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' library(mvtnorm)
#' data <- rbind(
#'   mvtnorm::rmvnorm(100, mean = c(0, 0)),
#'   mvtnorm::rmvnorm(100, mean = c(4, 4)),
#'   matrix(stats::runif(20, -5, 10), ncol = 2)  # Outliers
#' )
#'
#' # Run clustering with outlier detection
#' result <- gmm_outliers(data, K = 2, n_init = 5)
#'
#' # Check if outliers were detected
#' print(result$has_uniform)
#'
#' # View cluster assignments
#' table(result$map_assignments)
#'
#' # Plot results
#' if (interactive()) {
#'   plot(result, data = data)
#'   plot(result, type = "convergence")
#' }
gmm_outliers <- function(data, K = 2, n_init = 20, max_iter = 100, tol = 1e-6) {

  # Validate inputs
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  if (ncol(data) < 1) {
    stop("Data must have at least one dimension")
  }

  if (K < 1) {
    stop("K must be at least 1")
  }

  # Run model selection
  result <- select_best_model(data, K, n_init)

  # Format output
  output <- list(
    has_uniform = result$has_uniform_selected,
    responsibilities = result$responsibilities,
    map_assignments = result$map_assignments,
    outliers = result$outliers,
    log_likelihood_trace = result$log_likelihood_trace,
    BIC = result$BIC,
    parameters = result$params,
    K = K,
    n_iter = result$n_iter,
    data = data  # Stocker les données pour le plot
  )

  class(output) <- "gmm_outliers"

  return(output)
}

#' Plot method for gmm_outliers results
#'
#' Creates visualizations of clustering results.
#' For 2D data: scatter plot with cluster assignments and outliers.
#' For higher dimensions: first two principal components.
#'
#' @param x gmm_outliers object
#' @param type Character, type of plot: "clusters" or "convergence"
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns the input object
#' @export
#' @method plot gmm_outliers
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(200), ncol = 2)
#' result <- gmm_outliers(data, K = 2, n_init = 5)
#' plot(result)
#' plot(result, type = "convergence")
#' }
plot.gmm_outliers <- function(x, type = "clusters", ...) {

  if (type == "convergence") {
    # Plot convergence of EM algorithm
    plot_convergence(x, ...)
  } else if (type == "clusters") {
    # Plot clustering results
    plot_clusters(x, ...)
  } else {
    stop("type must be 'clusters' or 'convergence'")
  }

  invisible(x)
}

# Helper function for convergence plot
plot_convergence <- function(x, main = NULL, ...) {
  ll_trace <- x$log_likelihood_trace

  # Définir un titre par défaut si non fourni
  if (is.null(main)) {
    main <- "EM Algorithm Convergence"
  }

  plot(seq_along(ll_trace), ll_trace,
       type = "b",
       xlab = "Iteration",
       ylab = "Log-Likelihood",
       main = main,  # Utiliser le titre fourni ou le défaut
       pch = 19,
       col = "blue",
       ...)

  # Add final value
  points(length(ll_trace), tail(ll_trace, 1),
         pch = 17, col = "red", cex = 1.5)

  legend("bottomright",
         legend = c("Iterations", "Final"),
         pch = c(19, 17),
         col = c("blue", "red"))
}

# Helper function for cluster plot
plot_clusters <- function(x, data = NULL, ...) {
  # Check if we have the original data
  if (is.null(data)) {
    message("Original data not provided. Cannot plot clusters.")
    return(invisible(x))
  }

  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  d <- ncol(data)
  n <- nrow(data)

  if (d == 1) {
    # 1D plot
    plot_1d_clusters(x, data, ...)
  } else if (d == 2) {
    # 2D scatter plot
    plot_2d_clusters(x, data, ...)
  } else {
    # PCA for higher dimensions
    plot_pca_clusters(x, data, ...)
  }
}

# 1D visualization
plot_1d_clusters <- function(x, data, ...) {

  if (is.null(main)) {
    main <- "Cluster Assignments (1D)"
  }

  plot(data[, 1], rep(0, nrow(data)),
       col = x$map_assignments,
       pch = if (!is.null(x$outliers) && any(x$outliers)) {
         ifelse(x$outliers, 4, 19)
       } else {
         19
       },
       xlab = "Feature 1",
       ylab = "",
       main = main,
       yaxt = "n",
       ...)

  # Add cluster centers
  points(x$parameters$mu[, 1], rep(0, x$K),
         pch = 17, col = "red", cex = 2)

  legend("topright",
         legend = c(paste("Cluster", 1:x$K),
                    if (!is.null(x$outliers) && any(x$outliers)) "Outliers" else NULL,
                    "Centers"),
         pch = c(rep(19, x$K),
                 if (!is.null(x$outliers) && any(x$outliers)) 4 else NULL,
                 17),
         col = c(1:x$K,
                 if (!is.null(x$outliers) && any(x$outliers)) "black" else NULL,
                 "red"),
         bg = "white")
}

# 2D visualization
plot_2d_clusters <- function(x, data, ...) {

  if (is.null(main)) {
    main <- "Cluster Assignments"
  }

  # Colors for clusters
  colors <- 1:x$K

  # Plot all points
  plot(data[, 1], data[, 2],
       col = colors[x$map_assignments],
       pch = if (!is.null(x$outliers) && any(x$outliers)) {
         ifelse(x$outliers, 4, 19)
       } else {
         19
       },
       xlab = "Feature 1",
       ylab = "Feature 2",
       main = main,
       ...)

  # Add cluster centers
  points(x$parameters$mu[, 1], x$parameters$mu[, 2],
         pch = 17, col = "red", cex = 2)

  # Add legend
  legend_items <- c(paste("Cluster", 1:x$K))
  legend_pch <- rep(19, x$K)
  legend_col <- colors

  if (!is.null(x$outliers) && any(x$outliers)) {
    legend_items <- c(legend_items, "Outliers")
    legend_pch <- c(legend_pch, 4)
    legend_col <- c(legend_col, "black")
  }

  legend_items <- c(legend_items, "Centers")
  legend_pch <- c(legend_pch, 17)
  legend_col <- c(legend_col, "red")

  legend("topright",
         legend = legend_items,
         pch = legend_pch,
         col = legend_col,
         bg = "white")
}

# PCA visualization for high dimensions
plot_pca_clusters <- function(x, data, ...) {

  if (is.null(main)) {
    main <- "Cluster Assignments (PCA)"
  }

  # Perform PCA
  pca_result <- prcomp(data, scale = TRUE)

  # Plot first two principal components
  plot(pca_result$x[, 1], pca_result$x[, 2],
       col = x$map_assignments,
       pch = if (!is.null(x$outliers) && any(x$outliers)) {
         ifelse(x$outliers, 4, 19)
       } else {
         19
       },
       xlab = paste("PC1 (", round(100 * pca_result$sdev[1]^2 / sum(pca_result$sdev^2), 1), "%)", sep = ""),
       ylab = paste("PC2 (", round(100 * pca_result$sdev[2]^2 / sum(pca_result$sdev^2), 1), "%)", sep = ""),
       main = main,
       ...)

  legend("topright",
         legend = c(paste("Cluster", 1:x$K),
                    if (!is.null(x$outliers) && any(x$outliers)) "Outliers" else NULL),
         pch = c(rep(19, x$K),
                 if (!is.null(x$outliers) && any(x$outliers)) 4 else NULL),
         col = c(1:x$K,
                 if (!is.null(x$outliers) && any(x$outliers)) "black" else NULL),
         bg = "white")
}
