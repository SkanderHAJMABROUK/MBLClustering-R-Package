#' Calculate Validation Metrics for GMM
#'
#' Compute various metrics to evaluate the quality of
#' Gaussian Mixture Model clustering.
#'
#' @param x gmm_outliers object
#' @param true_labels Optional vector of true labels for supervised metrics
#' @param data Original data (optional)
#'
#' @return List of validation metrics
#' @export
#'
#' @examples
#' \dontrun{
#' result <- gmm_outliers(data, K = 3)
#' metrics <- calculate_metrics(result)
#' print(metrics)
#' }
calculate_metrics <- function(x, true_labels = NULL, data = NULL) {

  metrics <- list()

  # 1. Basic model metrics
  metrics$basic <- calculate_basic_metrics(x)

  # 2. Clustering quality metrics
  metrics$clustering <- calculate_clustering_metrics(x, data)

  # 3. Supervised metrics (if true_labels provided)
  if (!is.null(true_labels)) {
    metrics$supervised <- calculate_supervised_metrics(x, true_labels)
  }

  # 4. Component-specific metrics
  metrics$components <- calculate_component_metrics(x)

  # 5. Overall score (weighted combination)
  metrics$overall_score <- calculate_overall_score(metrics)

  class(metrics) <- "gmm_metrics"

  return(metrics)
}

#' Calculate basic model metrics
#' @noRd
calculate_basic_metrics <- function(x) {
  list(
    log_likelihood = max(x$log_likelihood_trace),
    BIC = x$BIC,
    has_uniform = x$has_uniform,
    K = x$K,
    n_iterations = x$n_iter,
    convergence_rate = if (x$n_iter < 100) 1.0 else 0.5  # simplifie
  )
}

#' Calculate clustering quality metrics
#' @noRd
calculate_clustering_metrics <- function(x, data = NULL) {
  metrics <- list()

  # 1. Classification uncertainty
  if (!is.null(x$responsibilities)) {
    uncertainties <- 1 - apply(x$responsibilities, 1, max)
    metrics$mean_uncertainty <- mean(uncertainties)
    metrics$max_uncertainty <- max(uncertainties)
    metrics$entropy <- -sum(x$responsibilities * log(x$responsibilities + 1e-10))
  }

  # 2. Separation between components (si donnees disponibles)
  if (!is.null(data) && !is.null(x$parameters)) {
    separation <- calculate_component_separation(x$parameters, data)
    metrics$component_separation <- separation
  }

  # 3. Outlier metrics
  if (x$has_uniform && !is.null(x$outliers)) {
    metrics$outlier_proportion <- mean(x$outliers)
    metrics$n_outliers <- sum(x$outliers)
  }

  # 4. Cluster balance
  if (!is.null(x$map_assignments)) {
    if (x$has_uniform && !is.null(x$outliers)) {
      # Exclure les outliers
      cluster_counts <- table(x$map_assignments[!x$outliers])
    } else {
      cluster_counts <- table(x$map_assignments)
    }

    if (length(cluster_counts) > 0) {
      metrics$cluster_balance <- sd(cluster_counts) / mean(cluster_counts)
      metrics$min_cluster_size <- min(cluster_counts)
      metrics$max_cluster_size <- max(cluster_counts)
    }
  }

  return(metrics)
}

#' Calculate component separation
#' @noRd
calculate_component_separation <- function(params, data) {
  K <- length(params$pi)
  d <- ncol(params$mu)

  if (K < 2) return(NA)

  # Calculer les distances entre les centres
  center_distances <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        center_distances[i, j] <- sqrt(sum((params$mu[i, ] - params$mu[j, ])^2))
      }
    }
  }

  # Normaliser par la dispersion des donnees
  data_range <- apply(data, 2, function(col) diff(range(col)))
  avg_data_range <- mean(data_range)

  list(
    min_separation = min(center_distances[center_distances > 0]),
    mean_separation = mean(center_distances[center_distances > 0]),
    normalized_separation = mean(center_distances[center_distances > 0]) / avg_data_range
  )
}

#' Calculate supervised metrics
#' @noRd
calculate_supervised_metrics <- function(x, true_labels) {
  metrics <- list()

  # Predictions
  pred_labels <- x$map_assignments

  # Si outliers, les traiter comme classe supplementaire
  if (x$has_uniform && !is.null(x$outliers)) {
    pred_labels <- ifelse(x$outliers, max(pred_labels) + 1, pred_labels)
  }

  # 1. Accuracy
  metrics$accuracy <- sum(pred_labels == true_labels) / length(true_labels)

  # 2. Adjusted Rand Index (si package mclust disponible)
  if (requireNamespace("mclust", quietly = TRUE)) {
    metrics$adjusted_rand_index <- mclust::adjustedRandIndex(pred_labels, true_labels)
  }

  # 3. Confusion matrix
  conf_matrix <- table(Predicted = pred_labels, True = true_labels)
  metrics$confusion_matrix <- conf_matrix

  # 4. Precision, Recall, F1 pour chaque classe
  n_classes <- length(unique(true_labels))
  precision <- numeric(n_classes)
  recall <- numeric(n_classes)
  f1 <- numeric(n_classes)

  for (i in 1:n_classes) {
    tp <- sum(pred_labels == i & true_labels == i)
    fp <- sum(pred_labels == i & true_labels != i)
    fn <- sum(pred_labels != i & true_labels == i)

    precision[i] <- if (tp + fp > 0) tp / (tp + fp) else 0
    recall[i] <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1[i] <- if (precision[i] + recall[i] > 0)
      2 * precision[i] * recall[i] / (precision[i] + recall[i]) else 0
  }

  metrics$precision <- precision
  metrics$recall <- recall
  metrics$f1_score <- f1
  metrics$macro_f1 <- mean(f1)

  return(metrics)
}

#' Calculate component-specific metrics
#' @noRd
calculate_component_metrics <- function(x) {
  K <- x$K
  metrics <- list()

  # Pour chaque composante
  for (k in 1:K) {
    comp_metrics <- list()

    # Taille effective
    if (!is.null(x$responsibilities)) {
      comp_metrics$effective_size <- sum(x$responsibilities[, k])
    }

    # Regularite de la covariance
    sigma <- x$parameters$sigma[, , k]
    eigen_vals <- eigen(sigma, symmetric = TRUE)$values
    comp_metrics$condition_number <- max(eigen_vals) / min(eigen_vals)
    comp_metrics$sphericity <- (prod(eigen_vals)^(1/length(eigen_vals))) / mean(eigen_vals)

    metrics[[paste0("component_", k)]] <- comp_metrics
  }

  return(metrics)
}

#' Calculate overall score
#' @noRd
calculate_overall_score <- function(metrics) {
  score <- 0
  weights <- list(
    log_likelihood = 0.2,
    BIC = 0.2,
    uncertainty = 0.15,
    separation = 0.15,
    balance = 0.1,
    convergence = 0.1,
    other = 0.1
  )

  # Score base sur la log-vraisemblance (normalisee)
  if (!is.null(metrics$basic$log_likelihood)) {
    # Juste un exemple simple
    score <- score + weights$log_likelihood *
      (1 - exp(-abs(metrics$basic$log_likelihood) / 1000))
  }

  # Score base sur l'incertitude (plus bas = mieux)
  if (!is.null(metrics$clustering$mean_uncertainty)) {
    score <- score + weights$uncertainty * (1 - metrics$clustering$mean_uncertainty)
  }

  # Score base sur la separation
  if (!is.null(metrics$clustering$component_separation)) {
    sep <- metrics$clustering$component_separation$normalized_separation
    if (!is.na(sep)) {
      score <- score + weights$separation * min(sep, 1)
    }
  }

  # Score base sur l'equilibre des clusters
  if (!is.null(metrics$clustering$cluster_balance)) {
    score <- score + weights$balance * (1 - metrics$clustering$cluster_balance)
  }

  return(min(1, max(0, score)))  # Normaliser entre 0 et 1
}

#' Print method for metrics
#' @param x gmm_metrics object
#' @param digits Number of digits to display
#' @param ... Additional arguments
#' @export
#' @method print gmm_metrics
print.gmm_metrics <- function(x, digits = 3, ...) {
  cat("=== GMM VALIDATION METRICS ===\n\n")

  # Basic metrics
  cat("BASIC METRICS:\n")
  cat(sprintf("  Log-likelihood: %.2f\n", x$basic$log_likelihood))
  cat(sprintf("  BIC: %.2f\n", x$basic$BIC))
  cat(sprintf("  K: %d\n", x$basic$K))
  cat(sprintf("  Uniform component: %s\n", x$basic$has_uniform))
  cat(sprintf("  EM iterations: %d\n", x$basic$n_iterations))
  cat("\n")

  # Clustering metrics
  cat("CLUSTERING QUALITY:\n")
  if (!is.null(x$clustering$mean_uncertainty)) {
    cat(sprintf("  Mean uncertainty: %.3f\n", x$clustering$mean_uncertainty))
  }
  if (!is.null(x$clustering$entropy)) {
    cat(sprintf("  Classification entropy: %.3f\n", x$clustering$entropy))
  }
  if (!is.null(x$clustering$component_separation)) {
    cat(sprintf("  Component separation: %.3f\n",
                x$clustering$component_separation$normalized_separation))
  }
  if (!is.null(x$clustering$outlier_proportion)) {
    cat(sprintf("  Outlier proportion: %.1f%%\n",
                x$clustering$outlier_proportion * 100))
  }
  if (!is.null(x$clustering$cluster_balance)) {
    cat(sprintf("  Cluster balance (CV): %.3f\n", x$clustering$cluster_balance))
  }
  cat("\n")

  # Supervised metrics
  if (!is.null(x$supervised)) {
    cat("SUPERVISED METRICS:\n")
    cat(sprintf("  Accuracy: %.3f\n", x$supervised$accuracy))
    if (!is.null(x$supervised$adjusted_rand_index)) {
      cat(sprintf("  Adjusted Rand Index: %.3f\n", x$supervised$adjusted_rand_index))
    }
    if (!is.null(x$supervised$macro_f1)) {
      cat(sprintf("  Macro F1-score: %.3f\n", x$supervised$macro_f1))
    }
    cat("\n")
  }

  # Overall score
  cat(sprintf("OVERALL SCORE: %.1f/1.0\n", x$overall_score))

  invisible(x)
}

#' Summary method for metrics
#' @param object gmm_metrics object
#' @param ... Additional arguments
#' @export
#' @method summary gmm_metrics
summary.gmm_metrics <- function(object, ...) {
  cat("=== GMM METRICS SUMMARY ===\n\n")

  # Tableau recapitulatif
  summary_df <- data.frame(
    Metric = c("Log-Likelihood", "BIC", "Uncertainty",
               "Separation", "Balance", "Overall Score"),
    Value = c(
      sprintf("%.2f", object$basic$log_likelihood),
      sprintf("%.2f", object$basic$BIC),
      if (!is.null(object$clustering$mean_uncertainty))
        sprintf("%.3f", object$clustering$mean_uncertainty) else "N/A",
      if (!is.null(object$clustering$component_separation))
        sprintf("%.3f", object$clustering$component_separation$normalized_separation)
      else "N/A",
      if (!is.null(object$clustering$cluster_balance))
        sprintf("%.3f", object$clustering$cluster_balance) else "N/A",
      sprintf("%.3f", object$overall_score)
    ),
    Interpretation = c(
      ifelse(object$basic$log_likelihood > -500, "Good", "Poor"),
      ifelse(object$basic$BIC < object$basic$log_likelihood * (-2), "Good", "OK"),
      ifelse(!is.null(object$clustering$mean_uncertainty) &&
               object$clustering$mean_uncertainty < 0.1, "Good", "OK"),
      ifelse(!is.null(object$clustering$component_separation) &&
               object$clustering$component_separation$normalized_separation > 0.5,
             "Good", "OK"),
      ifelse(!is.null(object$clustering$cluster_balance) &&
               object$clustering$cluster_balance < 0.5, "Good", "OK"),
      ifelse(object$overall_score > 0.7, "Excellent",
             ifelse(object$overall_score > 0.5, "Good", "Needs improvement"))
    )
  )

  print(summary_df, row.names = FALSE)
  cat("\n")

  # Recommandations
  cat("RECOMMENDATIONS:\n")
  if (object$clustering$mean_uncertainty > 0.2) {
    cat("  - High uncertainty: Consider increasing K or checking data quality\n")
  }
  if (object$clustering$component_separation$normalized_separation < 0.3) {
    cat("  - Low separation: Components may be overlapping\n")
  }
  if (object$basic$n_iterations >= 100) {
    cat("  - Slow convergence: Consider increasing max_iter or changing initialization\n")
  }

  invisible(object)
}
