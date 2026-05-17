#' Diagnostic Plots for GMM Results
#'
#' Generate comprehensive diagnostic plots to evaluate
#' Gaussian Mixture Model clustering quality.
#'
#' @param x gmm_outliers object
#' @param data Original data (optional, stored in x by default)
#' @param type Type of diagnostics: "basic", "advanced", or "all"
#' @param ... Additional arguments passed to plot functions
#'
#' @return Invisibly returns the input object
#' @export
#'
#' @examples
#' \dontrun{
#' result <- gmm_outliers(data, K = 3)
#' plot_diagnostics(result)
#' }
plot_diagnostics <- function(x, data = NULL, type = "basic", ...) {

  # Récupérer les données
  if (is.null(data)) {
    if (!is.null(x$data)) {
      data <- x$data
    } else {
      stop("Please provide data or run gmm_outliers with store_data = TRUE")
    }
  }

  # Validation
  if (!inherits(x, "gmm_outliers")) {
    stop("x must be a gmm_outliers object")
  }

  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  cat("=== GMM DIAGNOSTICS ===\n")
  cat(sprintf("K: %d, Has uniform: %s\n", x$K, x$has_uniform))
  cat(sprintf("BIC: %.2f, LogLik: %.2f\n", x$BIC, max(x$log_likelihood_trace)))

  # Sélectionner les plots selon le type
  if (type == "basic") {
    par(mfrow = c(2, 2))
    plot_convergence_diagnostic(x, ...)
    plot_responsibilities_diagnostic(x, ...)
    plot_mahalanobis_diagnostic(x, data, ...)
    plot_cluster_sizes(x, ...)
  } else if (type == "advanced") {
    par(mfrow = c(3, 2))
    plot_convergence_diagnostic(x, ...)
    plot_responsibilities_diagnostic(x, ...)
    plot_mahalanobis_diagnostic(x, data, ...)
    plot_cluster_sizes(x, ...)
    plot_covariance_diagnostic(x, ...)
    plot_bic_sensitivity(x, data, ...)
  } else if (type == "all") {
    # Créer un PDF avec tous les diagnostics
    pdf_file <- tempfile("gmm_diagnostics_", fileext = ".pdf")
    pdf(pdf_file, width = 11, height = 8.5)

    layout(matrix(c(1,2,3,4,5,6,7,7), 4, 2, byrow = TRUE))

    plot_convergence_diagnostic(x, main = "1. EM Convergence", ...)
    plot_responsibilities_diagnostic(x, main = "2. Responsibilities", ...)
    plot_mahalanobis_diagnostic(x, data, main = "3. Mahalanobis Distances", ...)
    plot_cluster_sizes(x, main = "4. Cluster Sizes", ...)
    plot_covariance_diagnostic(x, main = "5. Covariance Structure", ...)
    plot_bic_sensitivity(x, data, main = "6. Model Sensitivity", ...)

    # Résumé textuel
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
         xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.5, generate_diagnostic_summary(x), cex = 0.8, family = "mono")
    title("7. Diagnostic Summary")

    dev.off()
    cat(sprintf("\nDiagnostic report saved to: %s\n", pdf_file))

    # Ouvrir le fichier si possible
    if (interactive() && .Platform$OS.type == "windows") {
      shell.exec(pdf_file)
    }

    return(invisible(x))
  }

  par(mfrow = c(1, 1))
  invisible(x)
}

#' Convergence diagnostic plot
#' @noRd
plot_convergence_diagnostic <- function(x, ...) {
  ll_trace <- x$log_likelihood_trace

  # Plot de convergence
  plot(seq_along(ll_trace), ll_trace, type = "b",
       xlab = "Iteration", ylab = "Log-Likelihood",
       main = "EM Convergence",
       pch = 19, col = "blue", ...)

  # Ajouter des informations
  final_ll <- tail(ll_trace, 1)
  initial_ll <- ll_trace[1]
  improvement <- final_ll - initial_ll

  legend("bottomright",
         legend = c(sprintf("Iterations: %d", length(ll_trace)),
                    sprintf("Improvement: %.1f", improvement)),
         bg = "white")
}

#' Responsibilities diagnostic plot
#' @noRd
plot_responsibilities_diagnostic <- function(x, ...) {
  if (is.null(x$responsibilities)) {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Responsibilities\n(Not available)")
    return()
  }

  # Calculer l'incertitude
  uncertainties <- 1 - apply(x$responsibilities, 1, max)

  # Histogramme
  hist(uncertainties, breaks = 30, col = "lightblue",
       xlab = "1 - Max Probability",
       ylab = "Frequency",
       main = "Classification Uncertainty", ...)

  abline(v = mean(uncertainties), col = "red", lty = 2, lwd = 2)

  legend("topright",
         legend = c(sprintf("Mean: %.3f", mean(uncertainties)),
                    sprintf("SD: %.3f", sd(uncertainties))),
         bg = "white")
}

#' Mahalanobis distances diagnostic
#' @noRd
plot_mahalanobis_diagnostic <- function(x, data, ...) {
  # Calculer les distances de Mahalanobis
  distances <- calculate_mahalanobis_distances(x, data)

  # QQ-plot
  qqnorm(distances, main = "Mahalanobis Distances QQ-plot", ...)
  qqline(distances, col = "red")

  # Ajouter le test de normalité
  if (length(distances) > 3) {
    shapiro_test <- shapiro.test(distances)
    legend("topleft",
           legend = c(sprintf("Shapiro p-value: %.3f", shapiro_test$p.value)),
           bg = "white")
  }
}

#' Calculate Mahalanobis distances
#' @noRd
calculate_mahalanobis_distances <- function(x, data) {
  n <- nrow(data)
  distances <- numeric(n)

  for (i in 1:n) {
    # Trouver le cluster assigné
    cluster <- x$map_assignments[i]

    if (x$has_uniform && x$outliers[i]) {
      # Pour les outliers, distance NA
      distances[i] <- NA
    } else {
      # Calculer la distance de Mahalanobis
      mu <- x$parameters$mu[cluster, ]
      sigma <- x$parameters$sigma[, , cluster]

      # Ajouter une régularisation
      sigma <- sigma + diag(ncol(data)) * 1e-6

      diff <- data[i, ] - mu
      distances[i] <- sqrt(t(diff) %*% solve(sigma) %*% diff)
    }
  }

  # Retirer les NA (outliers)
  return(distances[!is.na(distances)])
}

#' Cluster sizes plot
#' @noRd
plot_cluster_sizes <- function(x, ...) {
  # Compter les tailles de clusters
  if (x$has_uniform && !is.null(x$outliers)) {
    # Séparer les clusters des outliers
    cluster_assignments <- x$map_assignments[!x$outliers]
    n_outliers <- sum(x$outliers)

    cluster_counts <- table(cluster_assignments)
    all_counts <- c(cluster_counts, Outliers = n_outliers)
  } else {
    all_counts <- table(x$map_assignments)
  }

  # Bar plot
  bar_colors <- if (x$has_uniform && !is.null(x$outliers)) {
    c(rainbow(length(cluster_counts)), "gray50")
  } else {
    rainbow(length(all_counts))
  }

  barplot(all_counts, col = bar_colors,
          xlab = "Cluster", ylab = "Number of points",
          main = "Cluster Sizes", ...)

  # Ajouter les proportions
  proportions <- round(all_counts / sum(all_counts) * 100, 1)
  text(1:length(all_counts), all_counts,
       labels = paste0(proportions, "%"),
       pos = 3, cex = 0.8)
}

#' Covariance structure diagnostic
#' @noRd
plot_covariance_diagnostic <- function(x, ...) {
  K <- x$K
  d <- ncol(x$parameters$mu)

  if (d == 1) {
    # Cas 1D : montrer les variances
    variances <- numeric(K)
    for (k in 1:K) {
      variances[k] <- x$parameters$sigma[, , k]
    }

    plot(1:K, variances, type = "b", pch = 19,
         xlab = "Component", ylab = "Variance",
         main = "Component Variances",
         ylim = c(0, max(variances) * 1.1), ...)
  } else if (d == 2) {
    # Cas 2D : montrer les ellipses de covariance
    plot(0, 0, type = "n",
         xlim = range(x$parameters$mu[, 1]) + c(-1, 1),
         ylim = range(x$parameters$mu[, 2]) + c(-1, 1),
         xlab = "Dimension 1", ylab = "Dimension 2",
         main = "Covariance Ellipses (95%)", ...)

    colors <- rainbow(K)
    for (k in 1:K) {
      # Ajouter l'ellipse
      draw_covariance_ellipse(x$parameters$mu[k, ],
                              x$parameters$sigma[, , k],
                              level = 0.95, col = colors[k])

      # Ajouter le centre
      points(x$parameters$mu[k, 1], x$parameters$mu[k, 2],
             pch = 19, col = colors[k], cex = 1.5)
    }

    legend("topright", legend = paste("Component", 1:K),
           col = colors, pch = 19, bg = "white")
  } else {
    # Cas multi-D : montrer les valeurs propres
    eigen_values <- list()
    max_eigen <- 0

    for (k in 1:K) {
      eigen_vals <- eigen(x$parameters$sigma[, , k], symmetric = TRUE)$values
      eigen_values[[k]] <- eigen_vals
      max_eigen <- max(max_eigen, max(eigen_vals))
    }

    plot(1:d, eigen_values[[1]], type = "b", pch = 19,
         xlab = "Eigenvalue Index", ylab = "Eigenvalue",
         main = "Covariance Eigenvalues",
         ylim = c(0, max_eigen * 1.1), col = 1, ...)

    for (k in 2:K) {
      lines(1:d, eigen_values[[k]], type = "b", pch = 19, col = k)
    }

    legend("topright", legend = paste("Component", 1:K),
           col = 1:K, pch = 19, bg = "white")
  }
}

#' Draw covariance ellipse
#' @noRd
draw_covariance_ellipse <- function(center, sigma, level = 0.95, ...) {
  angles <- seq(0, 2 * pi, length.out = 200)
  circle <- cbind(cos(angles), sin(angles))

  # Décomposition de la covariance
  eig <- eigen(sigma)
  axes <- eig$vectors
  radii <- sqrt(eig$values * qchisq(level, 2))

  # Transformer le cercle en ellipse
  ellipse <- t(center + t(circle %*% diag(radii) %*% t(axes)))

  lines(ellipse, ...)
}

#' BIC sensitivity analysis
#' @noRd
plot_bic_sensitivity <- function(x, data, n_repeats = 5, ...) {
  K_values <- max(1, x$K - 2):min(x$K + 2, 8)
  bic_matrix <- matrix(NA, n_repeats, length(K_values))
  colnames(bic_matrix) <- K_values

  cat("\nRunning sensitivity analysis...\n")
  for (i in 1:n_repeats) {
    for (j in seq_along(K_values)) {
      K <- K_values[j]
      model <- gmm_outliers(data, K = K, n_init = 5, max_iter = 50)
      bic_matrix[i, j] <- model$BIC
    }
  }

  # Boxplot des BIC
  boxplot(bic_matrix, xlab = "K", ylab = "BIC",
          main = "BIC Sensitivity (multiple runs)", ...)

  # Ajouter la ligne pour le modèle original
  points(which(K_values == x$K), x$BIC, pch = 17, col = "red", cex = 2)

  legend("topright", legend = "Original model",
         pch = 17, col = "red", bg = "white")
}

#' Generate diagnostic summary text
#' @noRd
generate_diagnostic_summary <- function(x) {
  summary_text <- c(
    "=== GMM DIAGNOSTIC SUMMARY ===",
    "",
    sprintf("Model configuration:"),
    sprintf("  K (Gaussian components): %d", x$K),
    sprintf("  Uniform component: %s", x$has_uniform),
    "",
    sprintf("Convergence:"),
    sprintf("  EM iterations: %d", x$n_iter),
    sprintf("  Log-likelihood: %.2f", max(x$log_likelihood_trace)),
    sprintf("  BIC: %.2f", x$BIC),
    "",
    sprintf("Clustering quality:"),
    sprintf("  Classification uncertainty: %.3f",
            mean(1 - apply(x$responsibilities, 1, max)))
  )

  if (x$has_uniform) {
    summary_text <- c(summary_text,
                      sprintf("  Outliers detected: %d (%.1f%%)",
                              sum(x$outliers),
                              mean(x$outliers) * 100))
  }

  return(paste(summary_text, collapse = "\n"))
}
