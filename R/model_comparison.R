#' Compare Multiple GMM Models
#'
#' Compare Gaussian Mixture Models with different numbers of components
#' using various information criteria.
#'
#' @param data Numeric matrix or data frame
#' @param K_range Integer vector of K values to compare
#' @param n_init Number of random initializations (default: 10)
#' @param max_iter Maximum EM iterations (default: 100)
#' @param criterion Criterion for comparison: "BIC", "AIC", or "ICL" (default: "BIC")
#'
#' @return A data frame with comparison results
#' @export
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(200), ncol = 2)
#' comparison <- compare_models(data, K_range = 1:5)
#' print(comparison)
#' plot(comparison)
#' }
compare_models <- function(data, K_range = 1:5, n_init = 10, max_iter = 100,
                           criterion = "BIC") {

  # Validation des entrées
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  if (min(K_range) < 1) {
    stop("K must be at least 1")
  }

  if (!criterion %in% c("BIC", "AIC", "ICL")) {
    stop("criterion must be 'BIC', 'AIC', or 'ICL'")
  }

  cat("=== COMPARING MODELS ===\n")
  cat(sprintf("K values: %s\n", paste(K_range, collapse = ", ")))
  cat(sprintf("Criterion: %s\n", criterion))
  cat(sprintf("Number of initializations: %d\n", n_init))
  cat("\n")

  # Initialiser les résultats
  results <- list()
  comparison_df <- data.frame(
    K = integer(),
    criterion_value = numeric(),
    log_likelihood = numeric(),
    has_uniform = logical(),
    n_clusters_detected = integer(),
    n_outliers = integer(),
    n_iter = integer(),
    convergence = logical(),
    stringsAsFactors = FALSE
  )

  # Évaluer chaque K
  for (K in K_range) {
    cat(sprintf("Fitting model with K = %d...\n", K))

    # Ajuster le modèle
    model <- gmm_outliers(data, K = K, n_init = n_init, max_iter = max_iter)

    # Calculer les critères
    n <- nrow(data)
    d <- ncol(data)
    log_lik <- max(model$log_likelihood_trace)

    # Nombre de paramètres
    n_params <- calculate_n_params(model$parameters, model$has_uniform)

    # Calcul des différents critères
    BIC_value <- -2 * log_lik + n_params * log(n)
    AIC_value <- -2 * log_lik + 2 * n_params

    # ICL (avec entropie)
    if (!is.null(model$responsibilities)) {
      entropy <- -sum(model$responsibilities * log(model$responsibilities + 1e-10))
      ICL_value <- BIC_value - 2 * entropy
    } else {
      ICL_value <- NA
    }

    # Sélectionner le critère demandé
    criterion_value <- switch(criterion,
                              "BIC" = BIC_value,
                              "AIC" = AIC_value,
                              "ICL" = ICL_value)

    # Nombre de clusters détectés (sans compter les outliers comme cluster)
    n_clusters <- length(unique(model$map_assignments))
    if (model$has_uniform) {
      # La composante uniforme est comptée dans map_assignments
      # On la retire du compte des clusters
      n_clusters <- n_clusters - 1
    }

    # Ajouter au dataframe
    new_row <- data.frame(
      K = K,
      criterion_value = criterion_value,
      log_likelihood = log_lik,
      has_uniform = model$has_uniform,
      n_clusters_detected = n_clusters,
      n_outliers = if (model$has_uniform) sum(model$outliers) else 0,
      n_iter = model$n_iter,
      convergence = model$n_iter < max_iter,
      BIC = BIC_value,
      AIC = AIC_value,
      ICL = ICL_value
    )

    comparison_df <- rbind(comparison_df, new_row)
    results[[as.character(K)]] <- model

    cat(sprintf("  %s: %.2f, LogLik: %.2f, Uniform: %s\n",
                criterion, criterion_value, log_lik, model$has_uniform))
  }

  # Trouver le meilleur modèle selon le critère
  best_idx <- which.min(comparison_df$criterion_value)
  best_K <- comparison_df$K[best_idx]

  cat("\n=== BEST MODEL ===\n")
  cat(sprintf("Best K: %d (according to %s)\n", best_K, criterion))
  cat(sprintf("Criterion value: %.2f\n", comparison_df$criterion_value[best_idx]))
  cat(sprintf("Log-likelihood: %.2f\n", comparison_df$log_likelihood[best_idx]))

  # Ajouter des attributs
  attr(comparison_df, "best_K") <- best_K
  attr(comparison_df, "criterion") <- criterion
  attr(comparison_df, "models") <- results
  attr(comparison_df, "data") <- data

  class(comparison_df) <- c("gmm_comparison", "data.frame")

  return(comparison_df)
}

#' Calculate number of parameters
#' @noRd
calculate_n_params <- function(params, has_uniform) {
  K <- length(params$pi)
  d <- ncol(params$mu)

  # Paramètres des gaussiennes
  n_params <- K * d + K * d * (d + 1) / 2 + (K - 1)

  # Ajouter la composante uniforme si présente
  if (has_uniform) {
    n_params <- n_params + 2 * d + 1
  }

  return(n_params)
}

#' Plot method for model comparison
#' @param x gmm_comparison object
#' @param ... Additional arguments passed to plot
#' @export
#' @method plot gmm_comparison
plot.gmm_comparison <- function(x, ...) {
  criterion <- attr(x, "criterion")
  best_K <- attr(x, "best_K")

  par(mfrow = c(2, 2))

  # 1. Critère principal
  plot(x$K, x$criterion_value, type = "b", pch = 19,
       xlab = "Number of components (K)", ylab = criterion,
       main = paste(criterion, "vs K"),
       col = ifelse(x$K == best_K, "red", "blue"))
  points(best_K, min(x$criterion_value), pch = 17, col = "red", cex = 2)

  # 2. Log-vraisemblance
  plot(x$K, x$log_likelihood, type = "b", pch = 19,
       xlab = "K", ylab = "Log-Likelihood",
       main = "Log-Likelihood vs K",
       col = ifelse(x$K == best_K, "red", "blue"))

  # 3. Tous les critères
  plot(x$K, x$BIC, type = "b", pch = 19, col = "blue",
       xlab = "K", ylab = "Criterion Value",
       main = "All Criteria", ylim = range(c(x$BIC, x$AIC, na.omit(x$ICL))))
  lines(x$K, x$AIC, type = "b", pch = 18, col = "green")
  if (!all(is.na(x$ICL))) {
    lines(x$K, x$ICL, type = "b", pch = 17, col = "purple")
    legend("topright", legend = c("BIC", "AIC", "ICL"),
           col = c("blue", "green", "purple"), pch = c(19, 18, 17), lty = 1)
  } else {
    legend("topright", legend = c("BIC", "AIC"),
           col = c("blue", "green"), pch = c(19, 18), lty = 1)
  }

  # 4. Nombre d'outliers
  if (any(x$n_outliers > 0)) {
    barplot(x$n_outliers, names.arg = x$K,
            xlab = "K", ylab = "Number of outliers",
            main = "Outliers detected",
            col = ifelse(x$K == best_K, "red", "lightblue"))
  } else {
    # Sinon, montrer la convergence
    plot(x$K, x$n_iter, type = "b", pch = 19,
         xlab = "K", ylab = "EM iterations",
         main = "Convergence speed",
         col = ifelse(x$K == best_K, "red", "blue"))
  }

  par(mfrow = c(1, 1))

  invisible(x)
}

#' Print method for model comparison
#' @param x gmm_comparison object
#' @param ... Additional arguments
#' @export
#' @method print gmm_comparison
print.gmm_comparison <- function(x, ...) {
  cat("=== GMM MODEL COMPARISON ===\n\n")

  # Tableau principal
  print_df <- x[, c("K", "criterion_value", "log_likelihood",
                    "has_uniform", "n_clusters_detected", "n_outliers")]
  colnames(print_df) <- c("K", attr(x, "criterion"), "LogLik",
                          "Uniform", "Clusters", "Outliers")

  print(print_df, row.names = FALSE)

  cat(sprintf("\nBest model: K = %d (lowest %s)\n",
              attr(x, "best_K"), attr(x, "criterion")))

  invisible(x)
}
