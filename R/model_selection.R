#' Calculate BIC for model selection
#'
#' @param log_likelihood Log-likelihood value
#' @param n Number of observations
#' @param params Model parameters
#' @param has_uniform Boolean for uniform component
#' @return BIC value
calculate_BIC <- function(log_likelihood, n, params, has_uniform) {
  K <- length(params$pi)
  d <- ncol(params$mu)

  # Vérifier la validité de la log-vraisemblance
  if (!is.finite(log_likelihood) || log_likelihood > 0) {
    # La log-vraisemblance devrait être négative
    warning("Invalid log-likelihood: ", log_likelihood)
    return(Inf)
  }

  # Number of parameters CORRIGÉ :
  # 1. Moyennes : K * d
  # 2. Matrices de covariance : K * d * (d + 1) / 2
  # 3. Proportions de mélange : K - 1 (car somme = 1)

  n_params <- K * d + K * d * (d + 1) / 2 + (K - 1)

  if (has_uniform) {
    # Composante uniforme : 2*d pour les bornes min/max
    # + 1 paramètre pour la proportion de la composante uniforme
    n_params <- n_params + 2 * d + 1
  }

  # BIC = -2 * log(L) + p * log(n)
  # Mais attention : votre log_likelihood est déjà la LOG-vraisemblance
  BIC <- -2 * log_likelihood + n_params * log(n)

  return(BIC)
}

#' Model selection using BIC
#'
#' @param data Data matrix
#' @param K Number of Gaussian components
#' @param n_init Number of random initializations
#' @return Best model according to BIC
select_best_model <- function(data, K, n_init = 20) {
  best_models <- list()
  BIC_values <- numeric(2)
  log_lik_values <- numeric(2)  # Pour stocker les log-vraisemblances

  # Model without uniform component
  cat("Fitting model without uniform component...\n")
  best_no_uniform <- run_multiple_initializations(data, K, n_init, has_uniform = FALSE)

  log_lik_no_uniform <- max(best_no_uniform$log_likelihood_trace)
  BIC_values[1] <- calculate_BIC(
    log_lik_no_uniform,
    nrow(data),
    best_no_uniform$params,
    has_uniform = FALSE
  )
  log_lik_values[1] <- log_lik_no_uniform
  best_models[[1]] <- best_no_uniform

  # Model with uniform component
  cat("Fitting model with uniform component...\n")
  best_with_uniform <- run_multiple_initializations(data, K, n_init, has_uniform = TRUE)

  log_lik_with_uniform <- max(best_with_uniform$log_likelihood_trace)
  BIC_values[2] <- calculate_BIC(
    log_lik_with_uniform,
    nrow(data),
    best_with_uniform$params,
    has_uniform = TRUE
  )
  log_lik_values[2] <- log_lik_with_uniform
  best_models[[2]] <- best_with_uniform

  # Afficher les valeurs pour débogage
  cat("\n=== DEBUG BIC CALCULATION ===\n")
  cat("Model without uniform:\n")
  cat("  Log-likelihood:", log_lik_values[1], "\n")
  cat("  BIC:", BIC_values[1], "\n")
  cat("Model with uniform:\n")
  cat("  Log-likelihood:", log_lik_values[2], "\n")
  cat("  BIC:", BIC_values[2], "\n")

  # Select best model based on BIC
  best_idx <- which.min(BIC_values)
  best_model <- best_models[[best_idx]]
  best_model$BIC <- BIC_values[best_idx]
  best_model$has_uniform_selected <- (best_idx == 2)

  return(best_model)
}

#' Run multiple initializations with error handling
#'
#' @param data Data matrix
#' @param K Number of Gaussian components
#' @param n_init Number of initializations
#' @param has_uniform Boolean indicating if uniform component is included
#' @return Best model from initializations
#' @keywords internal
run_multiple_initializations <- function(data, K, n_init, has_uniform) {
  best_log_likelihood <- -Inf
  best_result <- NULL
  successful_inits <- 0

  # Pour K=1, utiliser moins d'initialisations mais plus robustes
  if (K == 1) {
    n_init <- max(n_init, 10)  # Réduire à 10 pour K=1
    cat("Running initializations for K=1...\n")
  }

  for (init in 1:n_init) {
    if (init %% 5 == 0 || K == 1) {
      cat(sprintf("Initialization %d/%d\n", init, n_init))
    }

    result <- tryCatch({
      # Appeler EM_algorithm avec gestion d'erreurs
      res <- EM_algorithm(data, K, has_uniform = has_uniform)

      # Vérifier que le résultat est valide
      if (is.null(res) || !is.list(res) || length(res$log_likelihood_trace) == 0) {
        stop("Invalid EM result")
      }

      res
    }, error = function(e) {
      if (K == 1 || init %% 10 == 0) {
        cat(sprintf("Initialization %d failed: %s\n", init, e$message))
      }
      NULL
    })

    if (!is.null(result)) {
      successful_inits <- successful_inits + 1
      current_log_likelihood <- max(result$log_likelihood_trace)
      if (current_log_likelihood > best_log_likelihood) {
        best_log_likelihood <- current_log_likelihood
        best_result <- result
      }
    }
  }

  cat(sprintf("Successful initializations: %d/%d\n", successful_inits, n_init))

  if (successful_inits == 0) {
    # Dernière tentative avec paramètres très simples
    warning("All standard initializations failed. Using fallback initialization.")

    # Initialisation fallback TRÈS simple
    d <- ncol(data)

    # Paramètres de base
    params <- list(
      pi = 1,
      mu = matrix(colMeans(data, na.rm = TRUE), nrow = 1),
      sigma = array(diag(d) * 0.5, dim = c(d, d, 1)),
      uniform_min = apply(data, 2, min, na.rm = TRUE),
      uniform_max = apply(data, 2, max, na.rm = TRUE)
    )

    if (has_uniform) {
      params$eta <- 0.05
      params$pi <- 0.95
    }

    # Essayer EM avec ces paramètres
    best_result <- tryCatch({
      EM_algorithm_simple(data, K, params, has_uniform)
    }, error = function(e) {
      # Si même l'initialisation simple échoue, retourner un résultat minimal
      cat("Fallback initialization also failed. Returning minimal result.\n")
      return_minimal_result(data, K, has_uniform)
    })
  }

  return(best_result)
}

# Version simplifiée de EM pour l'initialisation fallback
EM_algorithm_simple <- function(data, K, init_params, has_uniform = FALSE, max_iter = 50, tol = 1e-4) {
  n <- nrow(data)
  d <- ncol(data)

  params <- init_params
  log_likelihood_trace <- numeric(max_iter)

  for (iter in 1:max_iter) {
    # E-step très simple
    if (K == 1 && !has_uniform) {
      # Cas simple : une seule composante gaussienne
      resp <- matrix(1, nrow = n, ncol = 1)
    } else {
      # Calculer les responsabilités
      log_resp <- E_step(data, params, has_uniform)
      resp <- exp(log_resp)
    }

    # M-step simplifié
    if (K == 1 && !has_uniform) {
      # Mise à jour simple pour une seule gaussienne
      params$mu <- matrix(colMeans(data), nrow = 1)
      centered <- sweep(data, 2, params$mu[1, ], "-")
      params$sigma[, , 1] <- (t(centered) %*% centered) / n
      params$sigma[, , 1] <- params$sigma[, , 1] + diag(d) * 1e-6
    }

    # Calculer la log-vraisemblance
    log_likelihood_trace[iter] <- calculate_log_likelihood(data, params, has_uniform)

    # Convergence simple
    if (iter > 1 && abs(log_likelihood_trace[iter] - log_likelihood_trace[iter-1]) < tol) {
      break
    }
  }

  # Calculer les affectations finales
  if (K == 1 && !has_uniform) {
    map_assignments <- rep(1, n)
    outliers <- rep(FALSE, n)
  } else {
    log_resp <- E_step(data, params, has_uniform)
    resp <- exp(log_resp)
    map_assignments <- apply(resp, 1, which.max)

    outliers <- NULL
    if (has_uniform) {
      outliers <- map_assignments == (K + 1)
    }
  }

  return(list(
    params = params,
    responsibilities = if (K == 1 && !has_uniform) matrix(1, n, 1) else resp,
    map_assignments = map_assignments,
    outliers = outliers,
    log_likelihood_trace = log_likelihood_trace[1:iter],
    n_iter = iter,
    has_uniform = has_uniform
  ))
}

# Résultat minimal si tout échoue
return_minimal_result <- function(data, K, has_uniform) {
  n <- nrow(data)
  d <- ncol(data)

  # Paramètres minimaux
  params <- list(
    pi = 1,
    mu = matrix(colMeans(data, na.rm = TRUE), nrow = 1),
    sigma = array(diag(d), dim = c(d, d, 1)),
    uniform_min = apply(data, 2, min, na.rm = TRUE),
    uniform_max = apply(data, 2, max, na.rm = TRUE)
  )

  if (has_uniform) {
    params$eta <- 0
  }

  return(list(
    params = params,
    responsibilities = matrix(1, n, 1),
    map_assignments = rep(1, n),
    outliers = rep(FALSE, n),
    log_likelihood_trace = c(-1000),  # Valeur arbitraire
    n_iter = 1,
    has_uniform = has_uniform
  ))
}
#' Model selection using BIC
select_best_model <- function(data, K, n_init = 20) {
  best_models <- list()
  BIC_values <- c(NA, NA)

  # Model without uniform component
  cat("Fitting model without uniform component...\n")
  best_no_uniform <- run_multiple_initializations(data, K, n_init, has_uniform = FALSE)

  if (!is.null(best_no_uniform)) {
    BIC_values[1] <- calculate_BIC(
      max(best_no_uniform$log_likelihood_trace),
      nrow(data),
      best_no_uniform$params,
      has_uniform = FALSE
    )
    best_models[[1]] <- best_no_uniform
  }

  # Model with uniform component
  cat("Fitting model with uniform component...\n")
  best_with_uniform <- run_multiple_initializations(data, K, n_init, has_uniform = TRUE)

  if (!is.null(best_with_uniform)) {
    BIC_values[2] <- calculate_BIC(
      max(best_with_uniform$log_likelihood_trace),
      nrow(data),
      best_with_uniform$params,
      has_uniform = TRUE
    )
    best_models[[2]] <- best_with_uniform
  }

  # Remove NA values
  valid_BIC <- !is.na(BIC_values)

  if (sum(valid_BIC) == 0) {
    stop("Both models failed to fit. Check your data and parameters.")
  }

  # Select best model based on BIC
  best_idx <- which.min(BIC_values[valid_BIC])
  # Map back to original indices (1 or 2)
  original_indices <- which(valid_BIC)
  best_original_idx <- original_indices[best_idx]

  best_model <- best_models[[best_original_idx]]
  best_model$BIC <- BIC_values[best_original_idx]
  best_model$has_uniform_selected <- (best_original_idx == 2)

  return(best_model)
}
