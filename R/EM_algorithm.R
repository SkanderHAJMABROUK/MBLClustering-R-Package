#' EM Algorithm for Gaussian Mixture Model
#'
#' Performs Expectation-Maximization algorithm for Gaussian Mixture Models
#' with optional uniform component for outlier detection.
#'
#' @param data Data matrix (n x d)
#' @param K Number of Gaussian components
#' @param max_iter Maximum iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param has_uniform Include uniform component (default: FALSE)
#' @param init_params Optional initial parameters list. If NULL, parameters
#'   are initialized automatically (default: NULL)
#'
#' @return List with EM results including:
#' \item{params}{Estimated parameters}
#' \item{responsibilities}{Component responsibilities matrix}
#' \item{map_assignments}{Cluster assignments by MAP}
#' \item{outliers}{Boolean vector indicating outliers}
#' \item{log_likelihood_trace}{Log-likelihood values during iterations}
#' \item{n_iter}{Number of iterations performed}
#' \item{has_uniform}{Boolean indicating if uniform component was used}
#'
#' @keywords internal
#' @noRd
EM_algorithm <- function(data, K, max_iter = 100, tol = 1e-6, has_uniform = FALSE, init_params = NULL) {
  n <- nrow(data)
  d <- ncol(data)

  # Initialize parameters
  if (!is.null(init_params)) {
    params <- init_params
  } else {
    params <- initialize_parameters(data, K)
  }

  if (has_uniform) {
    params$eta <- 0.1  # Initial proportion for uniform component
    params$pi <- rep((1 - params$eta) / K, K)
  }

  log_likelihood_trace <- numeric(max_iter)

  for (iter in 1:max_iter) {
    # E-step: Calculate responsibilities
    log_resp <- E_step(data, params, has_uniform)

    # M-step: Update parameters
    params <- M_step(data, log_resp, params, has_uniform)

    # Calculate log-likelihood
    log_likelihood_trace[iter] <- calculate_log_likelihood(data, params, has_uniform)

    # Check convergence
    if (iter > 1 && abs(log_likelihood_trace[iter] - log_likelihood_trace[iter-1]) < tol) {
      log_likelihood_trace <- log_likelihood_trace[1:iter]
      break
    }
  }

  # Calculate final responsibilities and MAP assignments
  log_resp <- E_step(data, params, has_uniform)
  resp <- exp(log_resp)
  map_assignments <- apply(resp, 1, which.max)

  # Identify outliers if uniform component exists
  outliers <- NULL
  if (has_uniform) {
    outliers <- map_assignments == (K + 1)
  }

  return(list(
    params = params,
    responsibilities = resp,
    map_assignments = map_assignments,
    outliers = outliers,
    log_likelihood_trace = log_likelihood_trace,
    n_iter = iter,
    has_uniform = has_uniform
  ))
}

#' E-step: Calculate log-responsibilities
#'
#' @param data Data matrix
#' @param params Model parameters
#' @param has_uniform Boolean indicating if uniform component is included
#' @return Matrix of log-responsibilities
#' @keywords internal
E_step <- function(data, params, has_uniform) {
  n <- nrow(data)
  K <- length(params$pi)
  d <- ncol(data)

  # Gaussian components log-densities
  log_dens_gaussian <- matrix(0, n, K)
  for (k in 1:K) {
    # Vérifier que la matrice de covariance est valide
    sigma_k <- params$sigma[, , k]

    # S'assurer que la matrice est définie positive
    if (any(is.na(sigma_k)) || any(is.infinite(sigma_k))) {
      sigma_k <- diag(d) * 0.1
    }

    # Vérifier que la matrice est symétrique
    sigma_k <- (sigma_k + t(sigma_k)) / 2

    # Ajouter une régularisation pour éviter les singularités
    sigma_k <- sigma_k + diag(d) * 1e-6

    log_dens_gaussian[, k] <- mvtnorm::dmvnorm(
      data,
      mean = params$mu[k, ],
      sigma = sigma_k,
      log = TRUE
    ) + log(params$pi[k])
  }

  if (has_uniform) {
    # Uniform component - calculer la densité uniforme
    uniform_dens <- 1 / prod(params$uniform_max - params$uniform_min)

    # Vérifier que la densité uniforme est valide
    if (!is.finite(uniform_dens) || uniform_dens <= 0) {
      uniform_dens <- 1e-10
    }

    log_dens_uniform <- log(uniform_dens) + log(params$eta)
    log_dens <- cbind(log_dens_gaussian, log_dens_uniform)
  } else {
    log_dens <- log_dens_gaussian
  }

  # Normalize with log-sum-exp trick for numerical stability
  max_log_dens <- apply(log_dens, 1, max)
  log_resp <- sweep(log_dens, 1, max_log_dens, "-")
  resp <- exp(log_resp)

  # Normalize les responsabilités
  resp <- resp / rowSums(resp)

  # Retourner en log pour la stabilité
  log_resp <- log(resp)

  return(log_resp)
}

#' M-step: Update parameters
#'
#' @param data Data matrix
#' @param log_resp Log-responsibilities from E-step
#' @param params Model parameters
#' @param has_uniform Boolean indicating if uniform component is included
#' @return Updated parameters
#' @keywords internal
M_step <- function(data, log_resp, params, has_uniform) {
  n <- nrow(data)
  d <- ncol(data)
  resp <- exp(log_resp)
  K <- ncol(resp) - as.integer(has_uniform)

  # Update mixing proportions
  if (has_uniform) {
    resp_gaussian <- resp[, 1:K, drop = FALSE]
    resp_uniform <- resp[, K + 1, drop = FALSE]

    N_gaussian <- colSums(resp_gaussian)
    N_uniform <- sum(resp_uniform)

    params$pi <- N_gaussian / n
    params$eta <- N_uniform / n

    # Adjust Gaussian proportions
    params$pi <- params$pi * (1 - params$eta)
  } else {
    N <- colSums(resp)
    params$pi <- N / n
  }

  # Update Gaussian means
  for (k in 1:K) {
    params$mu[k, ] <- colSums(resp[, k] * data) / sum(resp[, k])
  }

  # Update Gaussian covariances
  for (k in 1:K) {
    centered <- sweep(data, 2, params$mu[k, ], "-")
    weighted_centered <- sqrt(resp[, k]) * centered
    params$sigma[, , k] <- (t(weighted_centered) %*% weighted_centered) / sum(resp[, k])

    # Add regularization for numerical stability
    params$sigma[, , k] <- params$sigma[, , k] + diag(d) * 1e-6
  }

  return(params)
}

