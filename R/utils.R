#' Initialize Gaussian Mixture Model parameters
#'
#' @param data Matrix of data (n x d)
#' @param K Number of Gaussian components
#' @return List with initialized parameters
initialize_parameters <- function(data, K) {
  n <- nrow(data)
  d <- ncol(data)

  # Calculer la covariance globale
  global_cov <- cov(data)

  # Régulariser la covariance globale
  global_cov <- global_cov + diag(d) * 0.1

  if (K == 1) {
    # Pour K=1, utiliser la covariance globale
    params <- list(
      pi = 1,
      mu = matrix(colMeans(data), nrow = 1),
      sigma = array(global_cov, dim = c(d, d, 1))
    )
  } else {
    # Initialisation k-means
    km <- kmeans(data, centers = K, nstart = 10)

    params <- list(
      pi = rep(1/K, K),
      mu = km$centers,
      sigma = array(0, dim = c(d, d, K))
    )

    for (k in 1:K) {
      cluster_data <- data[km$cluster == k, , drop = FALSE]

      if (nrow(cluster_data) > d + 1) {
        # Assez de points pour estimer la covariance
        cluster_cov <- cov(cluster_data)
      } else {
        # Pas assez de points, utiliser la covariance globale
        cluster_cov <- global_cov
      }

      # Régulariser
      cluster_cov <- cluster_cov + diag(d) * 0.1

      # Vérifier que la matrice est définie positive
      eigvals <- eigen(cluster_cov, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigvals) < 1e-6) {
        cluster_cov <- cluster_cov + diag(d) * (1e-6 - min(eigvals))
      }

      params$sigma[, , k] <- cluster_cov
    }
  }

  # Add uniform component parameters
  # Ajouter une marge pour éviter les densités infinies
  params$uniform_min <- apply(data, 2, min) - 1
  params$uniform_max <- apply(data, 2, max) + 1

  return(params)
}

#' Calculate log-likelihood
#'
#' @param data Data matrix
#' @param params Model parameters
#' @param has_uniform Boolean for uniform component
#' @return Log-likelihood value
calculate_log_likelihood <- function(data, params, has_uniform = FALSE) {
  n <- nrow(data)
  K <- length(params$pi)
  d <- ncol(data)

  # 1. Calculer les densités gaussiennes (non log)
  densities <- matrix(0, n, K)

  for (k in 1:K) {
    sigma_k <- params$sigma[, , k]

    # Assurer que la matrice est valide
    sigma_k <- (sigma_k + t(sigma_k)) / 2
    sigma_k <- sigma_k + diag(d) * 1e-6

    densities[, k] <- mvtnorm::dmvnorm(
      data,
      mean = params$mu[k, ],
      sigma = sigma_k,
      log = FALSE  # FALSE pour obtenir les densités, pas les log-densités
    )
  }

  # 2. Ajouter la composante uniforme si nécessaire
  if (has_uniform) {
    uniform_dens <- 1 / prod(params$uniform_max - params$uniform_min)
    if (!is.finite(uniform_dens) || uniform_dens <= 0) {
      uniform_dens <- 1e-10
    }

    densities <- cbind(densities, rep(uniform_dens, n))
    weights <- c(params$pi, params$eta)
  } else {
    weights <- params$pi
  }

  # Calculer les densités du mélange pondérées
  weighted_densities <- sweep(densities, 2, weights, "*")

  # Somme sur les composantes pour chaque point
  mixture_densities <- rowSums(weighted_densities)

  # Éviter les zéros numériques
  mixture_densities <- pmax(mixture_densities, 1e-100)

  # Calculer la log-vraisemblance
  log_likelihood <- sum(log(mixture_densities))

  return(log_likelihood)
}
