simulate_covariance <- function(n, rate = 1,
	weak_range = c(-0.1, 0.1),
	strong_range = c(0.4, 0.8),
	split_quantile = 0.5) {
	
	# 1. Draw standard deviations (σ_i ~ Exp(rate))
	sigma <- rexp(n, rate = rate)
	#sigma <- rep(1/rate, n)
    sigma[length(sigma)] <- sigma[length(sigma)]
	
	# 2. Initialize correlation matrix
	#Rho <- matrix(runif(n * n, -1, 1), n, n)
	Rho <- matrix(runif(n * n, 0.5, 1), n, n)
	Rho <- (Rho + t(Rho)) / 2
	diag(Rho) <- 1
	
	# 3. Quantile-based split index
	split_index <- ceiling(n * split_quantile)
	
	weak_idx <- 1:split_index
	strong_idx <- min((split_index+1),n):n
	
	# 4. Adjust correlations for last column
	Rho[weak_idx, n] <- runif(length(weak_idx), weak_range[1], weak_range[2])
	Rho[strong_idx, n] <- runif(length(strong_idx), strong_range[1], strong_range[2])
	Rho[n, ] <- Rho[, n]
    write.table(round(Rho,3), file.path(outdir, "Rho_old.tbl"))
	
	# 5. Build covariance
	D <- diag(sigma)
	Sigma <- D %*% Rho %*% D
	Sigma <- as.matrix(Matrix::nearPD(Sigma)$mat)
    Rho <- cov2cor(Sigma)
	
	return(list(Sigma=Sigma, Rho=Rho))
}


simulate_covariance2 <- function(n, rate = 1,
                                weak_range = c(-0.1, 0.1),
                                strong_range = c(0.4, 0.8),
                                split_quantile = 0.5) {
  # 1. Split indices
  split_idx <- ceiling(n * split_quantile)
  weak_idx   <- 1:split_idx
  strong_idx <- (split_idx + 1):n
  
  # 2. Build factor loadings that produce desired correlations
  loadings <- rep(0, n)
  loadings[weak_idx]   <- runif(length(weak_idx), weak_range[1], weak_range[2])
  loadings[strong_idx] <- runif(length(strong_idx), strong_range[1], strong_range[2])
    if(length(loadings)>n){
        loadings <- loadings[-length(loadings)]
    }
  
  # 3. Correlation from single common factor + unique variance
  Rho <- outer(loadings, loadings)
  diag(Rho) <- 1
  # unique variance so total correlation structure realistic
  Rho <- Rho + diag(1 - diag(Rho))
  
  # slight random jitter to avoid perfect structure
  Rho <- Rho + matrix(rnorm(n^2, 0, 0.02), n, n)
  Rho <- (Rho + t(Rho))/2
  diag(Rho) <- 1
  write.table(round(Rho,3), file.path(outdir, "Rho_old.tbl"))
  
  # 4. Force PD (small shrink only)
  Rho <- as.matrix(Matrix::nearPD(Rho, corr = TRUE)$mat)
  
  # 5. Create covariance with arbitrary SDs
  sigma <- rexp(n, rate)
    print("Sigma:")
    print(sigma)
  Sigma <- diag(sigma) %*% Rho %*% diag(sigma)
  
  return(list(Rho = Rho, Sigma = Sigma))
}
