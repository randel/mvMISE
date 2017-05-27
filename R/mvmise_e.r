#' Multivariate mixed-effects selection model with correlated outcome-specific error terms
#' 
#' This function fits a multivariate mixed-effects selection model with potential missing values in the outcome and correlated outcome-specific error terms.
#' It can shrink the error precision matrix with a graphical lasso penalty 
#' for high-dimensional outcomes.
#' 
#' The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
#' is a multivariate mixed-effects model, with correlations among multivariate outcomes modelled via outcome-specific error terms. For the ith cluster,
#' the outcome \eqn{\mathbf{Y}_{i}} is a matrix of \eqn{n_i} samples (rows) and \eqn{K} outcomes (columns). 
#' Let \eqn{\mathbf{y}_{i} = \mathrm{vec}\left( \mathbf{Y}_{i} \right)}. 
#' The outcome vector \eqn{\mathbf{y}_{i}} can be modelled as
#' \deqn{\mathbf{y}_{i}  = \mathbf{X}_{i}\boldsymbol{\alpha}+\mathbf{Z}_{i}\mathbf{b}_{i}+\mathbf{e}_{i},}
#' where the random effects (\eqn{\mathbf{b}_{i}}) follow a normal distribution \eqn{\mathbf{b}_{i}\sim N(\mathbf{0},\mathbf{D})}; 
#' and the error term \eqn{\mathbf{e}_{i}=\mathrm{vec}\left(\mathbf{E}_{i}\right) \sim N(\mathbf{0},\boldsymbol{\Sigma}\otimes\mathbf{S}_{i})}.  
#' The matrix \eqn{\mathbf{S}_{i}} is an \eqn{n_i\times n_i} diagonal matrix with diagonal elements corresponding to the error variances of the \eqn{n_i} samples 
#' within the ith cluster. 
#' The variances for the first and other samples can be different if sigma_diff = TRUE. 
#' The matrix \eqn{\boldsymbol{\Sigma}} captures the error (or unexplained) covariances among \eqn{K} outcomes. 
#' To facilitate the computation for high-dimensional outcomes, the off-diagonal elements of the inverse of \eqn{\boldsymbol{\Sigma}} can be shrinked
#' by a graphical lasso penalty. If admm = TRUE (the default), the alternating direction method of multipliers (ADMM) is used to estimate \eqn{\boldsymbol{\Sigma}}.
#' The last fixed effect in \eqn{\boldsymbol{\alpha}} can be outcome-specific, if specific_eff
#' is specified as TRUE.
#' 
#' The missing-data model can be written as
#' \deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}\mathbf{1}_{n_{i}}^{'}\mathbf{y}_{ik}  + 
#' \phi_{2}\mathbf{1}_{n_{i}}^{'}\mathbf{x}_{i} \right),}
#' where \eqn{r_{ik}} is the missing indicator for the kth outcome in the ith cluster. If missing, the kth outcome in the ith cluster \eqn{\mathbf{y}_{ik}} 
#' is missing altogether.
#' The estimation is implemented via an EM algorithm. Parameters in the missing-data models can be specified via the argument miss_mechanism. If miss_mechanism 
#' = "y" or "yx", i.e., the missingness depends on the outcome, the missing-data mechanism is missing not at random (MNAR), 
#' otherwise it is missing at random (MAR).
#' 
#' It works for fully observed data if miss_mechanism = "none". It also works for univariate outcome with potential missing values, if the outcome Y is a matrix
#' with one column.
#' 
#' @param Y an outcome matrix, each row is an observation, each column is an outcome variable, with potential missing values (NAs).
#' @param X a covariates matrix, each row is an observation, each column is a covariate. Now covariates are assumed to be common for outcomes.
#' @param Zidx column indexes of matrix X used as the design matrix of random effects. The default is 1, i.e., a random intercept is included 
#' if the first column of X is a vector of 1s.
#' @param id a vector for cluster/grouping index, matching with the rows of Y and X.
#' @param maxIter maximum number of iterations for the EM algorithm.
#' @param admm logical. If TRUE (the default), the alternating direction method of multipliers (ADMM) is 
#' used to estimate the error precision matrix with a graphical lasso penalty. This works for multivariate outcomes. 
#' For an univariate outcome, it should be set as FALSE.
#' @param lambda tuning parameter for the graphical lasso penalty of the error precision matrix. It can be selected by AIC (an output).
#' @param tol tolerance level for the relative change in the observed-data log-likelihood function.
#' @param verbose logical. If TRUE, the iteration history of each step of the EM algorithm will be printed. The default is FALSE.
#' @param specific_eff logical. If TRUE, outcome-specific fixed-effects are estimated for the last covariate in X. The default is FALSE.
#' @param miss_mechanism one of "y" (the default), "x", "yx", and "none", indicating the missingness of outcome k in cluster i 
#' depends on the mean of the outcome, the mean of the covariate of interest, both, or none. The missing probability is modelled as exp(phi0 +
#' phi1*mean(y) + phi2*mean(x)). If there is no missing values in Y, it should be set as "none".
#' @param sigma_diff logical. If TRUE, the sample error variance of the first sample is different from that for the rest of samples within each cluster.
#' This is the case for the reference sample in the iTRAQ proteomics data. The default is FALSE.

#' @return A list containing
#' \item{beta}{the estimated fixed effects.}
#' \item{se}{the standard errors for the estimated fixed effects.}
#' \item{Sigma}{the estimated error covariance matrix for the outcomes.}
#' \item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
#'  the variances for the first sample and the rest of samples within each cluster.}
#' \item{D}{the estimated covariance matrix for the random effects.}
#' \item{phi}{the estimated parameters for the missing-data mechanism. The missing probability is modelled as exp(phi0 +
#' phi1*mean(y) + phi2*mean(x)). A zero value implies that parameter is ignored via the specification of miss_mechanism.}
#' \item{loglikelihood}{the observed-data log-likelihood values.}
#' \item{iter}{the number of iterations for the EM algorithm.}
#' \item{AIC}{The Akaike information criterion (AIC) calculated for selecting the tuning parameter lambda.}
#' 
#' @references Jiebiao Wang, Pei Wang, Donald Hedeker, and Lin S. Chen. A multivariate mixed-effects selection model framework for 
#' labelling-based proteomics data with non-ignorable missingness. (In preparation).
#' 
#' @export
#' @import lme4
#' @importFrom stats coef glm poisson
#' @examples
#' data(sim_dat)
#'
#' fit0 = mvmise_e(Y=sim_dat$Y, X=sim_dat$X, id=sim_dat$id)


mvmise_e = function(Y, X, Zidx = 1, id, maxIter = 100, tol = 0.001, lambda = 0.05, 
                    admm = TRUE, verbose = FALSE, specific_eff = FALSE, 
                    miss_mechanism = "y", sigma_diff = FALSE) {
  
  N = length(unique(id))  # no. clusters
  nY = ncol(Y) # no. of outcomes
  id = as.character(id)
  
  miss_cluster = NULL
  # number of observations/samples within each cluster
  named_vec = rep(NA, N)
  names(named_vec) = unique(id)
  ni = named_vec
  for (i in unique(id)) {
    ni[i] = sum(id == i)
    # within each cluster, ordered by clusters
    miss_cluster = c(miss_cluster, apply(Y[id == i, , drop = FALSE], 2, function(x) mean(is.na(x)) == 1))
  }
  
  ## covariate mean for missing data mechanism
  cov_mean = tapply(X[, ncol(X)], id, mean)
  cov_mean_rep = rep(cov_mean, rep(nY, length(cov_mean)))
  
  # stack y, X, and id by outcome variables
  y = as.vector(Y)
  X_mat = X[rep(1:nrow(X), nY), ]
  # assume the last covariate is of interest to allow outcome-specific effects
  if(specific_eff) {
    # outcome indicator
    Y_ind = kronecker(diag(nY), rep(1, nrow(Y)))
    X_mat = cbind(X_mat[, - ncol(X_mat)], X_mat[, ncol(X_mat)] * Y_ind)
  }
  id = rep(id, nY)
  
  # no. of covariates
  nX = ncol(X_mat) 
  
  # missing indicator: long vector
  miss = is.na(y)
  
  E_b_yobs = var_b_yobs = var_b_y = named_vec
  named_ls = vector("list", N)
  names(named_ls) = unique(id)
  Ei = var_y_yobs_lis = named_ls
  
  
  ## used to calculate starting values
  
  lmefit = lmer(y ~ -1 + X_mat + (1 | id))
  vc = as.data.frame(VarCorr(lmefit))
  
  Sigma = diag(nY)
  sigma2_0 = sigma2 = vc[vc$grp == "Residual", "vcov"]
  D = diag(vc[vc$grp == "id", "vcov"], length(Zidx))
  if (D == 0) D = 1
  beta = fixef(lmefit)
  
  # in case fixed-effect model matrix is rank deficient, dropping columns
  if (length(beta) != nX) {
    x_id = which(!paste0("X_mat", 1:nX) %in% names(fixef(lmefit)))
    X_mat = X_mat[, - x_id]
    nX = nX - length(x_id)
  }
  
  phi = rep(0, 3)
  
  iter = 0
  likelihood = NULL
  cond = TRUE
  
  while (cond) {
    iter = iter + 1
    
    ## E step
    
    obs_likelihood = 0
    XRX = SE = matrix(0, nX, nX)
    XRE = rep(0, nX)
    Y_ik_mean = matrix(NA, N, nY)
    rownames(Y_ik_mean) = unique(id)
    
    SS_sigma = matrix(0, nY * ni[i], nY * ni[i])
    
    for (i in unique(id)) {
      Xi = X_mat[id == i, , drop = F]
      Xio = Xi[!miss[id == i], , drop = F]
      Xim = Xi[miss[id == i], , drop = F]
      Zi = Xi[, Zidx, drop = F]
      
      Ri = kronecker(Sigma, diag(c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i]))
      # efficient solve(Ri)
      Ri_inv = kronecker(solve(Sigma), diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i]))
      
      Vi = Zi %*% D %*% t(Zi) + Ri
      if (sum(!miss[id == i]) > 0) 
        Vi_obs_inv = solve(Vi[!miss[id == i], !miss[id == i], drop = F]) else Vi_obs_inv = matrix(0, 0, 0)
      
      var_b_y[i] = solve(t(Zi) %*% Ri_inv %*% Zi + solve(D))
      
	  # calculating observed log-likelihood
      obs_phi = 0
      var_ymis_yobs = Vi[miss[id == i], miss[id == i], drop = F] - 
        Vi[miss[id == i], !miss[id == i], drop = F] %*% Vi_obs_inv %*% 
        Vi[!miss[id == i], miss[id == i], drop = F]
      E_ymis_yobs = Xim %*% beta + Vi[miss[id == i], !miss[id == i], drop = F] %*% Vi_obs_inv %*% 
        (y[!miss & id == i] - Xio %*% beta) + phi[2] * var_ymis_yobs %*% rep(1, sum(miss[id == i]))/ni[i]
      
      E_y_yobs = y[id == i]
      E_y_yobs[miss[id == i]] = E_ymis_yobs
      # for missing-data mechanism
      Y_ik_mean[i, ] = colMeans(matrix(E_y_yobs, nrow = ni[i]))
      
      var_y_yobs = matrix(0, ni[i] * nY, ni[i] * nY)  # n_i[i],n_i[i]
      var_y_yobs[miss[id == i], miss[id == i]] = var_ymis_yobs
      var_y_yobs_lis[[i]] = var_y_yobs
      
      E_b_yobs[i] = var_b_y[i] %*% t(Zi) %*% Ri_inv %*% (E_y_yobs - Xi %*% beta)
      var_b_yobs[i] = var_b_y[i] %*% t(Zi) %*% Ri_inv %*% var_y_yobs %*% 
        Ri_inv %*% Zi %*% var_b_y[i]
      
      XRX = XRX + t(Xi) %*% Ri_inv %*% Xi
      XRE = XRE + t(Xi) %*% Ri_inv %*% (E_y_yobs - Zi %*% E_b_yobs[i])
      
      Ei[[i]] = matrix(E_y_yobs - Xi %*% beta - Zi %*% E_b_yobs[i], ncol = nY)
      
      ## observed-data likelihood for missing data
      
	  # observed Yi as a matrix
      Y_k = matrix(y[id == i & !miss], nrow = ni[i])
      
      obs_phi = ifelse(sum(miss[id == i]) == ni[i], 0, sum(apply(Y_k, 2, function(x) 
        log(1 - min(exp(phi[1] + phi[2] * mean(x) + phi[3] * cov_mean[id]), 1 - .Machine$double.eps)))) + 
          ifelse(sum(miss[id == i]) == 0, 0, phi[1] * sum(miss[id == i])/ni[i] + 
                   phi[3] * cov_mean[id] * sum(miss[id == i])/ni[i] + 
                   phi[2] * (sum(E_ymis_yobs) + phi[2] * sum(var_ymis_yobs)/ni[i]/2)/ni[i]))
      
      ## observed-data likelihood (observed y; r)
      ## overflow if log(det(Vi[!miss[id==i],!miss[id==i],drop=F]))
      obs_likelihood = obs_likelihood - 
        0.5 * (determinant(x = Vi[!miss[id == i], !miss[id == i], drop = F], logarithm = TRUE)$modulus + 
                                                 t(y[!miss & id == i] - Xio %*% beta) %*% Vi_obs_inv %*% 
                                                 (y[!miss & id == i] - Xio %*% beta)) + obs_phi
      SE = SE + t(Xio) %*% Vi_obs_inv %*% Xio
    }
    
    Sigma_inv = solve(Sigma)
    # truncated at 1e-5
    AIC = -2 * obs_likelihood + 2 * sum(abs(Sigma_inv[lower.tri(Sigma_inv)]) > 1e-05)
	# penalized log-likelihood
    obs_likelihood = obs_likelihood - lambda * N * (sum(abs(Sigma_inv)) - 
                                                      sum(abs(diag(Sigma_inv))))/2
    likelihood = c(likelihood, obs_likelihood)
    
    # relative change of log-likelihood
    likelihood_change = ifelse(iter >= 2, (likelihood[iter] - likelihood[iter - 1])/abs(likelihood[iter - 1]), NA)
    cond = (iter < maxIter & ifelse(iter >= 2, abs(likelihood_change) > tol, TRUE))  
    
    ## M step
    
    # phi's
    if (miss_mechanism == 'y') {
      mylog <- glm(miss_cluster ~ as.vector(t(Y_ik_mean)), family = poisson(link = "log"))
      # phi's index to be updated
      phi_idx = -3
    }
    if (miss_mechanism == 'x') {
      mylog <- glm(miss_cluster ~ cov_mean_rep, family = poisson(link = "log"))
      phi_idx = -2
    }
    if (miss_mechanism == 'yx') {
      mylog <- glm(miss_cluster ~ as.vector(t(Y_ik_mean)) + cov_mean_rep, family = poisson(link = "log"))
      phi_idx = 1:3
    }
	# make the algorithm stable
    if (miss_mechanism != 'none') if (mylog$converged) phi[phi_idx] = as.numeric(coef(mylog))
    
    D_old = D
    D = (t(E_b_yobs) %*% E_b_yobs + sum(var_b_yobs + var_b_y))/N
	# make the algorithm stable
    if (D[1] > 1e+05) 
      D = D_old  
    beta = solve(XRX) %*% XRE
    
    # random intercept model
    
    SS_sigma = matrix(0, nY, nY)
    SS_s = SS_s0 = 0
    Sigma_inv = solve(Sigma)
    
    for (i in unique(id)) {
      Si_inv = diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i])
      Xi = X_mat[id == i, , drop = F]
      Zi = matrix(Xi[, Zidx, drop = F], ncol = nY)
      # efficient solve(Ri)
      Ri_inv = kronecker(solve(Sigma), diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i]))
      
      B = diag(ni[i] * nY) - Xi[, Zidx, drop = F] %*% var_b_y[i] %*% 
        t(Xi[, Zidx, drop = F]) %*% Ri_inv
      tmp = kpsvd(B %*% var_y_yobs_lis[[i]] %*% t(B), nY, ni[i])
      
      # for Sigma
      A0 = matrix(0, nY, nY)
      for (k in 1:length(tmp$d)) A0 = A0 + tmp$d[k] * t(matrix(tmp$u[, k], nY, nY)) * 
        sum(diag(Si_inv %*% matrix(tmp$v[, k], ni[i], ni[i])))
      
      SS_sigma = SS_sigma + (var_b_y[i]) * t(Zi) %*% Si_inv %*% Zi + 
        t(Ei[[i]]) %*% Si_inv %*% Ei[[i]] + A0
      
      # for Si
      A0 = matrix(0, ni[i], ni[i])
      for (k in 1:length(tmp$d)) A0 = A0 + tmp$d[k] * t(matrix(tmp$v[, k], ni[i], ni[i])) * 
        sum(diag(Sigma_inv %*% matrix(tmp$u[, k], nY, nY)))
      
      Si = (var_b_y[i]) * Zi %*% Sigma_inv %*% t(Zi) + Ei[[i]] %*% 
        Sigma_inv %*% t(Ei[[i]]) + A0
      SS_s0 = SS_s0 + Si[1, 1]
      SS_s = SS_s + sum(diag(Si)[-1])
    }
    
    ### Sigma
    
    if (admm) 
      Sigma = admm(N, ni, SS_sigma, lambda, rho = 1, maxIter = 500) else {
        # unpenalized covariance matrix for outcomes
        Sigma = SS_sigma/sum(ni)
      }
    
    ### sigma
    
    if (sigma_diff) {
      
      ### sigma2_0, sigma2
      
      sigma2_0 = SS_s0/(N * nY)
      sigma2 = SS_s/((sum(ni) - N) * nY)
      
      # identification issue for kronecker product: set sigma2_0 as 1
      sigma2 = sigma2/sigma2_0
      sigma2_0 = 1
    } else {
      sigma2_0 = sigma2 = (SS_s0 + SS_s) / (sum(ni) * nY)
    }
    
    if (verbose) {
      print(round(c(iter = iter, logLike_change = likelihood_change, 
                    beta = beta, sigma2 = sigma2, 
                    Sigma = mean(diag(Sigma)), D = D, phi = phi), 3))
    }
  }
  
  if (sigma_diff) sigma2 = c(sigma2_0, sigma2)
  
  ## standard errors for fixed-effects
  se = sqrt(diag(solve(SE)))
  
  return(list(iter = iter, beta = beta, se = se, sigma2 = sigma2, 
              Sigma = Sigma, D = D, phi = phi, AIC = AIC, loglikelihood = obs_likelihood))
}


## Kronecker product SVD: write matrix A in Kronecker product thus we can
## have closed-form derivations w.r.t variance matrix (Van Loan, 2000)

kpsvd = function(A, nY, ni_i) {
  # permutation matrix
  RA = matrix(NA, nY^2, ni_i^2)  
  k = 1
  for (i in 1:nY) {
    for (j in 1:nY) {
      RA[k, ] = as.vector(A[1:ni_i + (j - 1) * ni_i, 1:ni_i + (i - 1) * ni_i])
      k = k + 1
    }
  }
  
  return(svd(RA))
}


# ADMM algorithm for estimating error covariance, Sigma: var(errors)

admm = function(N, ni, SS_sigma, lambda, rho = 1, maxIter) {
  nY = nrow(SS_sigma)
  Theta = diag(nY)
  U = Z = matrix(0, nY, nY)
  
  Theta_old = matrix(100, nY, nY)
  iter = 1
  while (iter < maxIter & sum(abs(Theta_old - Theta))/sum(abs(Theta_old)) >= 1e-05) {
    eigen0 = eigen((N * rho * (U - Z) + SS_sigma)/sum(ni))
    D = sum(ni) * (-eigen0$values + sqrt(eigen0$values^2 + 4 * N * 
                                           rho/sum(ni)))/(2 * N * rho)
    Theta = eigen0$vectors %*% diag(D) %*% t(eigen0$vectors)
    
    A = Theta + U
    Z = sign(A) * pmax(abs(A) - lambda/rho, 0)
    diag(Z) = diag(A)
    U = U + Theta - Z
    iter = iter + 1
  }
  
  return(solve(Theta))
}

