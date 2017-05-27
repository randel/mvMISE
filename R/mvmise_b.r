#' Multivariate mixed-effects selection model with correlated outcome-specific random intercepts
#' 
#' This function fits a multivariate mixed-effects selection model with potential missing values in the outcome
#'  and correlated outcome-specific random intercepts.
#' 
#' The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
#' is a multivariate mixed-effects model, with correlations among multivariate outcomes modelled via outcome-specific random intercepts with 
#' a factor-analytic structure
#' \deqn{\mathbf{y}_{i} =  \mathbf{X}_{i}\boldsymbol{\alpha}+\mathbf{Z}_{i}\boldsymbol{\tau}b_{i}+\mathbf{e}_{i},}
#' where \eqn{i} denotes a cluster, \eqn{\boldsymbol{\tau}} is a \eqn{K\times 1} vector for the outcome-specific variance components corresponding to 
#' the random effect \eqn{b_i} (a standard normal random variable), and \eqn{K} is the number of outcomes. 
#' The factor-analytic structure is used to facilitate the computation.
#' It assumes that the random effects 
#' are derived from a latent variable \eqn{b_i} with a loading vector \eqn{\boldsymbol{\tau}}.
#' In this way, only \eqn{K} rather parameters are needed 
#' in the estimation for the covariance matrix of random effects. The last fixed effect in \eqn{\boldsymbol{\alpha}} can be outcome-specific, if specific_eff
#' is specified as TRUE.
#' 
#' The missing-data model can be written as
#' \deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}\mathbf{1}_{n_{i}}^{'}\mathbf{y}_{ik}  + 
#' \phi_{2}\mathbf{1}_{n_{i}}^{'}\mathbf{x}_{i} \right),}
#' where \eqn{r_{ik}} is the missing indicator for the kth outcome in the ith cluster. If missing, the kth outcome in the ith cluster  \eqn{\mathbf{y}_{ik}} 
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
#' @param Z a design matrix for random effects, each row is an observation, each column is a random effect.
#'    If it is NULL (the default), a matrix with each column as an indicator for each outcome is generated.
#' @param id a vector for cluster/grouping index, matching with the rows of Y, X, Z (if specified).
#' @param maxIter maximum number of iterations for the EM algorithm.
#' @param tol tolerance level for the relative change in the observed-data log-likelihood function.
#' @param verbose logical. If TRUE, the iteration history of each step of the EM algorithm will be printed. The default is FALSE.
#' @param specific_eff logical. If TRUE, outcome-specific fixed-effects are estimated for the last covariate in X. The default is FALSE.
#' @param miss_mechanism one of "y" (the default), "x", "yx", and "none", indicating the missingness of outcome k in cluster i 
#' depends on the mean of the outcome, the mean of the covariate of interest, both, or none. The missing probability is modelled as exp(phi0 +
#' phi1*mean(y) + phi2*mean(x)). If there is no missing values in Y, it should be set as "none".
#' @param sigma_diff logical. If TRUE, the sample error variance of the first sample is different from that for the rest of samples within each cluster.
#' This is the case for the reference sample in the iTRAQ proteomics data. The default is FALSE.
#' 
#' @return A list containing
#' \item{beta}{the estimated fixed effects.}
#' \item{se}{the standard errors for the estimated fixed effects.}
#' \item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
#'  the variances for the first sample and the rest of samples within each cluster.}
#' \item{tau}{the estimated variance components for the outcome-specific factor-analytic random effects.}
#' \item{phi}{the estimated parameters for the missing-data mechanism. The missing probability is modelled as exp(phi0 +
#' phi1*mean(y) + phi2*mean(x)). A zero value implies that parameter is ignored via the specification of miss_mechanism.}
#' \item{loglikelihood}{the observed-data log-likelihood values.}
#' \item{iter}{the number of iterations for the EM algorithm.}
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
#' fit0 = mvmise_b(Y=sim_dat$Y, X=sim_dat$X, id=sim_dat$id)


mvmise_b = function(Y, X, Z = NULL, id, maxIter = 100, tol = 0.001, verbose = FALSE, specific_eff = FALSE,
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
  
  # starting values of missing data parameters, default is all 0
  phi = rep(0, 3)
  phi_old = phi
  
  named_ls = vector("list", N)
  names(named_ls) = unique(id)
  E_e_yobs = var_e_yobs = named_ls
  
  ## used to calculate starting values
  
  lmefit = lmer(y ~ -1 + X_mat + (1 | id), REML = FALSE)
  tau = rep(1, nY)
  
  sigma2_0 = sigma2 = 1
  beta = fixef(lmefit)
  
  # in case fixed-effect model matrix is rank deficient, dropping columns
  if (length(beta) != nX) {
    x_id = which(!paste0("X_mat", 1:nX) %in% names(fixef(lmefit)))
    X_mat = X_mat[, - x_id]
    nX = nX - length(x_id)
  }
  
  iter = 0
  likelihood = NULL
  cond = TRUE
  
  while (cond) {
    iter = iter + 1
    
    ## E step
    
    obs_likelihood = 0
    XRX = SE = matrix(0, nX, nX)
    
    XRE = rep(0, nX)
    ZRZ = matrix(0, nY, nY)
    ZRE = rep(0, nY)
    Y_ik_mean = matrix(NA, N, nY)
    rownames(Y_ik_mean) = unique(id)
    
    for (i in unique(id)) {
      
      Xi = X_mat[id == i, , drop = F]
      Xio = Xi[!miss[id == i], , drop = F]
      Xim = Xi[miss[id == i], , drop = F]
      if (is.null(Z)) 
        Zi = kronecker(diag(nY), rep(1, ni[i])) else Zi = Z[id==i, ]
      
      sigma = c(sigma2_0, rep(sigma2, ni[i] - 1))
      Ri = kronecker(diag(nY), diag(sigma, ni[i]))
      Ri_inv = kronecker(diag(nY), diag(1/sigma, ni[i]))
      
      Ti = Zi %*% tau
      TT = Ti %*% t(Ti)
      Vi = TT + Ri
      Ri_inv_Ti = Ri_inv %*% Ti
      var_b_y = as.numeric(1/(1 + t(Ti) %*% Ri_inv_Ti))
      Vi_inv = Ri_inv - var_b_y * Ri_inv_Ti %*% t(Ri_inv_Ti)
      
      Vi_mis_obs = Vi[miss[id == i], !miss[id == i], drop = F]
      Vi_obs = Vi[!miss[id == i], !miss[id == i], drop = F]
      Tio = Ti[!miss[id == i]]
      Ri_obs_inv = Ri_inv[!miss[id == i], !miss[id == i]]
      if (sum(!miss[id == i]) > 0) 
        Vi_obs_inv = Ri_obs_inv - 1/(1 + as.numeric(t(Tio) %*% Ri_obs_inv %*% Tio)) * 
        Ri_obs_inv %*% Tio %*% t(Tio) %*% 
        Ri_obs_inv else Vi_obs_inv = matrix(0, 0, 0)
      
	  # observed Yi as a matrix
      Y_k = matrix(y[id == i & !miss], nrow = ni[i])  
      var_ymis_yobs = Vi[miss[id == i], miss[id == i]] - Vi_mis_obs %*% 
        Vi_obs_inv %*% t(Vi_mis_obs)
      E_ymis_yobs = Xim %*% beta + Vi_mis_obs %*% Vi_obs_inv %*% 
        (y[!miss & id == i] - Xio %*% beta) + phi[2] * var_ymis_yobs %*% 
        rep(1, sum(miss[id == i]))/ni[i]
      
      E_y_yobs = y[id == i]
      E_y_yobs[miss[id == i]] = E_ymis_yobs
      Y_ik_mean[i, ] = colMeans(matrix(E_y_yobs, nrow = ni[i]))
      
      var_y_yobs = matrix(0, ni[i] * nY, ni[i] * nY)
      var_y_yobs[miss[id == i], miss[id == i]] = var_ymis_yobs
      
      Ri_Vi_inv = diag(nY * ni[i]) - var_b_y * TT %*% Ri_inv
      var_e_yobs[[i]] = var_b_y * TT + Ri_Vi_inv %*% var_y_yobs %*% 
        t(Ri_Vi_inv)
      
      Bi = var_b_y %*% t(Ri_inv_Ti)
      E_b_yobs = as.numeric(Bi %*% (E_y_yobs - Xi %*% beta))
      
      XRX = XRX + t(Xi) %*% Ri_inv %*% Xi
      XRE = XRE + t(Xi) %*% Ri_inv %*% (E_y_yobs - Ti %*% E_b_yobs)
      ZRZ = ZRZ + E_b_yobs^2 * t(Zi) %*% Ri_inv %*% Zi
      ZRE = ZRE + E_b_yobs * t(Zi) %*% Ri_inv %*% (E_y_yobs - Xi %*% beta)
      
      E_e_yobs[[i]] = E_y_yobs - Xi %*% beta - Ti %*% E_b_yobs
      
      ## observed-data likelihood related to missing data
      
      obs_phi = ifelse(sum(miss[id == i]) == ni[i], 0, sum(apply(Y_k, 2, function(x) 
        log(1 - min(exp(phi[1] + phi[2] * mean(x) + phi[3] * cov_mean[id]), 1 - .Machine$double.eps)))) + 
          ifelse(sum(miss[id == i]) == 0, 0, phi[1] * sum(miss[id == i])/ni[i] + 
                   phi[3] * cov_mean[id] * sum(miss[id == i])/ni[i] + 
                   phi[2] * (sum(E_ymis_yobs) + phi[2] * sum(var_ymis_yobs)/ni[i]/2)/ni[i]))
      
      ## observed-data likelihood (observed y; r)
      obs_likelihood = obs_likelihood - 0.5 * (ifelse(ncol(Vi_obs) == 0, 0, sum(log(eigen(Vi_obs)$values))) + 
                                                 t(y[!miss & id == i] - Xio %*% beta) %*% Vi_obs_inv %*% 
                                                 (y[!miss & id == i] - Xio %*% beta)) + obs_phi
      
      SE = SE + t(Xio) %*% Vi_obs_inv %*% Xio
    }
    
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
    
    beta = solve(XRX) %*% XRE
    tau = solve(ZRZ) %*% ZRE
    
    SS_s = SS_s0 = 0
    for (i in unique(id)) {
      Si = E_e_yobs[[i]] %*% t(E_e_yobs[[i]]) + var_e_yobs[[i]]
      SS_s0 = SS_s0 + sum(diag(Si)[1 + (1:nY - 1) * ni[i]])
      SS_s = SS_s + sum(diag(Si)[-(1 + (1:nY - 1) * ni[i])])
    }
   
	### sigma
    
    if (sigma_diff) {
      
      ### sigma2_0, sigma2
      
      sigma2_0 = SS_s0/(N * nY)
      sigma2 = SS_s/((sum(ni) - N) * nY)

    } else {
      sigma2_0 = sigma2 = (SS_s0 + SS_s) / (sum(ni) * nY)
    }
    
    if (verbose) 
      print(round(c(iter = iter, logLike_change = likelihood_change, 
                    beta = beta, sigma2_0 = sigma2_0, sigma2 = sigma2, tau = tau[1], 
                    phi = phi), 3))
  }
  
  if (sigma_diff) sigma2 = c(sigma2_0, sigma2)
  
  ## standard errors for fixed-effects
  se = sqrt(diag(solve(SE)))
  
  return(list(iter = iter, beta = beta, se = se, sigma2 = sigma2, 
              tau = tau, phi = phi, loglikelihood = likelihood[length(likelihood)]))
}

