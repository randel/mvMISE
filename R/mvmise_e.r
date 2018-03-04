#' A multivariate mixed-effects selection model with correlated outcome-specific error terms
#' 
#' This function fits a multivariate mixed-effects selection model with correlated outcome-specific error terms and potential missing values in the outcome.
#' Here an outcome refers to a response variable, for example, a genomic feature. The proposed model and function jointly analyze multiple outcomes/features.
#' For high-dimensional outcomes, the model can regularize the estimation by shrinking the error precision matrix with a graphical lasso penalty.
#' Given the introduction of the penalty and the choice of tuning parameter often being data-dependant, we recommend using permutation to calculate
#' p-values for testing with the mvMISE_e model. Please see mvMISE_e_perm for calculating the permutation-based p-values.
#' 
#' The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
#' is a multivariate mixed-effects model. The correlations among multivariate outcomes are modeled via outcome-specific error terms with an unstructured covariance matrix. 
#' For the i-th cluster, the outcome matrix \eqn{\mathbf{Y}_{i}} is a matrix of \eqn{n_i} samples (rows) and \eqn{K} outcomes (columns). 
#' Let \eqn{\mathbf{y}_{i} = \mathrm{vec}\left( \mathbf{Y}_{i} \right)}. 
#' The outcome vector \eqn{\mathbf{y}_{i}} can be modelled as
#' \deqn{\mathbf{y}_{i}  = \mathbf{X}_{i}\boldsymbol{\beta}+\mathbf{Z}_{i}\mathbf{b}_{i}+\mathbf{e}_{i},}
#' where the random effects (\eqn{\mathbf{b}_{i}}) follow a normal distribution \eqn{\mathbf{b}_{i}\sim N(\mathbf{0},\mathbf{D})}; 
#' and the error term \eqn{\mathbf{e}_{i}=\mathrm{vec}\left(\mathbf{E}_{i}\right) \sim N(\mathbf{0},\boldsymbol{\Sigma}\otimes\mathbf{S}_{i})}.  
#' The matrix \eqn{\mathbf{S}_{i}} is an \eqn{n_i\times n_i} diagonal matrix with diagonal elements corresponding to the error variances of the \eqn{n_i} samples 
#' within the i-th cluster. 
#' The variances for the first and other samples can be different if sigma_diff = TRUE. 
#' The matrix \eqn{\boldsymbol{\Sigma}} captures the error (or unexplained) covariances among the \eqn{K} outcomes. 
#' For high-dimensional outcomes, if admm = TRUE (the default), the off-diagonal elements of the inverse of \eqn{\boldsymbol{\Sigma}} will be shrinked
#' by a graphical lasso penalty and the alternating direction method of multipliers (ADMM) is used to estimate \eqn{\boldsymbol{\Sigma}}. 
#' If admm = FALSE, no penalty is used to estimate the unstructured error covariance matrix, and that is 
#' only applicable to low-dimensional multivariate outcomes.
#' 
#' The missing-data model can be written as
#' \deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}/n_{i}\cdot \mathbf{1}^{'}\mathbf{y}_{ik}  + 
#' \phi_{2}/n_{i}\cdot \mathbf{1}^{'}\mathbf{c}_{i} \right),}
#' where \eqn{r_{ik}} is the missing indicator for the k-th outcome in the i-th cluster. If missing \eqn{r_{ik}=1}, the k-th outcome in the i-th cluster \eqn{\mathbf{y}_{ik}} 
#' is missing altogether.
#' The estimation is implemented within an EM algorithm framework. Parameters in the missing-data models can be specified via the arguments miss_y and cov_miss. If miss_y 
#' = TURE, the missingness depends on the outcome values. 
#' If cov_miss is specified, the missingness can (additionally) depend on the specified covariates (cov_miss).
#' 
#' The model also works for fully observed data if miss_y = FALSE and cov_miss = NULL. It would also work for an univariate outcome with potential missing values, if the outcome Y is a matrix
#' with one column.
#' 
#' @param Y an outcome matrix. Each row is a sample, and each column is an outcome variable, with potential missing values (NAs).
#' @param X a covariate matrix. Each row is a sample, and each column is a covariate. The covariates can be common among all of the outcomes (e.g., age, gender) or outcome-specific.
#'    If a covariate is specific for the k-th outcome, one may set all the values corresponding to the other outcomes to be zero. If X is common across outcomes, the row number of X equals 
#'    the row number of Y. Otherwise if X is outcome-specific, the row number of X equals the number of elements in Y, i.e., outcome-specific X is stacked across outcomes within
#'    each cluster. See the Examples for demonstration.
#' @param id a vector for cluster/batch index, matching with the rows of Y, and X if it is not outcome specific.
#' @param Zidx the column indices of matrix X used as the design matrix of random effects. The default is 1, i.e., a random intercept is included 
#' if the first column of X is a vector of 1s. If Zidx=c(1,2), then the model would estimate the random intercept and the random effects of the 2nd column in the covariate matrix X.
#' The random-effects in this model are assumed to be independent.
#' @param maxIter the maximum number of iterations for the EM algorithm.
#' @param admm logical. If TRUE (the default), we impose a L1 graphical lasso penalty on the error precision (inverse of covariance) matrix, and the alternating direction method of multipliers (ADMM) is 
#' used to estimate the error precision and the error covariance matrix. If FALSE, no penalty is used to estimate the unstructured error covariance matrix, and that is 
#' only applicable to low-dimensional multivariate outcomes.
#' For an univariate outcome, it should be set as FALSE.
#' @param lambda the tuning parameter for the graphical lasso penalty of the error precision matrix. It can be selected by AIC (an output).
#' @param tol the tolerance level for the relative change in the observed-data log-likelihood function.
#' @param verbose logical. If TRUE, the iteration history of each step of the EM algorithm will be printed. The default is FALSE.
#' @param miss_y logical. If TRUE, the missingness depends on the outcome Y (see the Details). The default is TRUE.
#'      This outcome-dependent missing data pattern was motivated by and was observed in the mass-spectrometry-based quantitative proteomics data.  
#' @param cov_miss the covariate that can be used in the missing-data model. If it is NULL, 
#'    the missingness is assumed to be independent of the covariates.
#'   Check the Details for the missing-data model.
#'   If it is specified and the covariate is not outcome specific, its length equals the length of id. If it is outcome specific, the outcome-specific covariate is stacked across outcomes within
#'   each cluster.
#' @param sigma_diff logical. If TRUE, the sample error variance of the first sample is different from that for the rest of samples within each cluster.
#' This option is designed and used when analyzing batch-processed proteomics data with the first sample in each cluster/batch being the common reference sample. The default is FALSE.

#' @return A list containing
#' \item{beta}{the estimated fixed-effects.}
#' \item{stat}{the parametric Wald statistics for testing non-zero fixed-effects. It is used in permutation tests.}
#' \item{Sigma}{the estimated error covariance matrix for the outcomes.}
#' \item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
#'  the variances for the first sample and the rest of samples within each cluster.}
#' \item{D}{the estimated covariance matrix for the random-effects.}
#' \item{phi}{the estimated parameters for the missing-data mechanism. Check the Details for the missing-data model. 
#' A zero value implies that parameter is ignored via the specification of miss_y and cov_miss.}
#' \item{loglikelihood}{the observed-data log-likelihood values.}
#' \item{iter}{the number of iterations for the EM algorithm when reaching the convergence.}
#' \item{AIC}{The Akaike information criterion (AIC) calculated for selecting the tuning parameter lambda of the graphical lasso penalty.}
#' 
#' @references Jiebiao Wang, Pei Wang, Donald Hedeker, and Lin S. Chen. A multivariate mixed-effects selection model framework for 
#' labelling-based proteomics data with non-ignorable missingness. (In preparation).
#' 
#' @export
#' @import lme4
#' @importFrom stats coef glm poisson pnorm
#' @examples
#' data(sim_dat)
#' 
#' # Covariates X common across outcomes with common coefficients
#'
#' fit0 = mvMISE_e(Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id)
#' 
#' \dontrun{
#' 
#' # In the example below, we showed how to estimate outcome-specific coefficients
#' # for a common covariate. The second column of sim_dat$X matrix is a
#' # common covariate. But it has different effects/coefficients
#' # on different outcomes.
#' 
#' nY = ncol(sim_dat$Y)
#' # stack X across outcomes
#' X_mat = sim_dat$X[rep(1:nrow(sim_dat$X), nY), ]
#' # Y_ind is the indicator matrix corresponding to different outcomes
#' Y_ind = kronecker(diag(nY), rep(1, nrow(sim_dat$Y)))
#' # generate outcome-specific covariates
#' cidx = 2 # the index for the covariate with outcome-specific coefficient
#' X_mat = cbind(1, X_mat[, cidx] * Y_ind) 
#' 
#' # X_mat is a matrix of 460 (92*5) by 6, the first column is intercept and the
#' # next 5 columns are covariate for each outcome
#'
#' fit1 = mvMISE_e(Y=sim_dat$Y, X=X_mat, id=sim_dat$id)
#' 
#' 
#' # A covariate only specific to the first outcome
#' 
#' X_mat1 = X_mat[, 1:2]
#' 
#' fit2 = mvMISE_e(Y=sim_dat$Y, X=X_mat1, id=sim_dat$id)
#' 
#' 
#' ## An example to allow missingness to depend on both a covariate and the outcome
#' 
#' fit3 = mvMISE_e(Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id, cov_miss = sim_dat$X[,2])
#' 
#' }


mvMISE_e = function(Y, X, id, Zidx = 1, maxIter = 100, tol = 0.001, lambda = 0.05, 
                    admm = TRUE, verbose = FALSE, cov_miss = NULL, miss_y = TRUE,
                    sigma_diff = FALSE) {
  
  N = length(unique(id))  # no. clusters
  nY = ncol(Y) # no. of outcomes
  id = as.character(id) # as names for list/vectors
  
  miss_cluster = NULL
  # number of samples within each cluster
  named_vec = rep(NA, N)
  names(named_vec) = unique(id)
  ni = named_vec
  for (i in unique(id)) {
    ni[i] = sum(id == i)
    # within each cluster, ordered by clusters
    miss_cluster = c(miss_cluster, apply(Y[id == i, , drop = FALSE], 2, function(x) mean(is.na(x)) == 1))
  }
  
  ## covariate mean for missing data mechanism
  if(!is.null(cov_miss)) {
    if(length(cov_miss) == length(Y)) {
      # mean covariate for each cluster
      cov_mean = tapply(cov_miss, rep(id, nY), mean)
      # mean covariate for each cluster and each outcome
      cov_mean_rep = as.vector(sapply(unique(id), function(i) tapply(cov_miss[rep(id, nY)==i], rep(1:nY, rep(ni[i], nY)), mean)))
    } else {
      cov_mean = tapply(cov_miss, id, mean)
      cov_mean_rep = rep(cov_mean, rep(nY, length(cov_mean)))
    }
  } else {
    cov_mean = rep(0, length(unique(id)))
    names(cov_mean) = unique(id)
  }
  
  # stack y, X, and id by outcome variables
  y = as.vector(Y)
  
  # if X is not stacked by outcomes, stack it
  if(nrow(X)/nrow(Y) == nY) X_mat = X else X_mat = X[rep(1:nrow(X), nY), ]
  
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
      # efficient inverse of Ri
      Ri_inv = kronecker(ginv(Sigma), diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i]))
      
      Vi = Zi %*% D %*% t(Zi) + Ri
      if (sum(!miss[id == i]) > 0) 
        Vi_obs_inv = ginv(Vi[!miss[id == i], !miss[id == i], drop = F]) else Vi_obs_inv = matrix(0, 0, 0)
      
      var_b_y[i] = ginv(t(Zi) %*% Ri_inv %*% Zi + ginv(D))
      
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
        log(1 - min(exp(phi[1] + phi[2] * mean(x) + phi[3] * cov_mean[i]), 1 - .Machine$double.eps)))) + 
          ifelse(sum(miss[id == i]) == 0, 0, phi[1] * sum(miss[id == i])/ni[i] + 
                   phi[3] * cov_mean[i] * sum(miss[id == i])/ni[i] + 
                   phi[2] * (sum(E_ymis_yobs) + phi[2] * sum(var_ymis_yobs)/ni[i]/2)/ni[i]))
      
      ## observed-data likelihood (observed y; r)
      ## overflow if log(det(Vi[!miss[id==i],!miss[id==i],drop=F]))
      obs_likelihood = obs_likelihood - 
        0.5 * (determinant(x = Vi[!miss[id == i], !miss[id == i], drop = F], logarithm = TRUE)$modulus + 
                 t(y[!miss & id == i] - Xio %*% beta) %*% Vi_obs_inv %*% 
                 (y[!miss & id == i] - Xio %*% beta)) + obs_phi
      SE = SE + t(Xio) %*% Vi_obs_inv %*% Xio
    }
    
    Sigma_inv = ginv(Sigma)
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
    
    # phi's: Poisson working model (Lumley et al., 2006, http://biostats.bepress.com/uwbiostat/paper293/)
    
    if (miss_y & is.null(cov_miss)) {
      mylog <- glm(miss_cluster ~ as.vector(t(Y_ik_mean)), family = poisson(link = "log"))
      # phi's index to be updated
      phi_idx = -3
    }
    if (!miss_y & !is.null(cov_miss)) {
      mylog <- glm(miss_cluster ~ cov_mean_rep, family = poisson(link = "log"))
      phi_idx = -2
    }
    if (miss_y & !is.null(cov_miss)) {
      mylog <- glm(miss_cluster ~ as.vector(t(Y_ik_mean)) + cov_mean_rep, family = poisson(link = "log"))
      phi_idx = 1:3
    }
    # make the algorithm stable
    if (miss_y | !is.null(cov_miss)) if (mylog$converged) phi[phi_idx] = as.numeric(coef(mylog))

    
    D_old = D
    D = (t(E_b_yobs) %*% E_b_yobs + sum(var_b_yobs + var_b_y))/N
    # make the algorithm stable
    if (D[1] > 1e+05) 
      D = D_old  
    beta = ginv(XRX) %*% XRE
    
    # random intercept model
    
    SS_sigma = matrix(0, nY, nY)
    SS_s = SS_s0 = 0
    Sigma_inv = ginv(Sigma)
    
    for (i in unique(id)) {
      Si_inv = diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i])
      Xi = X_mat[id == i, , drop = F]
      Zi = matrix(Xi[, Zidx, drop = F], ncol = nY)
      # efficient inverse of Ri
      Ri_inv = kronecker(ginv(Sigma), diag(1/c(sigma2_0, rep(sigma2, ni[i] - 1)), ni[i]))
      
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
      
      # identification issue for Kronecker product: set sigma2_0 as 1
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
  se = sqrt(diag(ginv(SE)))
  
  return(list(iter = iter, beta = beta, stat = beta/se, sigma2 = sigma2, 
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
  
  return(ginv(Theta))
}

