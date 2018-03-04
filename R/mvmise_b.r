#' A multivariate mixed-effects selection model with correlated outcome-specific random intercepts
#' 
#' This function fits a multivariate mixed-effects selection model with correlated outcome-specific
#' random intercepts allowing potential ignorable or non-ignorable missing values in the outcome.
#' Here an outcome refers to a response variable, for example, a genomic feature. The proposed model and function jointly analyze multiple outcomes/features.
#' 
#' The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
#' is a multivariate mixed-effects model, with correlations among multivariate outcomes modeled via correlated outcome-specific random intercepts with 
#' a factor-analytic structure
#' \deqn{\mathbf{y}_{i} =  \mathbf{X}_{i}\boldsymbol{\beta} + \left(\mathbf{I}_{K}\otimes\mathbf{1}_{n_{i}}\right) \boldsymbol{\tau}b_{i}+\mathbf{e}_{i},}
#' where \eqn{i} denotes a cluster/batch, \eqn{n_{i}} is the number of samples/observations within each cluster,
#'  \eqn{\boldsymbol{\tau}} is a \eqn{K\times 1} vector for the outcome-specific variance components corresponding to 
#' the random effect \eqn{b_i} (a standard normal random variable), and \eqn{K} is the number of outcomes. 
#' By default, a matrix with each column as an indicator for each outcome is generated and is used as the random-effect design matrix (\eqn{\mathbf{I}_{K}\otimes\mathbf{1}_{n_{i}}}), 
#' and the model will estimate the outcome-specific random intercepts.
#' The factor-analytic structure assumes the outcome-specific random intercepts are identically correlated and this model 
#' is often used to capture the highly structured experimental or biological correlations among naturally related outcomes.
#' For example, the correlation among multiple phosphopeptides (i.e. phosphorylated segments) of a same protein.
#' The model assumes that the random effects are derived from a latent variable \eqn{b_i} with a loading vector \eqn{\boldsymbol{\tau}}.
#' With this model specification, only \eqn{K} parameters instead of \eqn{K(K+1)/2} are needed 
#' in the estimation for the covariance matrix of random-effects, and as such that greatly facilitates the computation. 
#' 
#' The missing-data model can be written as
#' \deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}/n_{i}\cdot \mathbf{1}^{'}\mathbf{y}_{ik}  + 
#' \phi_{2}/n_{i}\cdot \mathbf{1}^{'}\mathbf{c}_{i} \right),}
#' where \eqn{r_{ik}} is the missing indicator for the k-th outcome in the i-th cluster. If \eqn{r_{ik}=1}, the values of the k-th outcome in the i-th cluster  
#' \eqn{\mathbf{y}_{ik}} are missing altogether.
#' The estimation is implemented via an EM algorithm. Parameters in the missing-data models can be specified via the arguments miss_y and cov_miss. If miss_y 
#' = TURE, the missingness depends on the outcome values. 
#' If cov_miss is specified, the missingness can (additionally) depend on the specified covariate (cov_miss).
#' 
#' The model also works for fully observed data if miss_y = FALSE and cov_miss = NULL. It would also work for a univariate outcome with potential missing values, if the outcome Y is a matrix
#' with one column.
#' 
#' @param Y an outcome matrix. Each row is a sample, and each column is an outcome variable, with potential missing values (NAs).
#' @param X a covariate matrix. Each row is a sample, and each column is a covariate. The covariates can be common among all of the outcomes (e.g., age, gender) or outcome-specific.
#'    If a covariate is specific for the k-th outcome, one may set all the values corresponding to the other outcomes to be zero. 
#'    If X is common across outcomes, the row number of X equals 
#'    the row number of Y. Otherwise, if X is outcome-specific, the row number of X equals the number of elements in Y, i.e., outcome-specific X is stacked across outcomes within
#'    each cluster. See the Examples for demonstration.
#' @param id a vector of cluster/batch index, matching with the rows of Y, and X if it is not outcome specific.
#' @param maxIter the maximum number of iterations for the EM algorithm.
#' @param tol the tolerance level for the relative change in the observed-data log-likelihood function.
#' @param verbose logical. If TRUE, the iteration history of each step of the EM algorithm will be printed. The default is FALSE.
#' @param miss_y logical. If TRUE, the missingness depends on the outcome Y (see the Details). The default is TRUE.
#'      This outcome-dependent missing data pattern was motivated by and was observed in the mass-spectrometry-based quantitative proteomics data.  
#' @param cov_miss the covariate that can be used in the missing-data model. If it is NULL, 
#' the missingness is assumed to be independent of the covariates. 
#' Check the Details for the missing-data model.
#' If it is specified and the covariate is not outcome specific, its length equals the length of id. If it is outcome specific, the outcome-specific covariate is stacked across outcomes within
#' each cluster.
#' @param sigma_diff logical. If TRUE, the sample error variance of the first sample in each cluster/batch is different from that for the rest of samples within the same cluster/batch.
#' This option is designed and used when analyzing batch-processed proteomics data with the first sample in each cluster/batch being the common reference sample. The default is FALSE.
#' 
#' @return A list containing
#' \item{beta}{the estimated fixed-effects.}
#' \item{var}{the variance-covariance matrix of the estimated fixed effects. With the fixed effects and their covariance matrix estimates, 
#' one can obtain the Wald-statistics for testing fixed-effects beta/sqrt(diag(var)).}
#' \item{pval}{the parametric p-values for testing non-zero fixed-effects. It is obtained as the two-sided p-value based on the Wald statistics of beta/sqrt(diag(var)).}
#' \item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
#'  the variances for the first sample and for the rest of samples within each cluster.}
#' \item{tau}{the estimated variance components for the outcome-specific factor-analytic random-effects.}
#' \item{phi}{the estimated parameters for the missing-data mechanism. Check the Details for the missing-data model.
#'  A zero estimate implies that the parameter is ignored via the specification of miss_y and/or cov_miss.}
#' \item{loglikelihood}{the observed-data log-likelihood values.}
#' \item{iter}{the number of iterations for the EM algorithm when reaching the convergence.}
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
#' fit0 = mvMISE_b(Y=sim_dat$Y, X=sim_dat$X, id=sim_dat$id)
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
#' fit1 = mvMISE_b(Y=sim_dat$Y, X=X_mat, id=sim_dat$id)
#' 
#' 
#' # A covariate only specific to the first outcome
#' 
#' X_mat1 = X_mat[, 1:2]
#' 
#' fit2 = mvMISE_b(Y=sim_dat$Y, X=X_mat1, id=sim_dat$id)
#' 
#' 
#' ## An example that allows missingness depending on both a covariate and the outcome
#' 
#' fit3 = mvMISE_e(Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id, cov_miss = sim_dat$X[,2])
#' 
#' }


mvMISE_b = function(Y, X,  id, maxIter = 100, tol = 0.001, verbose = FALSE, cov_miss = NULL, miss_y = TRUE,
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
      
      Zi = kronecker(diag(nY), rep(1, ni[i]))
      
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
        log(1 - min(exp(phi[1] + phi[2] * mean(x) + phi[3] * cov_mean[i]), 1 - .Machine$double.eps)))) + 
          ifelse(sum(miss[id == i]) == 0, 0, phi[1] * sum(miss[id == i])/ni[i] + 
                   phi[3] * cov_mean[i] * sum(miss[id == i])/ni[i] + 
                   phi[2] * (sum(E_ymis_yobs) + phi[2] * sum(var_ymis_yobs)/ni[i]/2)/ni[i]))
      
      ## observed-data likelihood (observed y; r)
      obs_likelihood = obs_likelihood - 0.5 * (ifelse(ncol(Vi_obs) == 0, 0, determinant(Vi_obs, logarithm = TRUE)$modulus) + 
                                                 t(y[!miss & id == i] - Xio %*% beta) %*% Vi_obs_inv %*% 
                                                 (y[!miss & id == i] - Xio %*% beta)) + obs_phi
      
      SE = SE + t(Xio) %*% Vi_obs_inv %*% Xio
    }
    
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
    
    beta = ginv(XRX) %*% XRE
    tau = ginv(ZRZ) %*% ZRE
    
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
  se = sqrt(diag(ginv(SE)))
  
  return(list(iter = iter, beta = beta, var = ginv(SE), pval = 2*pnorm(-abs(beta/se)), sigma2 = sigma2, 
              tau = tau, phi = phi, loglikelihood = likelihood[length(likelihood)]))
}

