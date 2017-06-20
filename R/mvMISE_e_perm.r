#' A function to obtain permutation-based p-values for fixed effects estimates in mvMISE_e
#' 
#' This function calls mvMISE_e multiple times by permuting the row index (observations) of the covariate matrix X.
#' It may take a long time to permute high-dimensional outcomes, but can be run in parallel using multiple nodes.
#' 
#' @param nperm the number of permutations.
#' @param nnodes the number of nodes that will be used in parallel for permutations.
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
#' @param lambda the tuning parameter for the graphical lasso penalty of the error precision matrix. It can be selected by AIC (an output).
#' @param tol the tolerance level for the relative change in the observed-data log-likelihood function.
#' @param miss_y logical. If TRUE, the missingness depends on the outcome Y (see the Details). The default is TRUE.
#'      This outcome-dependent missing data pattern was motivated by and was observed in the mass-spectrometry-based quantitative proteomics data.  
#' @param cov_miss the covariate that can be used in the missing-data model. If it is NULL, 
#'    the missingness is assumed to be independent of the covariates.
#'   Check the Details for the missing-data model.
#'   If it is specified and the covariate is not outcome specific, its length equals the length of id. If it is outcome specific, the outcome-specific covariate is stacked across outcomes within
#'   each cluster.
#' @param sigma_diff logical. If TRUE, the sample error variance of the first sample is different from that for the rest of samples within each cluster.
#' This option is designed and used when analyzing batch-processed proteomics data with the first sample in each cluster/batch being the common reference sample. The default is FALSE.
#' 
#' @return The permutation based p-values for testing if fixed-effects (excluding the intercept) are zeros.
#' 
#' @references Jiebiao Wang, Pei Wang, Donald Hedeker, and Lin S. Chen. A multivariate mixed-effects selection model framework for 
#' labelling-based proteomics data with non-ignorable missingness. (In preparation).
#' 
#' @importFrom parallel makeCluster stopCluster parLapply
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' data(sim_dat)
#' 
#' pval_perm = mvMISE_e_perm(nperm=100, Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id)
#' }


mvMISE_e_perm = function(nperm = 100, nnodes = 2, Y, X, id, Zidx = 1, maxIter = 100, tol = 0.001, lambda = 0.05, 
                         cov_miss = NULL, miss_y = TRUE, sigma_diff = FALSE) {
  
  cl = makeCluster(nnodes)
  
  # test statistics
  
  stat = mvMISE_e(Y = Y, X = X, id = id, Zidx = Zidx, 
                  maxIter = maxIter, tol = tol, lambda = lambda,
                  cov_miss = cov_miss, miss_y = miss_y,
                  sigma_diff = sigma_diff)$stat
  
  # permute the row (sample index) of X
  # parLapply has a default X argument and does not allow a duplicated one
  
  mvMISE_e1 = function(...) mvMISE_e(Y = Y, X = X[sample(nrow(X)),], id = id, Zidx = Zidx, 
                   maxIter = maxIter, tol = tol, lambda = lambda,
                   cov_miss = cov_miss, miss_y = miss_y,
                   sigma_diff = sigma_diff, admm = TRUE, verbose = FALSE)
  
  fit0 = parLapply(cl, 1:nperm, mvMISE_e1)
  
  # parameters x nperm
  
  stat0 = sapply(fit0, function(x) x$stat)
  
  makeCluster(cl)
  
  pval_perm = sapply(1:length(stat), function(x) mean(abs(stat0[x,]) >= abs(stat[x])))
  
  return(pval_perm[-1])
}
