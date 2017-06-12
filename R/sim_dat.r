#' A Simulated Example data
#'
#' This simulated data list is for demonstration.
#' @docType data
#' @name sim_dat
#' @return A list containing
#' \item{Y}{a 92 by 5 outcome matrix, each row is a sample, and each 
#' column is an outcome variable, with potential missing values (NAs).}
#' \item{X}{a 92 by 2 covariate matrix, each row is a sample, and each column 
#' is a covariate with the first column being 1s for the intercept. In this example, we simulated the 
#' covariates to be common for all the outcomes and would estimate the common/averaged effects for all outcomes. 
#' If a covariate is specific for the k-th outcome, one may set all the values corresponding to the other outcomes to be zero. }
#' \item{id}{a vector of cluster/batch ID, matching with the rows of Y and X.}
#' @examples data(sim_dat)
NULL
