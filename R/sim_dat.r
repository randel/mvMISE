#' Example data
#'
#' This simulated data list is for demonstration.
#' @docType data
#' @name sim_dat
#' @return A list containing
#' \item{Y}{an outcome matrix, each row is an observation, each column is an outcome variable, with potential missing values (NAs).}
#' \item{X}{a covariates matrix, each row is an observation, each column is a covariate. Now covariates are assumed to be common for outcomes.}
#' \item{id}{a vector for cluster/grouping index, matching with the rows of Y and X.}
#' @examples data(sim_dat)
NULL
