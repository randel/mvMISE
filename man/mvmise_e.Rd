% Please edit documentation in R/mvmise_e.r
\name{mvmise_e}
\alias{mvmise_e}
\title{Multivariate mixed-effects selection model with correlated outcome-specific error terms}
\usage{
mvmise_e(Y, X, Zidx = 1, id, maxIter = 100, tol = 0.001, lambda = 0.05, admm = TRUE, 
    verbose = FALSE, specific_eff = FALSE, miss_mechanism = "y", sigma_diff = FALSE)
}
\arguments{
\item{Y}{an outcome matrix, each row is an observation, each column is an outcome variable, with potential missing values (NAs).}

\item{X}{a covariates matrix, each row is an observation, each column is a covariate. Now covariates are assumed to be common for outcomes.}

\item{Zidx}{column indexes of matrix X used as the design matrix of random effects. The default is 1, i.e., a random intercept is included 
if the first column of X is a vector of 1s.}

\item{id}{a vector for cluster/grouping index, matching with the rows of Y and X.}

\item{maxIter}{maximum number of iterations for the EM algorithm.}

\item{tol}{tolerance level for the relative change in the observed-data log-likelihood function.}

\item{lambda}{tuning parameter for the graphical lasso penalty of the error precision matrix. It can be selected by AIC (an output).}

\item{admm}{logical. If TRUE (the default), the alternating direction method of multipliers (ADMM) is 
used to estimate the error precision matrix with a graphical lasso penalty. This works for multivariate outcomes. 
For an univariate outcome, it should be set as FALSE.}

\item{verbose}{logical. If TRUE, the iteration history of each step of the EM algorithm will be printed. The default is FALSE.}

\item{specific_eff}{logical. If TRUE, outcome-specific fixed-effects are estimated for the last covariate in X. The default is FALSE.}

\item{miss_mechanism}{one of "y" (the default), "x", "yx", and "none", indicating the missingness of outcome k in cluster i 
depends on the mean of the outcome, the mean of the covariate of interest, both, or none. The missing probability is modelled as exp(phi0 +
phi1*mean(y) + phi2*mean(x)). If there is no missing values in Y, it should be set as "none".}

\item{sigma_diff}{logical. If TRUE, the sample error variance of the first sample is different from that for the rest of samples within each cluster.
This is the case for the reference sample in the iTRAQ proteomics data. The default is FALSE.}
}
\value{
A list containing
\item{beta}{the estimated fixed effects.}
\item{se}{the standard errors for the estimated fixed effects.}
\item{Sigma}{the estimated error covariance matrix for the outcomes.}
\item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
 the variances for the first sample and the rest of samples within each cluster.}
\item{D}{the estimated covariance matrix for the random effects.}
\item{phi}{the estimated parameters for the missing-data mechanism. The missing probability is modelled as exp(phi0 +
phi1*mean(y) + phi2*mean(x)). A zero value implies that parameter is ignored via the specification of miss_mechanism.}
\item{loglikelihood}{the observed-data log-likelihood values.}
\item{iter}{the number of iterations for the EM algorithm.}
\item{AIC}{The Akaike information criterion (AIC) calculated for selecting the tuning parameter lambda.}
}
\description{
This function fits a multivariate mixed-effects selection model with potential missing values in the outcome and correlated outcome-specific error terms.
It can shrink the error precision matrix with a graphical lasso penalty 
for high-dimensional outcomes.
}
\details{
The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
is a multivariate mixed-effects model, with correlations among multivariate outcomes modelled via outcome-specific error terms. For the ith cluster,
the outcome \eqn{\mathbf{Y}_{i}} is a matrix of \eqn{n_i} samples (rows) and \eqn{K} outcomes (columns). 
Let \eqn{\mathbf{y}_{i} = \mathrm{vec}\left( \mathbf{Y}_{i} \right)}. 
The outcome vector \eqn{\mathbf{y}_{i}} can be modelled as
\deqn{\mathbf{y}_{i}  = \mathbf{X}_{i}\boldsymbol{\alpha}+\mathbf{Z}_{i}\mathbf{b}_{i}+\mathbf{e}_{i},}
where the random effects (\eqn{\mathbf{b}_{i}}) follow a normal distribution \eqn{\mathbf{b}_{i}\sim N(\mathbf{0},\mathbf{D})}; 
and the error term \eqn{\mathbf{e}_{i}=\mathrm{vec}\left(\mathbf{E}_{i}\right) \sim N(\mathbf{0},\boldsymbol{\Sigma}\otimes\mathbf{S}_{i})}.  
The matrix \eqn{\mathbf{S}_{i}} is an \eqn{n_i\times n_i} diagonal matrix with diagonal elements corresponding to the error variances of the \eqn{n_i} samples 
within the ith cluster. 
The variances for the first and other samples can be different if sigma_diff = TRUE. 
The matrix \eqn{\boldsymbol{\Sigma}} captures the error (or unexplained) covariances among \eqn{K} outcomes. 
To facilitate the computation for high-dimensional outcomes, the off-diagonal elements of the inverse of \eqn{\boldsymbol{\Sigma}} can be shrinked
by a graphical lasso penalty. If admm = TRUE (the default), the alternating direction method of multipliers (ADMM) is used to estimate \eqn{\boldsymbol{\Sigma}}.
The last fixed effect in \eqn{\boldsymbol{\alpha}} can be outcome-specific, if specific_eff
is specified as TRUE.

The missing-data model can be written as
\deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}\mathbf{1}_{n_{i}}^{'}\mathbf{y}_{ik}  + 
\phi_{2}\mathbf{1}_{n_{i}}^{'}\mathbf{x}_{i} \right),}
where \eqn{r_{ik}} is the missing indicator for the kth outcome in the ith cluster. If missing, the kth outcome in the ith cluster \eqn{\mathbf{y}_{ik}} 
is missing altogether.
The estimation is implemented via an EM algorithm. Parameters in the missing-data models can be specified via the argument miss_mechanism. If miss_mechanism 
= "y" or "yx", i.e., the missingness depends on the outcome, the missing-data mechanism is missing not at random (MNAR), 
otherwise it is missing at random (MAR).

It works for fully observed data if miss_mechanism = "none". It also works for univariate outcome with potential missing values, if the outcome Y is a matrix
with one column.
}
\examples{
data(sim_dat)

fit0 = mvmise_e(Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id)
}
\references{
Jiebiao Wang, Pei Wang, Donald Hedeker, and Lin S. Chen. A multivariate mixed-effects selection model framework for 
labelling-based proteomics data with non-ignorable missingness. (In preparation).
}
