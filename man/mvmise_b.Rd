% Please edit documentation in R/mvmise_b.r
\name{mvmise_b}
\alias{mvmise_b}
\title{Multivariate mixed-effects selection model with correlated outcome-specific random intercepts}
\usage{
mvmise_b(Y, X, Z = NULL, id, maxIter = 100, tol = 0.001, verbose = FALSE, 
    specific_eff = FALSE, miss_mechanism = "y", sigma_diff = FALSE)
}
\arguments{
\item{Y}{an outcome matrix, each row is an observation, each column is an outcome variable, with potential missing values (NAs).}

\item{X}{a covariates matrix, each row is an observation, each column is a covariate. Now covariates are assumed to be common for outcomes.}

\item{Z}{a design matrix for random effects, each row is an observation, each column is a random effect.
If it is NULL (the default), a matrix with each column as an indicator for each outcome is generated.}

\item{id}{a vector for cluster/grouping index, matching with the rows of Y, X, Z (if specified).}

\item{maxIter}{maximum number of iterations for the EM algorithm.}

\item{tol}{tolerance level for the relative change in the observed-data log-likelihood function.}

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
\item{sigma2}{the estimated sample error variance(s). If sigma_diff is TRUE, it returns a vector of two elements,
 the variances for the first sample and the rest of samples within each cluster.}
\item{tau}{the estimated variance components for the outcome-specific factor-analytic random effects.}
\item{phi}{the estimated parameters for the missing-data mechanism. The missing probability is modelled as exp(phi0 +
phi1*mean(y) + phi2*mean(x)). A zero value implies that parameter is ignored via the specification of miss_mechanism.}
\item{loglikelihood}{the observed-data log-likelihood values.}
\item{iter}{the number of iterations for the EM algorithm.}
}
\description{
This function fits a multivariate mixed-effects selection model with potential missing values in the outcome
 and correlated outcome-specific random intercepts.
}
\details{
The multivariate mixed-effects selection model consists of two components, the outcome model and the missing-data model. Here the outcome model 
is a multivariate mixed-effects model, with correlations among multivariate outcomes modelled via outcome-specific random intercepts with 
a factor-analytic structure
\deqn{\mathbf{y}_{i} =  \mathbf{X}_{i}\boldsymbol{\alpha}+\mathbf{Z}_{i}\boldsymbol{\tau}b_{i}+\mathbf{e}_{i},}
where \eqn{i} denotes a cluster, \eqn{\boldsymbol{\tau}} is a \eqn{K\times 1} vector for the outcome-specific variance components corresponding to 
the random effect \eqn{b_i} (a standard normal random variable), and \eqn{K} is the number of outcomes. 
The factor-analytic structure is used to facilitate the computation.
It assumes that the random effects 
are derived from a latent variable \eqn{b_i} with a loading vector \eqn{\boldsymbol{\tau}}.
In this way, only \eqn{K} rather parameters are needed 
in the estimation for the covariance matrix of random effects. The last fixed effect in \eqn{\boldsymbol{\alpha}} can be outcome-specific, if specific_eff
is specified as TRUE.

The missing-data model can be written as
\deqn{\textrm{Pr}\left(r_{ik}=1|\mathbf{y}_{ik}\right)= \mathrm{exp}\left(\phi_{0} + \phi_{1}\mathbf{1}_{n_{i}}^{'}\mathbf{y}_{ik}  + 
\phi_{2}\mathbf{1}_{n_{i}}^{'}\mathbf{x}_{i} \right),}
where \eqn{r_{ik}} is the missing indicator for the kth outcome in the ith cluster. If missing, the kth outcome in the ith cluster  \eqn{\mathbf{y}_{ik}} 
is missing altogether.
The estimation is implemented via an EM algorithm. Parameters in the missing-data models can be specified via the argument miss_mechanism. If miss_mechanism 
= "y" or "yx", i.e., the missingness depends on the outcome, the missing-data mechanism is missing not at random (MNAR), 
otherwise it is missing at random (MAR).

It works for fully observed data if miss_mechanism = "none". It also works for univariate outcome with potential missing values, if the outcome Y is a matrix
with one column.
}
\examples{
data(sim_dat)

fit0 = mvmise_b(Y = sim_dat$Y, X = sim_dat$X, id = sim_dat$id)
}
\references{
Jiebiao Wang, Pei Wang, Donald Hedeker, and Lin S. Chen. A multivariate mixed-effects selection model framework for 
labelling-based proteomics data with non-ignorable missingness. (In preparation).
}
