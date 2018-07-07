# mvMISE

### A General Framework for Multivariate Mixed-Effects Selection Models with Potential Missing Data

This `R` package offers a general framework of multivariate mixed-effects models for the joint analysis of multiple correlated outcomes with clustered data structures and potential missingness proposed by Wang et al. (2018) <doi:10.1093/biostatistics/kxy022>. The missingness of outcome values may depend on the values themselves (missing not at random and non-ignorable), or may depend on only the covariates (missing at random and ignorable), or both. This package provides functions for two models: 
1) mvMISE_b allows correlated outcome-specific random intercepts with a factor-analytic structure;
2) mvMISE_e allows the correlated outcome-specific error terms with a graphical lasso penalty on the error precision matrix. 

Both functions are motivated by the multivariate data analysis on data with clustered structures from labelling-based quantitative proteomic studies. These models and functions can also be applied to univariate and multivariate analyses of clustered data with balanced or unbalanced design and no missingness.


### Installation
- For the stable version from [CRAN](https://cran.r-project.org/web/packages/mvMISE/index.html):
```r
install.packages('mvMISE')
```
- For the development version (requiring the `devtools` package):
```r
devtools::install_github('randel/mvMISE')
```

### Reference
Wang, J., Wang, P., Hedeker, D., & Chen, L. S. (2018). Using multivariate mixed-effects selection models for analyzing batch-processed proteomics data with non-ignorable missingness. *Biostatistics*.
