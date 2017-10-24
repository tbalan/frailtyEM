
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/frailtyEM)](https://cran.r-project.org/package=frailtyEM)

This is an R package for fitting semiparametric shared frailty models with the EM algorithm. You can check the "issues" section to see about known issues. For the gamma frailty model, the results are identical with those from the `survival` pacakage, although `frailtyEM` provides a more readable output, including confidence intervals for the frailty variance. Other supported distributions include the PVF, compound Poisson, inverse Gaussian, positive stable. Univariate and multivariate data with left truncation are supported, including recurrent events data in Andersen-Gill formulation.

The stable version may be installed from `CRAN`:

``` r
install.packages("frailtyEM")
```

and the development version from `GitHub`:

``` r
devtools::install_github("tbalan/frailtyEM")
```

The bulk of the documentation of the package can be found in the vignette. If the package is installed from `GitHub`, then the vignette is installed if the pacakge is installed like this:

``` r
devtools::install_github("tbalan/frailtyEM", build_vignettes = TRUE)
```

### Functions

The main fitting function is `emfrail()`, which in general returns an `emfrail()` object. Several plots can be produced from this objects, via the `plot.emfrail()` and `autoplot.emfrail()` methods (the latter using `ggplot2`).

Another useful tool is the Commenges-Andersen score test for heterogeneity. This test does not require estimating the shared frailty model. The `ca_test()` function may be used in conjunction with a `coxph` object to calculate this.
