# frailtyEM 0.6.2
- big overhaul of the `control` argument and the `emfrail_control()` function

# frailtyEM 0.6.1 
- removed some old dependencies in the documentation and DESCRIPTION

# frailtyEM 0.6.0 (release)
- overall, numerous improvements compared to the previous release. Key new features include likelihood based confidence interval for the frailty parameter, more measures of dependence calculated with `summary()`, plots using `ggplot2`, and numerous bug fixes. 

# frailtyEM 0.5.13
- now the call is printed also when the summary is printed

# frailtyEM 0.5.12
- performance improvements. Now the likelihood-based confidence intervals should take less time as they know better where to look. 

# frailtyEM 0.5.11
- moved from `optimize` + `numDeriv` to `nlm`

# frailtyEM 0.5.11
- added a number of dependence measures that can be compared across distributions such as Kendall's tau, median concordance. 
- changed quite a lot in the structure of the summary object and the print method to make it more consistent and easier to develop in the future

# frailtyEM 0.5.10
- added score test for dependent censoring

# frailtyEM 0.5.9
- `ggplot_emfrail()`  added! Now the same plots (and more) can be done with the good looking `ggplot2` engine. 

# frailtyEM 0.5.8
- `summary.emfrail()` now has a new argument `print_opts` that is used in `print.emfrail_summary()`; if the output becomes too big, then some parts of the output may be ommitted

# frailtyEM 0.5.7

- The optimization now is regulated by search intervals described in the `emfrail_control()` and the `.control` argument. 
- The parametrization of the stable distribution has been changed, just removed the $1-$ in the beginning (why did I have that there again?)
- There are different intervals for the gamma/pvf and stable distributions. That's roughly because the stable chokes with small values of `theta`. This should be tuned somehow in the future. The problem lies in the M step where `agreg.fit` can't deal with large offset values.
- Likelihood confidence based intervals now do the correct thing when the estimate is close to the parameter space but not quite there
- Eliminated the fast fit for the inverse gaussian, this also seems to choke (the fast E step, can't figure out why), while the slow fit in C++ works fine...
- A slight update in documentation. 

TODO: 
- recover lost features in this update: measures of dependence in `summary.emfrail`, first of all
- bring back the fast fit for inverse gaussian or... who knows, maybe now
- document `emfrail_control` properly
- update vignette

# frailtyEM 0.5.6
Likelihood based confidence intervals are here! 

# frailtyEM 0.5.5
Removed the maximization by `optimx` and doing it with `optimize()`, since it's one dimensional. 
A hessian estimate is obtained from `numDeriv()`.

# frailtyEM 0.5.4
Minor bug fixes 

# frailtyEM 0.5.3
Some big changes in how the confidence intervals are constructed in predict.emfrail. Now - they are first constructed with the delta method for the log(cumulative hazard) and then exponentiated, so they do not have to be truncated at 0 or 1 any more. 

# frailtyEM 0.5.2
Further improved compatibility with CRAN policies and added a bunch of stuff in the examples in `\dontrun` statements (now they should be less than 5 seconds runtime)

# frailtyEM 0.5.1
Improved compatibility with R-devel 3.4.0. Registered C++ files to get rid of an R CMD check NOTE. Small modifications in the C++ file - for some reason a segfault started happening out of nowhere, think it's fixed now.

# frailtyEM 0.5.0
Added vignette, fixed small things for R CMD check
R CMD check: PASS, 0 warnings, 1 note / about new developer, that's ok.

# frailtyEM 0.4.9
Added the Commenges-Andersen test for heterogeneity. 
The test is implemented in a pretty non-efficient way, and it can be skipped with a proper `emfrail_control()` call, see `?emfrail_control`. Also there I added an option to *just* calculate the test, instead of doing anything else, and then just that is returned. A nice idea would be to implement this as a post-hoc calculation for `coxph` objects but that seems like another project atm.

R CMD check: PASS, 0 warnings, 0 notes.

# frailtyEM 0.4.8
Changed name to the more professional `frailtyEM`.
Added CI and SE for Kendall's tau with gamma

bugfixes: CI for tau with stable is now ok

# frailtoys 0.4.7
Added a `newdata` option for the `predict` method and for the `plot` methods. 
This can be used instead of `lp`, and basically calculates the corresponding linear predictor for the 
given covariate values.

bugfixes

# frailtoys 0.4.6
There is an option now to calculate the unadjusted SE or no SE at all

# frailtoys 0.4.5

* There are now plot methods available! Check out `?plot_emfrail`\
* Documentation updated accordingly





# frailtoys 0.4.3

* Added a `NEWS.md` file to track changes to the package.



