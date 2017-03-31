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



