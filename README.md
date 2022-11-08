# Revisiting the Effect Size Analog for Nonlinear Regression: A Simpler Approach for Interpretability (GOALS)

When calculating the influence a varaible has on a model, there is a tension between accounting for all linear and nonlinear effects, interpretability and transparency of the algorithm, and computatational complexity and feasibility.

## R Packages for GOALS

The GOALS software requires the installation of the following R libraries:

[BAKR](https://github.com/lorinanthony/BAKR) (via GitHub) (Note that the latest BAKR functions are in the '''Software''' directory for GOALS.

[corpcor](https://cran.r-project.org/web/packages/corpcor/index.html)

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[MASS](https://cran.r-project.org/web/packages/MASS/index.html)

[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)

[RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)

[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)

Unless stated otherwise, the easiest method to install many of these packages is with the following example command entered in an R shell:

    install.packages("corpcor", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

## C++ Packages for GOALS

The code in this repository assumes that basic C++ functions and applications are already set up on the running personal computer or cluster. If not, the functions and necessary Rcpp packages to build nonlinear covariance matrices (e.g. BAKR) and fit a GP regression model will not work properly. A simple option is to use [gcc](https://gcc.gnu.org/). macOS users may use this collection by installing the [Homebrew package manager](http://brew.sh/index.html) and then typing the following into the terminal:

    brew install gcc

For macOS users, the Xcode Command Line Tools include a GCC compiler. Instructions on how to install Xcode may be found [here](http://railsapps.github.io/xcode-command-line-tools.html). For extra tips on how to run C++ on macOS, please visit [here](http://seananderson.ca/2013/11/18/rcpp-mavericks.html). For tips on how to avoid errors dealing with "-lgfortran" or "-lquadmath", please visit [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

## Demonstration and Tutorials for GOALS

The `Software` folder contains all necessary software to run GOALS. Note that '''RATE.R''' is only necessary for the '''PowerComparisons.R''' tutorial.

The `Tutorial` folder contains a tutorial on global variable selection, a tutorial on local variable selection, and a Power Comparison of GOALS to other variable selection methods. Here, we consider a simple (and small) genetics example where we simulate genotype data for n individuals with p measured genetic variants. We then randomly select a small number of these predictor variables to be causal and have true association with the generated (continuous) phenotype. These scripts are meant to illustrate proof of concepts and specifically demonstrate (1) the calculation of the Gaussian Kernel, (2) the precalculations and the calculation of the delta GOALS matrix Additional packages needed for the tutorials are in the README.md file in that folder. 

## Relevant Citations

Goals paper, RATE paper, BAKR?

## Questions and Feedback

Please send any questions or feedback to the corresponding authors Emily Winn-Nuñez (emily_winn@brown.edu) or Lorin Crawford (lorin_crawford@brown.edu).
