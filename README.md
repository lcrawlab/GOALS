# GOALS: A Simple Approach for Local and Global Variable Importance in Nonlinear Regression Models

The ability to interpret machine learning models has become increasingly important as their usage  in data science continues to rise. Most current interpretability methods are optimized to work on either (i) a global scale, where the goal is to rank features based on their contributions to overall variation in an observed population, or (ii) the local level, which aims to detail on how important a feature is to a particular individual in the dataset. In this work, we present the "GlObal And Local Score" (GOALS) operator: a simple _post hoc_ approach to simultaneously assess local and global feature variable importance in nonlinear models. Motivated by problems in statistical genetics, we demonstrate our approach using Gaussian process regression where understanding how genetic markers affect trait architecture both among individuals and across populations is of high interest. With detailed simulations and real data analyses, we illustrate the flexible and efficient utility of GOALS over state-of-the-art variable importance strategies.

## R Packages for GOALS and Tutorials

The GOALS software requires the installation of the following R libraries:

[adegenet](https://cran.r-project.org/web/packages/adegenet/index.html)

[BAKR](https://github.com/lorinanthony/BAKR) (via GitHub)

[corpcor](https://cran.r-project.org/web/packages/corpcor/index.html)

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[MASS](https://cran.r-project.org/web/packages/MASS/index.html)

[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)

[RATE](https://github.com/lorinanthony/RATE) (via GitHub)

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)

[RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)

[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)

[svd](https://cran.r-project.org/web/packages/svd/index.html)

Note that the latest BAKR and RATE functions are also included in the `Software` directory for this repo. Unless stated otherwise, the easiest method to install many of these packages is with the following example command entered in an R shell:

    install.packages("corpcor", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

## C++ Packages for GOALS and Tutorials

The code in this repository assumes that basic C++ functions and applications are already set up on the running personal computer or cluster. If not, the functions and necessary Rcpp packages to build nonlinear covariance matrices (e.g., BAKR and RATE) and fit a GP regression model will not work properly. A simple option is to use [gcc](https://gcc.gnu.org/). macOS users may use this collection by installing the [Homebrew package manager](http://brew.sh/index.html) and then typing the following into the terminal:

    brew install gcc

For macOS users, the Xcode Command Line Tools include a GCC compiler. Instructions on how to install Xcode may be found [here](http://railsapps.github.io/xcode-command-line-tools.html). For extra tips on how to run C++ on macOS, please visit [here](http://seananderson.ca/2013/11/18/rcpp-mavericks.html). For tips on how to avoid errors dealing with "-lgfortran" or "-lquadmath", please visit [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

## Demonstration and Tutorials for GOALS

The `Software` folder contains all necessary software to run GOALS. Note that `RATE.R` is only necessary for the `PowerComparisons.R` tutorial.

The `Tutorial` folder contains a tutorial on global variable selection, a tutorial on local variable selection, and a Power Comparison of GOALS to other variable selection methods. Here, we consider a simple (and small) genetics example where we simulate genotype data for n individuals with p measured genetic variants. We then randomly select a small number of these predictor variables to be causal and have true association with the generated (continuous) phenotype. These scripts are meant to illustrate proof of concepts and specifically demonstrate (1) the calculation of the Gaussian Kernel, (2) the precalculations and the calculation of the delta GOALS matrix Additional packages needed for the tutorials are in the README.md file in that folder. 

## Relevant Citations

E.T. Winn-Nuñez, M. Griffin, and L. Crawford (2022). A simple approach for local and global variable importance in nonlinear regression models. _arXiv_.

J. Ish-Horowicz*, D. Udwin*, S.R. Flaxman, S.L. Filippi, and L. Crawford (2019). Interpreting deep neural networks through variable importance. _arXiv_. 1901.09839.

L. Crawford, S.R. Flaxman, D.E. Runcie, and M. West (2019). Variable prioritization in nonlinear black box methods: a genetic association case study. _Annals of Applied Statistics_. **13**(2): 958-989.

L. Crawford, K.C. Wood, X. Zhou, and S. Mukherjee (2018). Bayesian approximate kernel regression with variable selection. _Journal of the American Statistical Association_. **113**(524): 1710-1721.

## Questions and Feedback

Please send any questions or feedback to the corresponding authors [Emily Winn-Nuñez](mailto:emily_winn@brown.edu) or [Lorin Crawford](mailto:lcrawford@microsoft.com).

We appreciate any feedback you may have with our repository and instructions.
