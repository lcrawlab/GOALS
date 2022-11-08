# Demonstrations and Tutorials on Running GOALS

In the file entitled "Global_Tutorial.R", we walk through the implementation of computing/utilizing distributional centrality measures. Specifically, we describe in detail: (1) how to compute a covariance matrix using the Gaussian kernel function; (2) how to fit a standard Bayesian Gaussian process (GP) regression model; and (3) calculate the delta GOALS matrix. This script is based on a simple (and small) genetics example where we simulate genotype data for some n individuals with p measured genetic variants. We randomly assume that three of the predictors j* = {23, 24, 25} are causal and have true association with the generated (continuous) phenotype y. We then assume that the j* predictor variables explain a fixed H2% (phenotypic variance explained; PVE) of the total variance in the response V(y). This parameter H2 can alternatively be described as a factor controlling the signal-to-noise ratio. The parameter rho represents the proportion of H2 that is contributed by additive effects versus interaction effects. Namely, the additive effects make up rho%, while the pairwise interactions make up the remaining (1 − rho)%.


In the file entitled "Local_Tutorial.R", we walk through the implementation of computing/utilizing distributional centrality measures. Specifically, we describe in detail: (1) how to compute a covariance matrix using the Gaussian kernel function; (2) how to fit a standard Bayesian Gaussian process (GP) regression model; and (3) calculate the delta GOALS matrix. This script is based on a simple (and small) genetics example where we simulate genotype data for some n individuals with p measured genetic variants. For half of individuals, we randomly assume that three of the predictors j* = {23, 24, 25} are causal and have true association with the generated (continuous) phenotype y. For the other half of individuals, we randomly assume that three of the predictors j* = {22, 23, 24, 25} are causal and have true association with the generated (continuous) phenotype y. We then assume that the j* predictor variables explain a fixed H2% (phenotypic variance explained; PVE) of the total variance in the response V(y). This parameter H2 can alternatively be described as a factor controlling the signal-to-noise ratio. The parameter rho represents the proportion of H2 that is contributed by additive effects versus interaction effects. Namely, the additive effects make up rho%, while the pairwise interactions make up the remaining (1 − rho)%.

In the "Power_Comparisons.R" file, we demonstrate the power of distributional centrality via RATE measures. Here, we focus on simulated genotype data for n = 2000 individuals with p = 10000 measured genetic variants. Next, we compare our approach to: (1) RATE, a KL-divergence based variable selection model also built on Gaussian processes; (2) L1- regularized lasso regression; (3) the combined regularization utilized by the elastic net; (4) a genome scan with individual single nucleotide polymorphisms (SNPs) fit via a univariate linear model (SCANONE). This script is based on the simulations from Crawford et al. (2018). We assess the association mapping ability of each method with (i) different signal-to-noise ratios in H2, (ii) varying levels of additive and interaction effects in rho, and (iii) with data affected by population stratification.
R Packages Required for GOALS Tutorials

The GOALS tutorial and examples require the additional installation of the following R libraries:

[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)

[varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)

[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

Once again, the easiest method to install these packages is with the following example command entered in an R shell:

  install.packages("glmnet", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).
