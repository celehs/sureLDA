## Overview

Surrogate-guided ensemble Latent Dirichlet Allocation (sureLDA) is a label-free multidimensional phenotyping method. It first uses the PheNorm or MAP algorithm to initialize probabilities based on two surrogate features for each target disease, and then leverages these probabilities to guide the LDA topic model to generate phenotype-specific topics. Finally, it combines phenotype-feature counts with surrogates via clustering ensemble to yield final phenotype probabilities. 

See Ahuja et al. JAMIA (2020) for details.

## Installation

If `devtools` is not installed, uncomment the code below and install it from CRAN.

``` r
# install.packages("devtools")
```

Run the code below to install `sureLDA` from GitHub:

``` r
devtools::install_github("celehs/sureLDA")
```

## Getting Started

Click [HERE](https://celehs.github.io/sureLDA/articles/example.html) to view a demo with a simulated example.

## References

Y. Ahuja, D. Zhou, Z. He, J. Sun, V. M. Castro, V. Gainer, S. N. Murphy, C. Hong, T. Cai (2020). sureLDA: A Multi-Disease Automated Phenotyping Method for the Electronic Health Record. JAMIA.
