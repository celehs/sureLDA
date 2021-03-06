
# sureLDA: A Novel Multi-Disease Automated Phenotyping Method for the EHR

[![CRAN](https://www.r-pkg.org/badges/version/sureLDA)](https://CRAN.R-project.org/package=sureLDA)

## Overview

Surrogate-guided ensemble Latent Dirichlet Allocation (sureLDA) is a
label-free multidimensional phenotyping method. It first uses the
PheNorm or MAP algorithm to initialize probabilities based on two
surrogate features for each target disease, and then leverages these
probabilities to guide the LDA topic model to generate
phenotype-specific topics. Finally, it combines phenotype-feature counts
with surrogates via clustering ensemble to yield final phenotype
probabilities.

![](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/jamia/27/8/10.1093_jamia_ocaa079/5/m_ocaa079f1.png?Expires=1616108034&Signature=MSKCug91r5EvS7JsKwW-fHGbbohlAUsbEPOGjhe0yGx7eN6kgSd2lTlR30wvWgnKia267mHEcSoJ~9AsmDS2XjIwu4DTkiXPnU29Glv25Exqo4WY697oJb2wFYYYsRi12ASSjnAhENyA9V4VfdVl104XUj~A8rdShqzv0ZqwkG8~k2qWeJ-SwlIEled5SqvDnuxIR0qG5g9B4gecUeEa87cHHmrSJh3~umMn8HTzLCp1x8QVNuGLLpcSIyUXqkvPNRz6u7bByS4SEYN-W~LYOkFe6YJPwcTaWNrlgkrsBkXVDahLYpO5FiG7ZZEoGstEw8by9AMv7sCh18m7Vpjp3A__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

## Installation

Install stable version from CRAN:

``` r
install.packages("sureLDA")
```

Install development version from GitHub:

``` r
# install.packages("remotes")
devtools::install_github("remotes/sureLDA")
```

## Citation

Ahuja Y, Zhou D, He Z, Sun J, Castro VM, Gainer V, Murphy SN, Hong C,
Cai T. sureLDA: A multidisease automated phenotyping method for the
electronic health record. J Am Med Inform Assoc. 2020 Aug
1;27(8):1235-1243. doi: 10.1093/jamia/ocaa079. PMID: 32548637; PMCID:
PMC7481024.
