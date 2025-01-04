# Dynamic Supervised Principal Component Analysis for Classification -- R Package
`DSPCA` implements the algorithm of *Dynamic Supervised Principal Component Analysis for Classification* ([arXiv](https://arxiv.org/abs/2411.01820)).

## Installation
Install from GitHub using the `devtools` package:

``` R
devtools::install_github("ywwry66/DSPCA")
```

## Usage
Use `cv_dspca` to fit a model and use `predict` to predict the class labels for new observations. Documentations can be accessed via `?DSPCA::cv_dspca` and `?DSPCA::predict.cv_dspca` in an interactive `R` session.
