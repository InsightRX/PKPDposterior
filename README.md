# PKPDposterior

## Introduction

This software package provides convenience functions around Stan (CmdStanR), Torsten, and
PKPDsim. The goal of the package is to make it easier to full-Bayesian inference
in the context of model-informed precision dosing.

The software is intended for use in education and research, and is **not
intended for clinical use**. By design, this package does not include
functionality that directly provides dosing advice.

## Installation and getting started

Install the latest version of PKPDposterior using remotes or devtools:

``` r
remotes::install_github("InsightRX/PKPDposterior")
devtools::install_github("InsightRX/PKPDposterior")
```

PKPDposterior depends on Torsten and Stan. If you don't already have these
software working, follow the instructions for installing
[Torsten](https://github.com/metrumresearchgroup/Torsten). Then, add the file
path to your Stan installation location to your environment variables.

```r
Sys.setenv(STAN_PATH = "path/to/stan")
```

To avoid running the above command every time you restart your R session,
consider adding this environment variable to your .Rprofile or other
configuration file.

Prepare your R session to use Stan by running the following line of code:

```r
cmdstanr::set_cmdstan_path(
    path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)
```

Great! Now you should be all ready to start modeling!

## Contributing

We welcome input from the community:

- If you think you have encountered a bug, please [submit an issue](https://github.com/InsightRX/PKPDposterior/issues)
on the GitHub page. Please include a reproducible example of the unexpected
behavior.

- Please [open a pull request](https://github.com/InsightRX/PKPDposterior/pulls) if
you have a fix or updates that would improve the package. If you're not sure if
your proposed changes are useful or within scope of the package, feel free to
contact one of the authors of this package.

## Disclaimer

The functionality in this R package is provided "as is". While its authors
adhere to software development best practices, the software may still contain
unintended errors.

InsightRX Inc. and the authors of this package can not be held liable for any
damages resulting from any use of this software. By the use of this software
package, the user waives all warranties, expressed or implied, including any
warranties to the accuracy, quality or suitability of InsightRX for any
particular purpose, either medical or non-medical.


<div align="right">
© <img src="man/figures/insightrx_logo_color.png" alt="InsightRX logo" width="120" />
</div>
