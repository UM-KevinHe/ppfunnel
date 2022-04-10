
# ppfunnel

## Overview

`ppfunnel` creates elegant funnel plots for profiling health care
providers.

## Installation

``` r
# install the development version from GitHub:
devtools::install_github("UM-KevinHe/ppfunnel")
```

## Usage

``` r
library(tidyverse)
library(ppfunnel)
mytheme <- theme(legend.justification=c(1,1), legend.position=c(1,1),
                 legend.box="horizontal", legend.text=element_text(size=14),
                 axis.title=element_text(size=14),
                 axis.text=element_text(size=14),
                 plot.title=element_text(hjust=0), text=element_text(size=13))
# Poisson outcome with method = "FE"
SHR %>% ppfunnel(indiv.data=F, legend.justification=c(1,1),
                 legend.position=c(1,1), legend.box="horizontal",
                 legend.text=element_text(size=14), axis.title=element_text(size=14),
                 axis.text=element_text(size=14), plot.title=element_text(hjust=0),
                 text=element_text(size=13))
# alternative use
SHR %>% ppfunnel(indiv.data=F) + mytheme
# Poisson outcome with method = "indivEN.0meanapprox"
SHR %>% ppfunnel(indiv.data=F, method="indivEN.0meanapprox") + mytheme
```
