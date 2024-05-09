[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="35%" src="vignettes/scgft_logo.png"> 
</p>

# scGFT 

scGFT (single-cell Generative Fourier Transformer) is a generative model built
upon the principles of the Fourier Transform. It employs a one-shot
transformation paradigm to synthesize single-cell gene expression profiles that
reflect the natural biological variability found in authentic datasets.



## Installation

**scGFT** can be installed directly from this github with:

```{r}
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Sanofi-GitHub/PMCB-scGFT", 
                         build_vignettes=FALSE)
```

## Usage

scGFT framework is designed to be compatible with the Seurat R analysis pipelines. 
To install, run:

```{r}
# Enter commands in R (or R studio, if installed)
install.packages("Seurat")
library("Seurat") 
```

Visit [Seurat](https://satijalab.org/seurat/articles/install_v5) for more details.
