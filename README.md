[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="42%" src="inst/doc/scgft_logo.png"> 
</p>

# scGFT 

scGFT (single-cell Generative Fourier Transformer) is a generative model built
upon the principles of the Fourier Transform. It employs a one-shot
transformation paradigm to synthesize single-cell gene expression profiles that
reflect the natural biological variability found in authentic datasets.

---


## Installation

<details>
<br>

**scGFT** can be installed directly from this github with:

```{r}
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("Sanofi-GitHub/PMCB-scGFT", 
                         build_vignettes=FALSE)
```

</details>

---


## Usage Requirements

<details>
<br>

scGFT framework is designed to be compatible with the Seurat R analysis pipelines. 
To install, please run:

```{r}
# Enter commands in R (or R studio, if installed)
install.packages("Seurat")
install.packages("SeuratObject")
```

Visit [Seurat](https://satijalab.org/seurat/articles/install_v5) for more details.

</details>

---


## Functions Introduction

<details>
<br>

The scGFT package comprises only two functions: one to synthesize cells and a
second to evaluate the synthesis quality.

```r
# to synthsize cells
RunScGFT(object, nsynth, ncpmnts = 1, groups, scale.factor, cells = NULL)
```

`RunScGFT` requires, at a minimum, a Seurat object (`object`), the number of
desired cells to be synthesized (`nsynth`), a metadata variable indicating
groups of cells (`groups`), and the scale factor used for log-normalization of
the original data (`scale.factor`).

```r
# to evaluate synthsized cells
statsScGFT(object, groups)
```

`statsScGFT` requires a Seurat object that includes synthesized cells (`object`)
and the same character variable from the original object metadata used for
synthesis (`groups`). It calculates the likelihood that synthesized cells will
have the same identity (belong to same `groups`) as their original counterparts.

</details>

---


## Use Case 

<details>
<br>

#### Get demo files

We provided the dataset PRJEB44878 (Wohnhaas 2021), which comprises 34,200
processed cells derived from primary small airway epithelial cells from 
healthy individuals and patients with chronic obstructive pulmonary
disease. To download this dataset please run:

```{r}
# Enter commands in R (or R studio, if installed)
data_url <- "https://zenodo.org/records/11166226/files/COPD-PRJEB44878.rds"
# Define the path where you want to save the file (correct destination path
including the filename)
data_path <- "/path-to-destination/COPD-PRJEB44878.rds"
download.file(data_url, destfile = data_path, method = "auto")
```


#### Read data into R
```{r}
data_obj <- readRDS(data_path)
cnts <- data_obj$counts
mtd <- data_obj$metadata
```


#### Perform Seurat standard pipeline including synthesis process
```{r}
set.seed(1234)
sobj_synt <- CreateSeuratObject(counts=cnts,
                                meta.data=mtd) %>%
  NormalizeData(., normalization.method="LogNormalize", scale.factor=1e4) %>%
  FindVariableFeatures(., nfeatures=2000) %>%
  ScaleData(.) %>%
  RunPCA(., seed.use=42) %>%
  RunHarmony(., group.by.vars="sample") %>% # sample-specific batch correction
  FindNeighbors(., reduction="harmony", dims=1:30) %>%
  FindClusters(., random.seed=42) %>%
  # ================================
  # synthesis 1x cells (34,200), through modification of 10 complex components.
  RunScGFT(., nsynth=1*dim(.)[2], ncpmnts=10, groups="seurat_clusters", scale.factor=1e4) %>%
  # The combined dataset of original and synthetic cells undergoes another round. 
  # Re-normalization is not necessary as the new cells are synthesized from already normalized data.
  # ================================
  FindVariableFeatures(., nfeatures=2000) %>%
  ScaleData(.) %>%
  RunPCA(., seed.use=42) %>%
  RunHarmony(., group.by.vars=c("sample", "synthesized")) %>% # sample- and synthsis-specific batch correction
  FindNeighbors(., reduction="harmony", dims=1:30) %>%
  FindClusters(., random.seed=42) %>%
  RunUMAP(., reduction="harmony", seed.use=42, dims=1:30)
```

`RunScGFT` console outputs:

```{r}
Discrete fourier transform...
Inverse fourier transform...
synthesizing 34,200 cells...
5,621 cells synthesized...
10,838 cells synthesized...
11,716 cells synthesized...
12,485 cells synthesized...
13,208 cells synthesized...
13,852 cells synthesized...
14,471 cells synthesized...
14,846 cells synthesized...
15,186 cells synthesized...
15,360 cells synthesized...
15,458 cells synthesized...
20,054 cells synthesized...
24,080 cells synthesized...
26,972 cells synthesized...
29,009 cells synthesized...
30,503 cells synthesized...
31,886 cells synthesized...
33,260 cells synthesized...
34,200 cells synthesized...
Synthesis completed in: 2.79 min
Integrating data (1/4)
  [==================================================] 100% in  3m
Integrating data (2/4)
  [==================================================] 100% in 34s
Integrating data (3/4)
  [==================================================] 100% in  1m
Integrating data (4/4)
  [==================================================] 100% in  1m
A Seurat object with 68,400 cells, including 34,200 synthesized.
```


#### Evaluate synthsized cells

```{r}
statsScGFT(object=sobj_synt, groups="seurat_clusters")
```

```{r}
Synthesized cells: 34,200
Matching groups: 34,026
Accuracy (%): 99.49
```

Utilizing UMAP for a qualitative evaluation, we project synthesized and
real cells onto the embedded manifold:

<p align="center" width="100%">
<img style="width: 65%; height: auto;" src="inst/doc/panel_1_demo.png">
</p>

<p align="center" width="100%">
<img style="width: 85%; height: auto;" src="inst/doc/panel_2_demo.png">
</p>

We note that depending on the operating system used for calculations, the
results can be slightly different from the projected ones.

</details>

---


## Contact

<details>
<br>

For help and questions please contact the [scgft's maintenance team](mailto:nima.nouri@sanofi.com).

</details>

---
