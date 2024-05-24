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
RunScGFT(object, nsynth, ncpmnts = 1, groups, cells = NULL)
```

`RunScGFT` requires, at a minimum, a Seurat object (`object`), the number of
desired cells to be synthesized (`nsynth`), and a metadata variable indicating
groups of cells (`groups`).

```r
# to evaluate synthsized cells
statsScGFT(object, groups)
```

`statsScGFT` requires a Seurat object that includes synthesized cells (`object`)
and the same character variable from the original object metadata used for
synthesis (`groups`). It calculates the likelihood that synthesized cells will
have the same identity (belong to same `groups`) as their original counterparts.
It also reports the relative deviation of synthesized gene expression profiles from
original cells.

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
data_path <- "~/COPD-PRJEB44878.rds" # correct destination path includes the filename
download.file(url=data_url, destfile=data_path, method="auto")
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
  NormalizeData(., normalization.method="LogNormalize", scale.factor=1e6) %>%
  FindVariableFeatures(., nfeatures=2000) %>%
  ScaleData(.) %>%
  RunPCA(., seed.use=42) %>%
  RunHarmony(., group.by.vars="sample") %>% # sample-specific batch correction
  FindNeighbors(., reduction="harmony", dims=1:30) %>%
  FindClusters(., random.seed=42) %>%
  # ================================
  # synthesis 1x cells (34,200), through modification of 10 complex components.
  RunScGFT(., nsynth=1*dim(.)[2], ncpmnts=10, groups="seurat_clusters") %>%
  # The combined dataset of original and synthetic cells undergoes another round. 
  # Re-normalization is not needed as the new cells are synthesized from already normalized data.
  # ================================
  FindVariableFeatures(., nfeatures=2000) %>%
  ScaleData(.) %>%
  RunPCA(., seed.use=42) %>%
  RunHarmony(., group.by.vars="sample") %>% # sample-specific batch correction
  FindNeighbors(., reduction="harmony", dims=1:30) %>%
  FindClusters(., random.seed=42) %>%
  RunUMAP(., reduction="harmony", seed.use=42, dims=1:30)
```

`RunScGFT` console outputs:

```{r}
Discrete fourier transform...
Inverse fourier transform...
synthesizing 34,200 cells...
4,902 cells synthesized...
9,712 cells synthesized...
14,389 cells synthesized...
17,116 cells synthesized...
19,708 cells synthesized...
21,988 cells synthesized...
24,052 cells synthesized...
25,812 cells synthesized...
27,264 cells synthesized...
28,563 cells synthesized...
29,474 cells synthesized...
30,334 cells synthesized...
31,123 cells synthesized...
31,898 cells synthesized...
32,570 cells synthesized...
33,226 cells synthesized...
33,600 cells synthesized...
33,882 cells synthesized...
34,057 cells synthesized...
34,150 cells synthesized...
34,200 cells synthesized...
Synthesis completed in: 5.18 min
Integrating data (1/4)
  [==================================================] 100% in  3m
Integrating data (2/4)
  [==================================================] 100% in 37s
Integrating data (3/4)
  [==================================================] 100% in  1m
Integrating data (4/4)
  [==================================================] 100% in  2m
A Seurat object with 68,400 cells, including 34,200 synthesized.
```


#### Evaluate synthsized cells

```{r}
statsScGFT(object=sobj_synt, groups="seurat_clusters")
```

`statsScGFT` console outputs:

```{r}
Synthesized cells: 34,200
Matching groups: 33,585
Accuracy (%): 98.2
Calculating deviation from originals...
  [==================================================] 100% in  3m
  [==================================================] 100% in  1m
Deviation (%): 1.55 +/- 0.38
```

Utilizing UMAP for a qualitative evaluation, we project synthesized and
real cells onto the embedded manifold:

<p align="center" width="100%">
<img style="width: 65%; height: auto;" src="inst/doc/panel_1_demo.png">
</p>

<p align="center" width="100%">
<img style="width: 85%; height: auto;" src="inst/doc/panel_2_demo.png">
</p>

Depending on the operating system used for calculations, the
results can be different from the projected ones.

#### Perform cell-based synthesis:

In this showcase, we expand rare epithelial subtypes, including aberrant
basaloid cells, PNECs, and ionocytes, each comprising less than 0.3% of the
population:

```{r}
set.seed(1234)
sobj_synt <- CreateSeuratObject(counts=cnts,
                                meta.data=mtd) %>%
  NormalizeData(., normalization.method="LogNormalize", scale.factor=1e6) %>%
  FindVariableFeatures(., nfeatures=2000) %>%
  # ================================
  # synthesis 1,000, through modification of 10 complex components, for each of given annotated rare epithelial subtypes
  RunScGFT(., nsynth=1000, ncpmnts=10, cells = "S2_ACGGAGAGTTCCCGAG-1") %>% # pre-annotated "Ionocyte" cell
  RunScGFT(., nsynth=1000, ncpmnts=10, cells = "S1_ATTACTCTCGTTGCCT-1") %>% # pre-annotated"Pnec" cell
  RunScGFT(., nsynth=1000, ncpmnts=10, cells = "S1_AAGCCGCGTGCCTGCA-1") %>% # pre-annotated "Aberrant_basaloid" cell
  # ================================
  FindVariableFeatures(., nfeatures=2000) %>%
  ScaleData(.) %>%
  RunPCA(., seed.use=42) %>%
  RunHarmony(., group.by.vars=c("sample")) %>% # sample-specific batch correction
  FindNeighbors(., reduction="harmony", dims=1:30) %>%
  FindClusters(., random.seed=42) %>%
  RunUMAP(., reduction="harmony", seed.use=42, dims=1:30)
```


#### Evaluate synthsized cells

Next, we evaluate the consistency of cell types in synthesized cells relative to
the originals. Cells goes through another round of cell type annotation using
[Sargent](https://github.com/Sanofi-Public/PMCB-Sargent), an automated, 
cluster-free, score-based annotation method that classifies cell types based on
distinct gene expression [markers](inst/doc/epithelial_subtypes_markers.xlsx) 
(a helper script can be found [here](inst/doc/sargent_helper.R)). 
Then, the annotations of the synthesized cells are evaluated by:

```{r}
statsScGFT(object=sobj_synt, groups="sargent_celltype")
```

`statsScGFT` console outputs:

```{r}
Synthesized cells: 3,000
Matching groups: 2,998
Accuracy (%): 99.93
Calculating deviation from originals...
  [==================================================] 100% in  3s
Deviation (%): 0.62 +/- 0.02
```

Utilizing UMAP for a qualitative evaluation, we project synthesized and
real cells onto the embedded manifold:

<p align="center" width="100%">
<img style="width: 65%; height: auto;" src="inst/doc/panel_3_demo.png">
</p>

Depending on the operating system used for calculations, the
results can be different from the projected ones.


</details>

---


## Contact

<details>
<br>

For help and questions please contact the [scgft's maintenance team](mailto:nima.nouri@sanofi.com).

</details>

---
