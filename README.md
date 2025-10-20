# Evidence Contributions in Component Network Meta-Analysis from the Shortest-Path Approach

This repository contains the data, code, and results accompanying the manuscript  
**“Evidence Contributions in Component Network Meta-Analysis from the Shortest-Path Approach.”**  
The materials support the identification of pseudo paths and the quantification of evidence contributions in CNMA networks.

---

## 📂 Repository Contents

| File | Description |
|------|--------------|
| **HF_dt.csv** | Dataset of randomized controlled trials on heart failure used as an illustrative example in the manuscript. |
| **contrib_cnma.R** | Contains the R functions for identifying pseudo paths and computing their proportional contributions. |
| **Illustration.R** | Includes the R code for constructing CNMA networks, running CNMA models, and computing path-derived estimates across the three datasets analyzed in the study. |

---

##  1. Load Required Packages

```r
library(lpSolve)
library(netmeta)
library(MASS)
library(tidyverse)
library(devtools)
````

## 2. Source the Functions
```r
source("contrib_cnma.R")
````

The function **`cnma_contrib (x, "component", "inactive", random, kmax)`**  is used to compute and output the results of pseudo-path identification and flow allocation in CNMA.

Arguments:

  - x: An object of class 'discomb' or 'netcomb', available for both additive and interaction models.
  - component: A character string defining the target component whose contribution is to be analyzed.
  - inactive: A character string defining the inactive (reference) treatment component.
  - random: A logical indicating whether a random-effects CNMA should be conducted.
  - kmax: Maximum length of pseudo paths to consider. Default is 4.

Returns three data frames:
  - Paths: pseudo paths and their assigned flows (φ) representing contributions to the target component.
  - Proportion_for_edges: proportional contributions of each edge in the CNMA network.
  - Comparisons_in_CNMA: original study-level comparisons, their CNMA model representation, 
    and whether they contributed to the integrated direct estimate.

## 3. Load dataset
```r
HF_dt <- read.csv("HF_dt.csv")
````

## 4. Run illustration
```r
source("Illustration.R")
````

