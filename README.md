# Evidence Contributions in Component Network Meta-Analysis from the Shortest-Path Approach

This repository contains the data, codes, and results accompanying the manuscript
“Evidence Contributions in Component Network Meta-Analysis from the Shortest-Path Approach.”
The materials support the identification of pseudo paths and the quantification of evidence contributions in CNMA networks.

| File               | Description                                                                                                                                                    |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **HF_dt.csv**      | Dataset of randomized controlled trials on heart failure used as an illustrative example in the manuscript.                                                    |
| **contrib_cnma.R** | Contains the R functions for identifying pseudo paths and computing their proportional contributions.       |
| **Illustration.R** | Includes the R code for constructing CNMA networks, running CNMA models, and computing path-derived estimates across the three datasets analyzed in the study. |


# Load required packages
library(lpSolve)
library(netmeta)
library(MASS)
library(netmeta)
library(tidyverse)
library(devtools)

# Source the functions
source("contrib_cnma.R")

# Run illustration
source("Illustration.R")

# Load dataset
HF_dt <- read.csv("HF_dt.csv")

