# Two-Stage HMREM Model for Analysing Relational Event Data

This repository contains the reproducible code and example scripts accompanying the article:

**"Two-Stage HMREM Model for Analysing Relational Event Data"**

---

## Overview

In this work, we introduce a **two-stage modelling framework (HMREM)** combining:

1️⃣ **Hidden Markov Model (HMM)** — to infer latent states governing the dynamics of relational events,  
2️⃣ **Relational Event Model (REM)** — to model event generation conditioned on the inferred HMM states.

This approach allows us to capture **dynamic regime shifts** and **state-dependent interaction mechanisms** in relational event sequences.

---

## Repository Contents

This repository includes the following code files:

| File                          | Description                                           |
|-------------------------------|-----------------------------------------------------|
| `2stateHMREM.Rmd`              | HMREM simulation and analysis for **2-state model** |
| `3stateHMREM.Rmd`              | HMREM simulation and analysis for **3-state model** |
| `Apollo2StateCase.Rmd`         | Real-world **Apollo 13 case study** using HMREM     |
| `predictive_performance.Rmd`   | Functions for **predictive performance evaluation** |

---

## Reproducibility

✅ All scripts can be run independently.  
✅ The case study uses the provided Apollo 13 event data (`Apollo13_completeData.RData`).  
✅ The repository is fully self-contained.

---

## Dependencies

The code relies on the following R packages:

- `remify`  
- `remstats`  
- `remstimate`  
- `remulate`  
- `momentuHMM`  
- `dplyr`  
- `ggplot2`  
- `plotly`  
- `readxl`  
- `fastDummies`  
- `abind`

You can install these packages using:

```r
install.packages(c("remify", "remstats", "remstimate", "momentuHMM", "dplyr", "ggplot2", "plotly", "readxl", "fastDummies", "abind"))
