# RSRF
## Description
Random Splitting Random Forest algorithms for mixed data including functional and scalar data. 
Response Types
- Continuous
- Categorical
- Survival

Version 1.0

## Origin 
It is an extension of the **Möller et al. (2016)** algorithm with the following new options:
- **Multiple** covariates
- Multiple **functioanl** covariates and **non-functioanl** covariates
- The respose type are **Survival** with censoring, **Categorical** and **Continuous**.
- The **Variable Importance Plot (VIP)** is enhanced.
- The random forest engine is from **randomForestSRC** by **Ishwaran et al (2022)**.

## Install in R

You can install the **RSRS** package from github and use it directly in *R* or *R-Studio*:
```
library(devtools)
install_github("mohammad-fayaz/RSRF")
```
or 
```
library(pak)
pak::pkg_install("mohammad-fayaz/RSRF")
```
## Vignettes
It has one vignette in this version:
- Random Splitting Random Forest for Categorical Response
You can run it step by step. 
> It is in the **Vignettes** folder. (PDF file **13 pages**)
[Vignettes - Link](https://github.com/mohammad-fayaz/RSRF/blob/master/vignettes/RSRF_Categorical_Output.pdf).

## References 
- [Oral Presentation] Fayaz M. and Abadi A.,**Functional Random Forest for Mixed Data**, The 41st Annual Conference of the International Society for Clinical Biostatistics (ISCB 2020)
- [Poster Presentation]  Fayaz M., Shakeri N., Abadi A. and Khodakarim S., **Random splitting random forest for survival analysis with non-functional and functional covariates in the EEG-fNIRS trial**, The 14th International Conference of the ERCIM WG on Computational and Methodological Statistics (CMStatistics 2021) 
- [Under Review] A Paper is submitted to the Peer-reviewed journal.
- Möller, Annette, Gerhard Tutz, and Jan Gertheiss. **"Random forests for functional covariates."** Journal of Chemometrics 30, no. 12 (2016): 715-725.
- Ishwaran, Hemant, Udaya B. Kogalur, and Maintainer Udaya B. Kogalur. "Package **‘randomForestSRC’**." breast 6 (2022): 1.

## Contacts
Mohammad Fayaz, PhD in Biostatistics 
ORCID: 0000-0002-5643-9763
