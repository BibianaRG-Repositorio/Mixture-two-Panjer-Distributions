

# Mixture-two-Panjer-Distributions

This repository contains the **R code used for the statistical analysis and results presented in the published research article**:

**“Mixture of Two Panjer Distributions: Optimizing the Insurance Claims Model with Factor Moment Estimation.”**

The scripts included in this repository reproduce the main steps of the study, including:

- Data preparation
- Estimation of distribution parameters using the **factorial moment method**
- Fitting several discrete distributions commonly used in **insurance claim modeling**
- Construction of a **mixture of two Panjer distributions**
- Goodness-of-fit evaluation using the **chi-square test**
- Visualization and comparison of fitted distributions

The analysis was performed using the **R statistical programming language**.

---

# Publication

This repository accompanies the following **peer-reviewed article**:

**Jiménez Moscoso, J. A., & Riaño-Gaona, N. B. (2025).**  
*Mixture of two Panjer Distributions: Optimizing the Insurance Claims Model with Factor Moment Estimation.*  

**Journal:** Ciencia en Desarrollo  
**Volume:** 16(2)  
**Pages:** 136–141  

**DOI:**  
https://doi.org/10.19053/uptc.01217488.v16.n2.2025.18930

**Publication date:** July 20, 2025  

**Journal issue:**  
Vol. 16, No. 2 (July–December 2025)

---

# Abstract

The main objective of this study is to model the number of events, such as **insurance claims reported during a given period**. Traditional models often use discrete distributions such as **Poisson, Geometric, Binomial, and Poisson–Gamma**, which belong to the **Panjer family of distributions**.

This research proposes a **mixture of two Panjer frequency distributions** to improve the modeling of claim counts. The parameters of the distributions are estimated using the **factorial moment method**.

The empirical results show, based on the **chi-square goodness-of-fit test**, that the proposed mixture model provides a **better fit than a single discrete distribution** for the analyzed dataset.

All computations and simulations were performed using **R**.

---

# Repository Structure
.
├── scripts
│ ├── distribution_models.R
│ └── mixture_panjer_model.R
├── data
│ └── claims_frequency_data.csv
├── figures
│ └── distribution_comparison_plot.png
└── README.md


