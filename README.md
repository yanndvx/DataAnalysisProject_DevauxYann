# DataAnalysisProject - OxyMove Abstract

[![Language](https://img.shields.io/badge/Language-R%20%2F%20Python-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-In%20Progress-orange.svg)]()

## 👤 General Information
* **Author:** Yann Devaux
* **Program:** M1 IDIL DigiMove
* **Academic Year:** 2025-2026
* **Supervisor:** Pr. Mottet

---

## 📝 Project Overview
This repository contains the processing and analysis pipeline for a data subset from the **OxyMove** project. The primary objective is to measure hemodynamic changes and compare physiological responses between **eccentric** and **concentric** contraction modes.

### Experimental Design
The study follows a repeated measures design including:
* **Independent Variables:**
    * `Contraction_Mode`: Concentric vs. Eccentric.
    * `Intensity`: 30% vs. 50% of MVC (Maximum Voluntary Contraction).
* **Dependent Variable:** `EnvelopeMean_HbR_uM` (Deoxygenated Hemoglobin, currently, might change).
* **Sample Size:** 23 subjects (repeated measures).

![Experimental protocol](C:\Users\dvxya\Downloads\Chronologie_MOXY(1).png)
---

## 🛠️ Analysis Pipeline
The repository is structured into two distinct analysis tools:

### 1. NIRS Analysis (Near-Infrared Spectroscopy)
The workflow for NIRS data consists of three stages:
1.  **Loading & Preprocessing:** Raw `.snirf` files are processed in a Jupyter Notebook environment using:
    * [`mNIRS_V0.0.1.ipynb`](main_scritps/NIRS/NIRS_Treatment/mNIRS_V0.0.1.ipynb): Main working environment.
    * [`NIRS_functions_Oxymove.py`](main_scritps/NIRS/NIRS_Treatment/NIRS_functions_Oxymove.py): Function library for signal processing.
2.  **Export:** Generation of `.xlsx` long-form files (primarily for the *Semaxone* device).
3.  **Statistics:** Final statistical analysis performed via the RMarkdown script [`NIRS_Statistical_Pipeline.Rmd`](main_scritps/NIRS/NIRS_Stats/NIRS_Statistical_Pipeline.Rmd).

### 2. RPE Analysis (Rate of Perceived Exertion)
RPE data is analyzed directly in RStudio, provided the data is pre-organized in long-form `.csv` or `.xlss` files.

---

## 🚀 Getting Started

### Prerequisites
* **Python 3.12** (required libraries: `pandas`, `mne`, `matplotlib`, `mne-bids`, `matplotlib`, `glob`, `os`, `warnings`, `logging`, `pyqt6`) -> follow instructions in script preambule.
* **R 4.0+** (required libraries: `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `ggplot2`)

### Instructions
1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/yanndvx/DataAnalysisProject_DevauxYann.git](https://github.com/yanndvx/DataAnalysisProject_DevauxYann.git)
    ```
2.  **NIRS Processing:** Run the Python notebook to convert `.snirf` files into tabular data.
3.  **Statistical Analysis:** Open the two R project and run `.Rmd` file.

---

## ⚖️ License
This project is licensed under the **MIT License**. You are free to use, copy, and modify the code for your own projects, provided that the original author is credited.

```text
Copyright (c) 2026 Yann Devaux

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions...