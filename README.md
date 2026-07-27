# Improving Parent Identification in Plant Breeding
*This repository accompanies the manuscript which is currenly under review*

**Improving parent identification in plant breeding using combined similarity and likelihood metrics**

Authors: *Alexandra Nirsha¹, John Riviere², Martin Spanoghe², Deborah Lanterbecq¹,²,³*

¹ Centre pour l’Agronomie et l’Agro-Industrie de la Province de Hainaut (asbl CARAH), 7800 Ath, Belgium

² Laboratoire de Biotechnologie et Biologie Appliquée, Haute Ecole Provinciale de Hainaut-CONDORCET, 7000 Mons, Belgium

³ Hainaut Analyses (HA), 7000 Mons, Belgium

Journal: *Discover Plants* (under review)


---

# Overview

This repository provides a complete workflow for evaluating parent–offspring relationships from binary molecular marker data.

The implemented pipeline combines:

- Jaccard similarity analysis
- LOD score computation
- comparative analysis between both approaches
- graphical visualization of the results

The workflow reproduces the analyses presented in the accompanying manuscript.

## Contents

- pedigree preparation functions
- Jaccard similarity analysis
- LOD score analysis
- combined Jaccard–LOD comparison
- plotting utilities
- scripts reproducing the manuscript figures
- anonymized example dataset

## Workflow

1. Load the project.
2. Load the example dataset.
3. Prepare pedigree subsets.
4. Compute similarity matrices.
5. Perform Jaccard analysis.
6. Perform LOD analysis.
7. Compare both methods.
8. Generate figures.

---

# Installation

## Clone the repository

```bash
git clone https://github.com/ale-nir/parent-identification-plant-breeding.git
```
or download the repository as a ZIP archive from GitHub and extract it.

---

## Repository structure

```bash
.
├── R/
├── src/
├── data/
├── docs/
├── scripts/
├── examples/
├── setup.R
├── README.md
└── LICENSE
```

---

## Open the project

Open the project in RStudio (recommended) or set the working directory to the repository root.

```r
setwd("path/to/repository")
```

---

## Quick start

Source the setup script. 

```r
source("setup.R")
```

This script automatically:

- installs missing packages (if necessary);
- loads all required functions;
- compiles the C++ source files;
- loads the example dataset.

---

# Documentation

Detailed documentation for each analysis step is available below.

- [Pedigree preparation](docs/pedigree_preparation.md)
- [Jaccard analysis](docs/jaccard_analysis.md)
- [LOD analysis](docs/LOD_analysis.md)
- [Combined analysis](docs/combined_analysis.md)
- [Plotting figures](docs/plotting.md)

---

## Reproducing the manuscript results

The complete workflow used to generate the manuscript figures is available in:

```r
source("examples/reproduce_manuscript_figures.R")
```

This script reproduces the analyses used for the publication figures and generates the intermediate analysis objects from which the Jaccard summary tables reported in the manuscript are also derived.

---

# Citation

If you use this repository in your work, please cite:

Nirsha A., Rivière J., Spanonghe M., Lanterbecq D. (2026) *Improving parent identification in plant breeding using combined similarity and likelihood metrics.* (under review)

The citation will be updated once the manuscript is published.

---

# License

This repository is released under the Academic Research License.

The source code is provided to ensure scientific transparency and reproducibility of the associated publication.

Commercial use or integration into proprietary software requires prior written authorization from CARAH asbl.

For commercial licensing opportunities or other permissions, please contact: deborah.lanterbecq@condorcet.be | a.nirsha@carah.be
