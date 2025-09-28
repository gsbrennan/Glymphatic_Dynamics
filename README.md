# GlymphaticDynamics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17220461.svg)](https://doi.org/10.5281/zenodo.17220461)

Pipeline for ROI-specific glymphatic clearance mapping and modeling proteinopathy cascades in the human brain


Code and data accompanying the paper:

> **Heterogeneity in human brain clearance adds resilience against tauopathy – a computational model informed by glymphatic MRI**  
> Georgia S. Brennan, Travis B. Thompson, Hadrien Oliveri, Vegard Vinje, Geir Ringstad, Per Kristian Eide, Alain Goriely, Marie E. Rognes  
> bioRxiv 2025.06.03.657596  
> https://doi.org/10.1101/2025.06.03.657596

---

## 📖 Overview

This repository provides code and example data for simulating how heterogeneity in human brain clearance (measured via glymphatic MRI) affects tauopathy progression in a network-based computational model.  


Key features:
- Extraction and analysis of subject-specific glymphatic clearance maps  
- Integration with a connectome-based neurodegeneration model; simulation of tau propagation and clearance under heterogeneous vs. homogeneous conditions  
- Visualization of progression patterns across cortical lobes and Braak regions  
- Analysis of both clearance maps and model outputs to evaluate heterogeneous clearance effects  

---


## License

- **Code**: Released under the [MIT License](LICENSE).  
- **Data and generated figures**: Released under the [CC-BY 4.0 License](LICENSE-DATA).  


## Quick Start (macOS/Linux)

We recommend running inside a Python virtual environment to ensure reproducibility.

```bash
git clone https://github.com/gsbrennan/GlymphaticDynamics.git
cd GlymphaticDynamics
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
jupyter lab

```
For Julia scripts and notebooks, an environment is provided via `Project.toml` (located in the scripts folder). To set it up, open Julia in the repo root and run:

```bash 
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```



## 📂 Repository Structure

```text
GlymphaticDynamics/
│
├── src/                        # Core source code
│   ├── clearance_pipeline/      # Full pipeline for clearance map processing
│   │   ├── 1-extract-field.py           # Step 1: extract MRI fields
│   │   ├── 2-drop-and-replace.py        # Step 2: clean data
│   │   ├── 3-cul-data.py                # Step 3: preprocessing
│   │   ├── 4-compute-clearance.py       # Step 4: fit clearance models
│   │   ├── 4b-average-computed-clearance.py
│   │   ├── 4c-identify-outliers.py
│   │   ├── 4d-replace-outliers.py
│   │   ├── 5a-normalize-patient-results.py
│   │   ├── 5b-amalgamate.py
│   │   ├── 6-average-cohort.py
│   │   ├── runall.sh                    # Shell script to run the full pipeline
│   │   ├── master-std33.graphml         # Standard connectome used in model
│   │   ├── raw-data/                    # Input raw clearance data
│   │   └── reformatted-data/            # Intermediate + processed outputs
│   │
│   ├── model/                           # Core computational model
│   │   └── bg_clearance_dynamics.py     # Network model of tau/clearance dynamics
│
├── scripts/                    
│   ├── patient_outputs/                 # Example patient-level simulation results
│   ├── plotting_cohorts.jl   # Julia script: cohort averages
│   ├── Project.toml                     # Julia environment for the scripts
│   └── Manifest.toml                    # Julia dependency lock file (auto-generated)
│
├── data/                                # Packaged data for reproducibility
│
├── notebooks/                           
│   └── data_paper_analysis.ipynb        # Main analysis notebook for paper
│
├── requirements.txt                     # Pinned dependencies (reproducible)
├── requirements_unpinned.txt            # Same packages, flexible versions
└── README.md                            # Project documentation

