# Amyloid-driven neuronal network collapse
Data coding  for simulations of amyloid-driven neuronal network collapse
# Amyloid-driven neuronal network collapse

This repository contains the code and data supporting the manuscript:

"Amyloid-Induced Network Resilience and Collapse in Alzheimer's Disease: Insights from Computational Modeling"

## Overview
We study the impact of progressive amyloid accumulation on neuronal network
structure and dynamics. Network topology is generated in MATLAB, while spiking
network simulations and analyses are performed in Python. However figures were plotted and processed in MATLAB.

## Repository structure
- `matlab/`: network construction and connectivity export
- `python/`: spiking simulations and analysis
- `data/`: connectivity matrices and simulation outputs
- `figures/`: graphical abstract and related figures

## Requirements
- MATLAB R2020a or later
- Python 3.9+
- numpy, scipy, matplotlib, networkx

## How to reproduce
1. Run MATLAB scripts in `matlab/` to generate network connectivity
2. Run Python scripts in `python/` using the exported matrices

## Contact
Corresponding author: Ediline Nguessap  
University of São Paulo  
Email: fonela@usp.br
