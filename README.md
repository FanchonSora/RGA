# RGA Project Usage Guide

## Overview
This repository contains the RGA reaction generator and analyzer. The primary execution directory is `RAG_initial/RGA`.
The main executable script is `isodesmic_multiple.py`, which reads a YAML config file (`config.yaml`), loads species data, generates isodesmic reactions, computes similarity, and performs HoF analysis.

## Requirements
Install Python dependencies before running the project.

Recommended packages:
- numpy
- pandas
- matplotlib
- rdkit
- scikit-learn
- sympy
- PuLP
- PyYAML

If you are using a Python virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

> Note: `rdkit` often requires a special install method depending on your Python distribution.

## Project Layout
- `RAG_initial/RGA/` - Main directory containing the application code and data.
- `RAG_initial/RGA/isodesmic_multiple.py` - Main CLI entry point for processing multiple reactions.
- `RAG_initial/RGA/config.yaml` - Configuration file defining input data and execution parameters.
- `RAG_initial/RGA/CBSQB3_2019.10.25.txt` - Species database file.
- `src/` - A cleaner, refactored version of the source code (currently retained for reference/future use).

## How to Run

1. Open your terminal and change the directory to the working folder:
```powershell
cd RAG_initial\RGA
```

2. Open `config.yaml` in your editor and configure the necessary inputs. The script `isodesmic_multiple.py` reads these values directly so you don't have to enter them interactively. Example configuration:
```yaml
files:
  species_file: CBSQB3_2019.10.25.txt
  res_file: ...\RGA_repo\output\'...'.out
  executor: Lam
  species_smiles: "C1=CC=C2C=CC=CC2=C1"
  calc_value: 38.7
```

3. Execute the script:
```powershell
python isodesmic_multiple.py
```

## Expected Output
The script will output its progress in the terminal and write the results to the path specified by `res_file` in your `config.yaml`. The output filename will be prefixed with the target SMILES string (e.g., `C1=CC=C2C=CC=CC2=C1_lam_CBSQB3.out`).

The result file will contain:
- Header block with executor and species details.
- Constraint parameters used for generation.
- List of reactions with similarity, HoF (Heat of Formation), and uncertainty.
- Execution details and performance summary.

## Troubleshooting
- If the script cannot find `config.yaml` or `CBSQB3_2019.10.25.txt`, make sure you have `cd` into the `RAG_initial\RGA` directory before executing the script.
- On Windows, installing `rdkit` via `pip` can sometimes fail; if so, consider using `conda install -c conda-forge rdkit` or a pre-compiled wheel.
