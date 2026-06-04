# RGA Project - Source (`src`) Usage Guide

## Overview
This repository contains a refactored version of the RGA reaction generator and analyzer in `src/`.
The main executable script is `src/isodesmic.py`, which reads a YAML config file, loads species data, generates isodesmic reactions, computes similarity, and performs HoF analysis.

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
- `src/` - main source code for the RGA application
- `src/isodesmic.py` - CLI entry point
- `src/config.yaml` - example configuration file for `src/isodesmic.py`
- `src/CBSQB3_2019.10.25.txt` - example species database file used by the app
- `src/lam_CBSQB3.out` - example output file from a previous run
- `data/species/` - original data files for species and thermochemistry

## Setup for Running `src`
The `src/isodesmic.py` script looks for `config.yaml` in the current working directory when it runs.
Therefore, run it from the repository root and make sure a `config.yaml` file exists there.

If you do not already have `config.yaml` in the root, copy it from `src/`:

```powershell
copy .\src\config.yaml .\config.yaml
copy .\data\species\CBSQB3_2019.10.25.txt .\CBSQB3_2019.10.25.txt
```

## How to Run
From the repository root, execute:

```powershell
python src\isodesmic.py
```

The app will prompt you for:
1. executor name
2. SMILES of the target species

Example interactive input:

```text
Please enter the name of the executor: Lam
Input the SMILES of the species: CC
```

You can also automate the input with a temporary file:

```powershell
Set-Content -Path input.txt -Value "Lam`nCC`n`
"
cmd /c "python src\isodesmic.py < input.txt"
```

## Configuration
Edit `config.yaml` to control:
- species database file (`species_file`)
- result file name (`res_file`)
- executor name
- target species smiles
- stochastic / exhaustive generation mode
- number of cores
- reaction constraints
- similarity and analysis toggles

## Expected Output
The result file named in `res_file` (for example `lam_CBSQB3.out`) will contain:
- header block with executor and species
- constraint dump
- summary statistics
- list of reactions with similarity, HoF, and uncertainty
- execution details

## Troubleshooting
- If the script reports `config.yaml or config.json does not exist`, confirm you are running from the repo root and `config.yaml` exists there.
- If the script cannot find the species file, copy the correct data file to the same folder or adjust `species_file` in `config.yaml`.
- On Windows, `rdkit` install may require a conda channel or wheel; use your Python environment's recommended install method.

## Notes
- `src/isodesmic.py` is the main entry point for the refactored source tree.
- The root repository now has a runnable `config.yaml` and species file to support direct execution.
- This README is intended to help new contributors understand how to run and validate the `src` project.
