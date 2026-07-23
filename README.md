# RGA — Reaction Generator & Analyzer

## Overview

RGA generates isodesmic reactions for a target molecule, computes molecular similarity, and performs Heat-of-Formation (HoF) analysis using thermochemical species data.

The main entry point is `isodesmic_multiple.py` at the repo root. It reads `config.yaml`, loads species data, generates reactions, computes similarity, and outputs results.

---

## Project Layout

```
RGA_repo/
├── isodesmic_multiple.py    # ← Main CLI entry point
├── config.yaml              # ← Configuration (edit before running)
├── requirements.txt
│
├── src/
│   ├── core/                # Reaction generation & similarity engine
│   │   ├── rxngenerator_complete_parallel.py
│   │   ├── rxngenerator_stochastic_parallel.py
│   │   ├── similarity_parallel.py
│   │   ├── balancing_module.py
│   │   ├── bond_finder.py
│   │   ├── rxngenconfig.py
│   │   └── utility_module.py
│   ├── analysis/            # Analysis, ML & dataset modules
│   │   ├── analysis_module_new.py
│   │   ├── dataset_generator.py
│   │   ├── gaussian_process.py
│   │   ├── coulomb_matrix.py
│   │   ├── rga_ml.py
│   │   └── rga_selective_ml.py
│   └── tools/               # Utility & validation scripts
│       ├── ATcT_checker_new.py
│       ├── smiles_checking.py
│       ├── rxn_writer.py
│       ├── rxn_overlap_finder.py
│       ├── check_rxns.py
│       └── test_smiles.py
│
├── data/
│   ├── species/             # Thermochemical species databases
│   │   ├── CBSQB3_2019.10.25.txt   # ← default species file
│   │   ├── CBS-QB3_new.csv
│   │   ├── CBS-QB3_full.csv
│   │   ├── AM1_full.csv
│   │   ├── AM1.data
│   │   └── HOF_CBSQB3_2019.10.25.csv
│   ├── atct/
│   │   └── ATCT_DATABASE-2018.06.14.xlsx
│   ├── literature/
│   │   └── literature_data.csv
│   └── ml/                  # ML-ready datasets
│       ├── datasets/
│       └── similarity/
│
└── output/                  # Generated result files (gitignored)
```

---

## Requirements

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Key packages: `numpy`, `pandas`, `matplotlib`, `rdkit`, `scikit-learn`, `sympy`, `PuLP`, `PyYAML`

> **Note:** On Windows, `rdkit` via pip can fail — use `conda install -c conda-forge rdkit` if needed.

---

## How to Run

1. Edit `config.yaml` at the repo root:

```yaml
files:
  species_file: data/species/CBSQB3_2019.10.25.txt
  res_file: output/my_result.out
  executor: YourName
  species_smiles: "C1=CC=C2C=CC=CC2=C1"
  calc_value: 38.7
```

2. Run from the repo root:

```powershell
python isodesmic_multiple.py
```

---

## Expected Output

Results are written to `output/` (path set in `config.yaml`). Each output file contains:
- Header block (executor, species, parameters)
- List of reactions with similarity score, HoF, and uncertainty
- Execution summary and performance stats

---

## Troubleshooting

- **`config.yaml` not found**: Run from the repo root, not a subdirectory.
- **Import errors**: Ensure `.venv` is activated and all packages are installed.
- **rdkit install fails on Windows**: Use `conda install -c conda-forge rdkit`.
