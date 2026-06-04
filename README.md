# RGA (Reaction Generator and Analyzer)

RGA is a Python-based tool for automatically generating isodesmic reactions and calculating the standard Heat of Formation (HoF) at 298 K of molecules using similarity checks.

## Project Structure (Modernized)

- `run_demo.py`: Non-interactive demo runner. Use this to quickly test the engine.
- `config.yaml`: Core configuration file for setting up the constraints, stochastic/complete mode, and species info.
- `rga/`: The core Python package containing the generation, balancing, and similarity modules.
- `data/`: Contains database files (e.g., CBSQB3_2019.10.25.txt).
- `output/`: Generated analysis outputs.

## Requirements

RGA requires Python 3.8+ and the following packages:
- `rdkit`
- `numpy`
- `sympy`
- `pulp`
- `scikit-learn`
- `pyyaml`
- `pandas`
- `matplotlib`

You can install all dependencies via:
```bash
pip install -r requirements.txt
```

## Usage

To run a quick demo using the pre-configured settings in `config.yaml`:

```bash
python run_demo.py
```

Results will be generated in the `output/` directory.

## License
Copyright (c) Triet Le & Lam Huynh