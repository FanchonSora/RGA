"""
RGA: Reaction Generator and Analyzer - Demo Runner
Non-interactive entry point for demonstration.

Usage:
    python run_demo.py
"""

import time
import os
from multiprocessing import freeze_support
from rdkit import Chem

from rga.config import Config
from rga.generator_complete import generate_rxns_complete
from rga.generator_stochastic import generate_rxns_stochastic
from rga.similarity import checkSimilarity
from rga.analysis import analyzeReactions
from rga.utils import printHeader, printExecutionDetails


def main():
    print("\n*** ------------------------------------- ***")
    print("*** RGA: Reaction Generator and Analyzer  ***")
    print("***   Copyright @ Triet Le & Lam Huynh    ***")
    print("***         Version 2.0 (2026)            ***")
    print("*** ------------------------------------- ***\n")

    # Parse config
    config = Config()
    if not config.parseConfig():
        print("ERROR: config file (config.yaml or config.json) not found or invalid!")
        return

    speciesFile = config.speciesFile
    resFile = config.resFile

    species = {}
    expts = {}
    uncerts = {}
    calcs = {}

    if not os.path.exists(speciesFile):
        print(f"ERROR: Species file '{speciesFile}' does not exist!")
        return

    # Read species database
    with open(speciesFile) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "":
                continue

            toks = line.split()

            if len(toks) == 3:
                uncert = 1
                calc = float(toks[2])
            else:
                uncert = float(toks[2])
                calc = float(toks[3])

            s = toks[0]
            expt = float(toks[1])

            species[s] = 1
            expts[s] = expt
            uncerts[s] = uncert
            calcs[s] = calc

    keys = list(species.keys())
    print(f"Number of species in the database: {len(keys)}\n")

    # Get input species from config
    executor = config.executor if config.executor else "Demo"
    input_species = str(config.speciesSmiles).split()
    calc_vals = str(config.calcValue).split()

    # Create output directory
    os.makedirs(os.path.dirname(resFile) if os.path.dirname(resFile) else ".", exist_ok=True)

    for index in range(len(input_species)):
        s1 = input_species[index]

        print(f"\n########################################\n")
        print(f"Processing: {s1}")

        m1 = Chem.MolFromSmiles(s1)
        if not m1:
            print(f"ERROR: Invalid SMILES '{s1}'. Skipping...")
            continue

        inputCalcData = float(calc_vals[index])
        curResFile = s1.replace("/", "'").replace("\\", "`") + "_" + os.path.basename(resFile)
        curResFile = os.path.join(os.path.dirname(resFile), curResFile)

        with open(curResFile, "w") as fout:
            fout = printHeader(fout, s1, executor, config)

            start = time.perf_counter()

            if config.stochastic:
                rxnCnt, rxnList = generate_rxns_stochastic(keys, s1, m1, fout, config)
            else:
                rxnCnt, rxnList = generate_rxns_complete(keys, s1, m1, fout, config)

            print(f"Reactions generated: {rxnCnt}")

            if config.similarityOn and rxnCnt > 0:
                reactions = checkSimilarity(rxnList=rxnList, config=config, fout=fout)

                if config.analysisOn:
                    analyzeReactions(
                        fout=fout, r1=s1, reactions=reactions,
                        inputCalcData=inputCalcData,
                        uncerts=uncerts, expts=expts, calcs=calcs,
                        config=config
                    )

            exe_time = time.perf_counter() - start
            print(f"Execution time: {exe_time:.2f} s.")

            printExecutionDetails(fout, exe_time)

        print(f"\nResults written to: {curResFile}")

    print("\n*** Demo completed successfully! ***")


if __name__ == '__main__':
    freeze_support()
    main()
