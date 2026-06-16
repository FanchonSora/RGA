from utility_module import time, Chem, os, freeze_support, printHeader, printExecutionDetails, np
from datetime import datetime
from rxngenconfig import Config
from rxngenerator_complete_parallel import generate_rxns_complete
from rxngenerator_stochastic_parallel import generate_rxns_stochastic
from similarity_parallel import checkSimilarity
from analysis_module_new import analyzeReactions

def main():

    print("\n*** ------------------------------------ ***")
    print("***           ATcT Checker               ***")
    print("***     Copyright @ Triet Le 2018        ***")
    print("*** ------------------------------------ ***\n")

    # Input species parameters
    config = Config()
    config.parseConfig()

    speciesFile = config.speciesFile
    resFile = config.resFile

    species = {}
    expts = {}
    uncerts = {}
    calcs = {}

    if not os.path.exists(speciesFile):
        print("The species file", speciesFile, "does not exist. Please check again !!!")
        input("Press Enter to continue...")
        exit(0)

    # speciesFile = "literature_data.csv"

    with open(speciesFile) as f:
        for line in f:

            line = line.rstrip("\n")

            if line == "":
                continue

            toks = line.split(",")

            if len(toks) == 3:
                uncert = 1
                calc = float(toks[2])
            else:
                uncert = float(toks[2])
                calc = float(toks[3])

            s = toks[0]
            expt = float(toks[1])

            # toks = line.split()
            # toks = line.split(",")
            #
            # s = toks[0]
            # expt = float(toks[1])
            # uncert = float(toks[2])
            # calc = float(toks[3])

            # s = toks[-1].replace("\"", "")
            # expt = float(toks[1].replace("\"", ""))
            # calc = float(toks[2].replace("\"", ""))
            # uncert = float(toks[3].replace("\"", ""))

            # species[line] = 1
            species[s] = 1
            expts[s] = expt
            uncerts[s] = uncert
            calcs[s] = calc

    keys = sorted(species.keys())

    print("Mean of calculated method:", np.mean([abs(calcs[key] - expts[key]) for key in species]))
    print("Standard deviation of calculated method:", np.std([abs(calcs[key] - expts[key]) for key in species]))
    print("Mean uncertainty:", np.mean(list(uncerts.values())))

    # print(keys.index("[CH2]C=C"))

    print("Number of species in the database:", len(keys))

    # ------------------------------------

    excludedSpecies = ["C", "C#C", "C#CC#C", "C1=CC1CC2C=C2", "C1=CC=CC=C1", "C=1=C=C=1", "C=C", "C=CC=C", "C=CCCC=C", "CC", "CC(C)(C)C", "CC(C)CC(C)(C)C", "CC1=CC=CC=C1",
                       "[CH3]", "[CH]1C=C1", "[CH]C[CH]", "[C]1C=C1", "[C]=C[CH]", "[C]C=[C]", "[C]C[CH]"]

    excludedSpecies_2 = ["C", "C#C", "C1=CC=CC=C1", "C=C", "C=C=C", "C=CC=C", "C=CCCC=C", "CC", "[CH3]"]

    sz = len(keys)

    oldExpts = expts.copy()

    usedSpecies = ['CC', 'C=CC', 'C#CC', 'C1CCCCC1', 'C1=CC=CC=C1']

    usedSpecies = ["CC"]

    oldKeys = keys[:]

    start = time.perf_counter()

    for iter in range(1):

        datasets = []

    # for split_point in range(5, sz, 100):
    #
    #     datasets = []
    #
    #     print(split_point)
    #
    #     keys = oldKeys[:split_point]

        print("\n\nIteration", iter + 1)
        # print("Size of key:", split_point)
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

        # keys = np.random.permutation(keys)

        for index in range(len(keys)):
        # for index in range(10):

            s1 = keys[index]

            # if s1 not in usedSpecies:
            #     continue

            # if s1 in excludedSpecies:
            #     continue

            # if s1 == "[HH]" or s1 == "C1=CC=CC=C1":
            #     continue

            print("\n\n###########################################\n")
            print("Considering", s1)
            print(index)

            inputCalcData = calcs[s1]

            m1 = Chem.MolFromSmiles(s1)

            curResFile = s1.replace("/", "'").replace("\\", "`") + "_" + resFile

            with open(curResFile, "w") as fout:

                cur_start = time.perf_counter()

                fout = printHeader(fout, s1, "Triet", config, datetime.now().year)

                if config.stochastic:
                    rxnCnt, rxnList = generate_rxns_stochastic(keys, s1, m1, fout, config)
                else:
                    rxnCnt, rxnList = generate_rxns_complete(keys, s1, m1, fout, config)

                # -------------------------------------------------

                if config.similarityOn:
                    if rxnCnt > 0:
                        reactions = checkSimilarity(rxnList=rxnList, config=config, fout=fout)

                        # -------------------------------------------------

                        if config.analysisOn:
                            datasets.append(analyzeReactions(fout=fout, r1=s1, reactions=reactions, inputCalcData=inputCalcData,
                                             uncerts=uncerts, expts=expts, calcs=calcs, config=config))

                exe_time = time.perf_counter() - cur_start
                print("Current execution time of", s1, ":", exe_time, "s.")

                printExecutionDetails(fout, exe_time)

        # --------------------------------------------------

        import pandas as pd

        indices = list(range(1, len(datasets) + 1))
        columns = ['No.', 'Species', 'No. rxns', 'True HoF', 'Pred. HoF', 'Absolute Difference', 'Std. HoF',
                   'True Uncertainty', 'Pred. Uncertainty']
        datasets = pd.DataFrame(np.c_[indices, np.asarray(datasets)], columns=columns)

        print("Average HoF error:", np.mean(datasets['Absolute Difference'].values.astype(np.float64)))
        print("Std HoF error:", np.std(datasets['Absolute Difference'].values.astype(np.float64)))
        print("Entropy mean:", np.mean(datasets['Std. HoF'].values.astype(np.float64)))
        print("Entropy std:", np.std(datasets['Std. HoF'].values.astype(np.float64)))
        print("Total no. of rxns:", np.sum(datasets['No. rxns'].values.astype(np.int32)))

    # diffs = 0
    # for s in oldExpts:
    #     print(s, abs(expts[s] - oldExpts[s]))
    #     diffs += abs(expts[s] - oldExpts[s])
    #
    # print("Mean difference", diffs / len(keys))

    # datasets.to_csv('HoF_sameCoef_all_species-CBSQB3.csv', index=False)
    # datasets.to_csv('HoF_sameCoef_cos_top10_FCFP_species-CBSQB3.csv', index=False)
    # datasets.to_csv('HoF_sameCoef_all_species-CBSQB3_noiso.csv', index=False)
    # datasets.to_csv('HoF_sameCoef_cos_top10_highres_noiso.csv', index=False)
    # datasets.to_csv('HoF_sameCoef_cos_top10_highres.csv', index=False)

    # --------------------------------------------------
    print("Execution time:", (time.perf_counter() - start), "s.")

    input("Press Enter to continue...")

if __name__ == '__main__':
    freeze_support()
    main()
