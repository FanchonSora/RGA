from utility_module import *
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
    # rxnFile = config.rxnFile
    # simFile = config.simFile
    resFile = config.resFile

    species = {}
    expts = {}
    uncerts = {}
    calcs = {}

    # Isodesmic reaction generation parameters
    sameCoef = config.sameCoef

    # Similarity checking parameters
    radius = config.radius
    chirality = config.chirality

    featureMethod = config.featureMethod
    if featureMethod.lower() == "ecfp":
        useFeatures = False
    else:
        useFeatures = True

    # Analysis parameters
    topReactions = config.top
    inputCalcData = None

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

            # toks = line.split()
            toks = line.split(",")

            if len(toks) == 3:
                uncert = 1
                calc = float(toks[2])
            else:
                uncert = float(toks[2])
                calc = float(toks[3])

            s = toks[0]
            expt = float(toks[1])

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

    print("Number of species in the database:", len(keys))

    # ------------------------------------
    start = time.clock()

    sz = len(keys)

    def generateIsodesmic(s1, keys, config):

        inputCalcData = calcs[s1]

        m1 = Chem.MolFromSmiles(s1)

        with open(resFile, "w") as fout:

            fout = printHeader(fout, s1, "Triet", config)

            if config.stochastic:
                # rxnCnt = generate_rxns_stochastic(keys, s1, m1, rxnFile, config)
                rxnCnt, rxnList = generate_rxns_stochastic(keys, s1, m1, fout, config)
            else:
                # rxnCnt = generate_rxns_complete(keys, s1, m1, rxnFile, config)
                rxnCnt, rxnList = generate_rxns_complete(keys, s1, m1, fout, config)

            # -------------------------------------------------

            if config.similarityOn:
                if rxnCnt > 0:
                    # reactions = checkSimilarity(rxnCnt=rxnCnt, rxnFile=rxnFile, simFile=simFile, config=config)
                    reactions = checkSimilarity(rxnList=rxnList, config=config, fout=fout)

                # -------------------------------------------------

                    if config.analysisOn:
                        # analyzeReactions(resFile=resFile, r1=s1, reactions=reactions, inputCalcData=inputCalcData, uncerts=uncerts, expts=expts, calcs=calcs, config=config)
                        res = analyzeReactions(fout=fout, r1=s1, reactions=reactions, inputCalcData=inputCalcData,
                                         uncerts=uncerts, expts=expts, calcs=calcs, config=config)

                    return res[3]

        return 0

    dataset_matrix = []

    print(np.mean([abs(calcs[key] - expts[key]) for key in species]))

    config.analysisOn = True
    config.similarityOn = True
    config.stochastic = False
    config.filter = False

    for index in range(sz):

        s1 = keys[index]

        # if s1 != "C=1=C=C=1":
        #     continue

        # if s1 != "CCC":
        #     continue

        # if s1 != "CCC" and s1 != "CC=C" and s1 != "CC#C" and s1 != "C1CCCCC1" and s1 != "C1=CC=CC=C1":
        #     continue

        print("\n\n###########################################\n")
        print("Considering", s1)
        print(index)

        curRow = [s1, expts[s1], calcs[s1]]

        print("All rxns")

        # HRCF

        # config.featureMethod = "ECFP"
        # config.radius = 0
        # config.top = -1
        # curRow.append(generateIsodesmic(s1, keys, config))

        config.filter = True

        config.featureMethod = "HRCF"
        config.radius = 4
        config.top = 0.1
        curRow.append(generateIsodesmic(s1, keys, config))

        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False
        #
        # # ECFP
        # config.featureMethod = "ECFP"
        # config.radius = 4
        # config.top = 0.1
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # # config.radical = True
        # # config.bondRadical = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False
        #
        # # FCFP
        # config.featureMethod = "FCFP"
        # config.radius = 4
        # config.top = 0.1
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # # config.radical = True
        # # config.bondRadical = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False
        #
        # # BoB
        # config.featureMethod = "BoB"
        # config.radius = 4
        # config.top = 0.1
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # # config.radical = True
        # # config.bondRadical = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False
        #
        # # CM
        # # config.featureMethod = "CM"
        # # config.radius = 4
        # # config.top = 0.1
        # # curRow.append(generateIsodesmic(s1, keys, config))
        # #
        # # config.ePair = True
        # # config.hydroBond = True
        # # config.normalBond = True
        # # # config.radical = True
        # # # config.bondRadical = True
        # # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # # config.hydroBond = False
        # # config.normalBond = False
        #
        # # AP
        # config.featureMethod = "AP"
        # config.radius = 4
        # config.top = 0.1
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # # config.radical = True
        # # config.bondRadical = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False
        #
        # # TT
        # config.featureMethod = "TT"
        # config.radius = 4
        # config.top = 0.1
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.ePair = True
        # config.hydroBond = True
        # config.normalBond = True
        # # config.radical = True
        # # config.bondRadical = True
        # curRow.append(generateIsodesmic(s1, keys, config))
        #
        # config.hydroBond = False
        # config.normalBond = False


        # Isogyric

        print("Isogyric")

        config.featureMethod = "ECFP"
        config.radius = 0
        config.top = -1
        config.ePair = True
        curRow.append(generateIsodesmic(s1, keys, config))

        # Isodesmic
        print("Isodesmic")

        # config.ePair = False
        config.hydroBond = True
        config.normalBond = True
        # config.radical = True
        # config.bondRadical = True
        curRow.append(generateIsodesmic(s1, keys, config))

        # Hypohomodesmotic
        print("Hypohomodesmotic")
        # config.hydroBond = False
        # config.normalBond = False
        config.atomHybridization = True
        config.atomWithHydro = True
        curRow.append(generateIsodesmic(s1, keys, config))

        # Homodesmotic
        print("Homodesmotic")
        config.atomWithHydro = False
        # config.atomHybridization = False
        config.atomHybridWithHydro = True
        config.bondHybrid = True
        curRow.append(generateIsodesmic(s1, keys, config))

        # Hyperhomodesmotic
        print("Hyperhomodesmotic")
        # config.hybridization = False
        config.bondHydro = True
        curRow.append(generateIsodesmic(s1, keys, config))
        #
        config.ePair = False
        config.atomHybridization = False
        config.atomHybridWithHydro = False
        config.hydroBond = False
        config.normalBond = False
        # config.radical = False
        # config.bondRadical = False
        config.bondHybrid = False
        config.bondHydro = False

        dataset_matrix.append(curRow)

    print(dataset_matrix)

    import pandas as pd
    results = pd.DataFrame(np.asarray(dataset_matrix))
    # results.to_csv('dataset_same_CBS_std.csv', index=False)
    # results.to_csv('dataset_same_AM1_std.csv', index=False)
    results.to_csv('testing.csv', index=False)
    # results.to_csv('dataset_same_AM1_no_filter.csv', index=False)
    # results.to_csv('dataset_same_AM1_v2.csv', index=False)
    # results.to_csv('similarity_same_CBS.csv', index=False)
    # results.to_csv('similarity_same_CBS.csv', index=False)
    # results.to_csv('similarity_same_AM1.csv', index=False)

    # --------------------------------------------------
    print("Execution time:", (time.clock() - start), "s.")

    raw_input("Press Enter to continue...")

if __name__ == '__main__':
    freeze_support()
    main()
