from utility_module import np, os, time, Chem
from rxngenconfig import Config
from rxngenerator_complete import generate_rxns_complete
from rxngenerator_stochastic import generate_rxns_stochastic
from similarity_module import checkSimilarity
from analysis_module import analyzeReactions

print("\n*** ------------------------------------ ***")
print("***           ATcT Checker               ***")
print("***     Copyright @ Triet Le 2018        ***")
print("*** ------------------------------------ ***\n")

# Input species parameters
config = Config()
config.parseConfig()

speciesFile = config.speciesFile
rxnFile = config.rxnFile
simFile = config.simFile
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
start = time.perf_counter()

excludedSpecies = ["C", "C#C", "C#CC#C", "C1=CC1CC2C=C2", "C1=CC=CC=C1", "C=1=C=C=1", "C=C", "C=CC=C", "C=CCCC=C", "CC", "CC(C)(C)C", "CC(C)CC(C)(C)C", "CC1=CC=CC=C1",
                   "[CH3]", "[CH]1C=C1", "[CH]C[CH]", "[C]1C=C1", "[C]=C[CH]", "[C]C=[C]", "[C]C[CH]"]

excludedSpecies_2 = ["C", "C#C", "C1=CC=CC=C1", "C=C", "C=C=C", "C=CC=C", "C=CCCC=C", "CC", "[CH3]"]

sz = len(keys)

oldExpts = expts.copy()

usedSpecies = ['CC', 'C=CC', 'C#CC', 'C1CCCCC1', 'C1=CC=CC=C1']

oldKeys = keys[:]

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

        if config.stochastic:
            rxnCnt = generate_rxns_stochastic(keys, s1, m1, rxnFile, config)
        else:
            rxnCnt = generate_rxns_complete(keys, s1, m1, rxnFile, config)

        if rxnCnt > 0:
            reactions = checkSimilarity(rxnCnt, rxnFile=rxnFile, simFile=simFile, config=config)
            datasets.append(analyzeReactions(resFile=resFile, r1=s1, reactions=reactions, inputCalcData=inputCalcData, uncerts=uncerts, expts=expts, calcs=calcs, config=config))

            # expts[s1] = datasets[-1][3]

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
