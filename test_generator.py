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

test_expts = {}

with open('test.csv') as f:
    for line in f:

        line = line.rstrip("\n")

        if line == "":
            continue

        toks = line.split(",")
        test_expts[toks[0]] = float(toks[1])
        # species[toks[0]] = 1
        expts[toks[0]] = float(toks[1])
        uncerts[toks[0]] = 1
        calcs[toks[0]] = float(toks[2])


keys = sorted(species.keys())

print(keys)

print("Number of species in the database:", len(keys))

# ------------------------------------
start = time.clock()

sz = len(keys)

def generateIsodesmic(s1, keys, config):
    inputCalcData = calcs[s1]

    m1 = Chem.MolFromSmiles(s1)

    if config.stochastic:
        rxnCnt = generate_rxns_stochastic(keys, s1, m1, rxnFile, config)
    else:
        rxnCnt = generate_rxns_complete(keys, s1, m1, rxnFile, config)

    if rxnCnt > 0:
        reactions = checkSimilarity(rxnCnt, rxnFile=rxnFile, simFile=simFile, config=config)
        res = analyzeReactions(resFile=resFile, r1=s1, reactions=reactions, inputCalcData=inputCalcData,
                                  uncerts=uncerts, expts=expts, calcs=calcs, config=config)

        return res[3]

    return 0

test_matrix = []

test_keys = test_expts.keys()

for s1 in test_keys:

    # if s1 != "C=1=C=C=1":
    #     continue

    print("\n\n###########################################\n")
    print("Considering", s1)

    curRow = [s1, expts[s1], calcs[s1]]

    config.featureMethod = "HRCF"
    config.radius = 4
    config.top = 0.1
    curRow.append(generateIsodesmic(s1, keys, config))

    # Isogyric
    config.featureMethod = "ECFP"
    config.radius = 0
    config.top = -1
    config.ePair = True
    curRow.append(generateIsodesmic(s1, keys, config))

    # Isodesmic
    # config.ePair = False
    config.hydroBond = True
    config.normalBond = True
    # config.radical = True
    # config.bondRadical = True
    curRow.append(generateIsodesmic(s1, keys, config))

    # Hypohomodesmotic
    # config.hydroBond = False
    # config.normalBond = False
    config.atomHybridization = True
    config.atomWithHydro = True
    curRow.append(generateIsodesmic(s1, keys, config))

    # Homodesmotic
    config.atomWithHydro = False
    # config.atomHybridization = False
    config.atomHybridWithHydro = True
    config.bondHybrid = True
    curRow.append(generateIsodesmic(s1, keys, config))

    # Hyperhomodesmotic
    # config.hybridization = False
    config.bondHydro = True
    curRow.append(generateIsodesmic(s1, keys, config))

    config.ePair = False
    config.atomHybridization = False
    config.atomHybridWithHydro = False
    config.hydroBond = False
    config.normalBond = False
    # config.radical = False
    # config.bondRadical = False
    config.bondHybrid = False
    config.bondHydro = False

    test_matrix.append(curRow)

print(test_matrix)

import pandas as pd
results = pd.DataFrame(np.asarray(test_matrix))
# results.to_csv('testset_same_CBS.csv', index=False)
results.to_csv('testset_same_AM1.csv', index=False)

# --------------------------------------------------
print("Execution time:", (time.clock() - start), "s.")

raw_input("Press Enter to continue...")
