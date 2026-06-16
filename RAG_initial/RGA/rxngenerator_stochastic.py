from utility_module import defaultdict, combine2Dicts, compare2Dicts, multiplyCoef, Chem, permutations, randrange
from balancing_module import BalanceEq
from bond_finder import findBonds
from rxn_writer import printResultDiff, printResultSame

def generate_same_rxn(keys, m1, r1, mols, atomR1, atomList, bondR1, bonds, rSize, pSize, maxCnt, fout):

    sCnt = 0

    reactants, products, rKeys, pKeys, atomBags, bondList = generate_permutations(keys=keys, m1=m1, r1=r1, mols=mols,
                                                                         atomR1=atomR1, atomList=atomList,
                                                                         bonds=bonds, rSize=rSize, pSize=pSize, maxCnt=maxCnt)

    rLen = len(reactants)
    rSize += 1

    for index in range(rLen):

        curRKeys = rKeys[index]
        curPKeys = pKeys[index]

        curAtomList = atomBags[index]
        curBonds = bondList[index]

        rList = {}

        for rIndex in range(rSize):
            rList = combine2Dicts(rList, curAtomList[rIndex])

        pList = {}
        for pIndex in range(pSize):
            pList = combine2Dicts(pList, curAtomList[pIndex + rSize])

        if compare2Dicts(rList, pList):

            rBonds = bondR1.copy()

            for rIndex in range(rSize - 1):
                rBonds = combine2Dicts(rBonds, curBonds[rIndex])

            pBonds = {}

            for pIndex in range(pSize):
                pBonds = combine2Dicts(pBonds, curBonds[pIndex + rSize - 1])

            if compare2Dicts(rBonds, pBonds):
                printResultSame(reactants=curRKeys, products=curPKeys, fout=fout)
                sCnt += 1

    return sCnt

def generate_diff_rxn(keys, m1, r1, mols, atomR1, atomList, bondR1, bonds, rSize, pSize, maxCnt, fout):

    sCnt = 0

    reactants, products, rKeys, pKeys, atomBags, bondList = generate_permutations(keys=keys, m1=m1, r1=r1, mols=mols,
                                                                         atomR1=atomR1, atomList=atomList,
                                                                         bonds=bonds, rSize=rSize, pSize=pSize, maxCnt=maxCnt)

    rLen = len(reactants)
    rSize += 1

    for index in range(rLen):

        curReactants = reactants[index]
        curProducts = products[index]

        curRKeys = rKeys[index]
        curPKeys = pKeys[index]

        curAtomList = atomBags[index]
        curBonds = bondList[index]

        coefList = BalanceEq().balanceEq(reactants=curReactants, products=curProducts, atomList=curAtomList)

        if coefList:

            for coefs in coefList:

                rBonds = multiplyCoef(bondR1, coefs[0])
                for rIndex in range(rSize - 1):
                    rBonds = combine2Dicts(rBonds, multiplyCoef(curBonds[rIndex], coefs[rIndex + 1]))

                pBonds = {}
                for pIndex in range(pSize):
                    pBonds = combine2Dicts(pBonds, multiplyCoef(curBonds[pIndex + rSize - 1], coefs[pIndex + rSize]))

                if compare2Dicts(rBonds, pBonds):
                    printResultDiff(reactants=curRKeys, products=curPKeys, coefs=coefs, fout=fout)
                    sCnt += 1

    return sCnt

def processSame_fast_3(r1, m1, keys, fout, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomList = {}

    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 1
    return generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, fout=fout)

def processSame_fast_4(r1, m1, keys, fout, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomList = {}

    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 2
    return generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                             bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, fout=fout)


def processSame_fast_5(r1, m1, keys, fout, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomList = {}

    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 3
    return generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                             bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, fout=fout)

def processSame_fast_6(r1, m1, keys, fout, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomList = {}

    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 3 - 3
    return generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                             bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, fout=fout)


def processSame(r1, m1, keys, fout, config):

    sCnt = 0

    noSpecies = config.noSpecies

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomList = {}

    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 1

    sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, fout=fout)


    if noSpecies > 3:

        # 2 - 2


        # start2 = time.clock()

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, fout=fout)

        # print "Time needed to generate rxns:", (time.clock() - start2), "s."


        # -----------------

        # 3 - 1
        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=1, maxCnt=maxCnt, fout=fout)

    if noSpecies > 4:

        # 2 - 3
        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=3, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 3 - 2

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 4 - 1

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=1, maxCnt=maxCnt, fout=fout)

    if noSpecies > 5:

        # 2 - 4
        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=4, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 3 - 3

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 4 - 2

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=2, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 5 - 1

        sCnt += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=4, pSize=1, maxCnt=maxCnt, fout=fout)

    return sCnt

def generate_permutations(keys, m1, r1, mols, atomR1, atomList, bonds, rSize=1, pSize=2, maxCnt=10):

    arr = keys[:]

    if rSize == 0 and pSize == 0:
        return

    if r1 in arr:
        arr.remove(r1)

    rList = []
    pList = []
    rKeys = []
    pKeys = []
    bondList = []
    atomBags = []

    dups = {}

    totalSz = rSize + pSize
    arrSz = len(arr)

    sNum = arrSz - totalSz + 1

    totalPermutations = 1

    for i in range(sNum, arrSz + 1):
        totalPermutations *= i

    if maxCnt > totalPermutations or maxCnt < 1:
        maxCnt = totalPermutations

    for cnt in range(maxCnt):

        remainingArr = arr[:]
        curReactants = []
        reactantsKey = []

        for index in range(rSize):
            chosenIndex = randrange(0, len(remainingArr))
            curReactants.append(mols[remainingArr[chosenIndex]])
            reactantsKey.append(remainingArr[chosenIndex])
            remainingArr.pop(chosenIndex)

        curProducts = []
        productsKey = []
        for index in range(pSize):
            chosenIndex = randrange(0, len(remainingArr))
            curProducts.append(mols[remainingArr[chosenIndex]])
            productsKey.append(remainingArr[chosenIndex])
            remainingArr.pop(chosenIndex)

        combinedKeys = reactantsKey + productsKey

        curKey = ','.join(combinedKeys)

        if curKey not in dups:

            curAtoms = [atomR1]
            curAtoms += [atomList[s] for s in combinedKeys]
            atomBags.append(curAtoms)

            curBonds = [bonds[s] for s in combinedKeys]
            bondList.append(curBonds)

            rPermutations = list(permutations(reactantsKey))
            pPermutations = list(permutations(productsKey))

            for rPermutation in rPermutations:
                for pPermutation in pPermutations:
                    curCombinedKeys = ','.join(rPermutation + pPermutation)
                    dups[curCombinedKeys] = 1

            reactantsKey = [r1] + reactantsKey
            curReactants = [m1] + curReactants

            rList.append(curReactants)
            pList.append(curProducts)

            rKeys.append(reactantsKey)
            pKeys.append(productsKey)

    return rList, pList, rKeys, pKeys, atomBags, bondList

def processDiff(r1, m1, keys, fout, config):

    sCnt = 0

    noSpecies = config.noSpecies

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = {}

    # for s in range(sz):
    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 1

    sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, fout=fout)


    if noSpecies > 3:

        # 2 - 2

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 3 - 1

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=1, maxCnt=maxCnt, fout=fout)

    # cnt = 0

    if noSpecies > 4:

        print("Starting with 2 - 3")

        # 2 - 3

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=3, maxCnt=maxCnt, fout=fout)

        # 3 - 2

        print("Starting with 3 - 2")

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, fout=fout)

        # 4 - 1

        print("Starting with 4 - 1")

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=1, maxCnt=maxCnt, fout=fout)

    if noSpecies > 5:

        print("Starting with 6 species")

        # 2 - 4

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=4, maxCnt=maxCnt, fout=fout)

        # 3 - 3

        print("Starting with 3 - 3")

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, fout=fout)

        # -----------------

        # 4 - 2

        print("Starting with 4 - 2")

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=2, maxCnt=maxCnt, fout=fout)

        # 5 - 1

        print("Starting with 5 - 1")

        sCnt += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=4, pSize=1, maxCnt=maxCnt, fout=fout)

    return sCnt

def processDiff_fast_3(r1, m1, keys, fout, config):

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = {}

    # for s in range(sz):
    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 1

    return generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, fout=fout)

def processDiff_fast_4(r1, m1, keys, fout, config):

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = {}

    # for s in range(sz):
    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 2 - 2

    return generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, fout=fout)

def processDiff_fast_5(r1, m1, keys, fout, config):

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = {}

    # for s in range(sz):
    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 3 - 2

    return generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, fout=fout)

def processDiff_fast_6(r1, m1, keys, fout, config):

    maxCnt = config.noTrials

    mols = {key: Chem.MolFromSmiles(key) for key in keys}
    bonds = {key: findBonds(Chem.MolFromSmiles(key), coef=1, config=config) for key in keys}

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = {}

    # for s in range(sz):
    for key in keys:
        m = mols[key]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList[key] = atoms

    bondR1 = findBonds(m1, coef=1, config=config)

    # 3 - 3

    return generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, fout=fout)

def generate_rxns_stochastic(keys, s1, m1, rxnFile, config):

    # The rxn generator parameters
    sameCoef, noSpecies, completeSet = config.sameCoef, config.noSpecies, config.completeSet

    print("Start generating reactions...")

    with open(rxnFile, "w") as fout:

        if sameCoef:
            if completeSet:
                sCnt = processSame(s1, m1, keys, fout, config)
            else:

                if noSpecies > 6:
                    noSpecies = 6

                fastSameProcessors = [processSame_fast_3, processSame_fast_4, processSame_fast_5, processSame_fast_6]

                index = noSpecies - 3

                sCnt = fastSameProcessors[index](s1, m1, keys, fout, config)
        else:
            if completeSet:
                sCnt = processDiff(s1, m1, keys, fout, config)
            else:
                if noSpecies > 6:
                    noSpecies = 6

                fastDiffProcessors = [processDiff_fast_3, processDiff_fast_4, processDiff_fast_5, processDiff_fast_6]

                index = noSpecies - 3

                sCnt = fastDiffProcessors[index](s1, m1, keys, fout, config)

        if sCnt == 0:
            fout.write("There is no reaction found for species " + s1 + "\n")

    print("End generating reactions...")
    print("Total isodesmic reactions found:", str(sCnt))
    print("--------------------------------\n")

    return sCnt
