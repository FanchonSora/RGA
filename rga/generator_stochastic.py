from .utils import defaultdict, combine2Dicts, compare2Dicts, multiplyCoef, Chem, permutations, randrange, \
    Pool, partial, adjustNoProcessors, divide_work, np, seed
from .balancing import BalanceEq
from .bond_finder import findBonds
from .rxn_writer import printResultDiffParallel, printResultSameParallel

def parallel_same_stoch(curRKeys, curPKeys, curAtomList, bondR1, curBonds, rSize, pSize):
# def parallel_same_stoch(rKeys, pKeys, atomBags, bondR1, bondList, rSize, pSize):

    # rxns = []

    rxn = None

    # sz = len(rKeys)
    #
    # for index in range(sz):
    #
    #     curRKeys, curPKeys, curAtomList, curBonds = rKeys[index], pKeys[index], atomBags[index], bondList[index]

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
            rxn = printResultSameParallel(reactants=curRKeys, products=curPKeys)
            # rxns.append(printResultSameParallel(reactants=curRKeys, products=curPKeys))

    return rxn

def generate_same_rxn(keys, m1, r1, mols, atomR1, atomList, bondR1, bonds, rSize, pSize, maxCnt, config):

    reactants, products, rKeys, pKeys, atomBags, bondList = generate_permutations(keys=keys, m1=m1, r1=r1, mols=mols,
                                                                         atomR1=atomR1, atomList=atomList,
                                                                         bonds=bonds, config=config, rSize=rSize, pSize=pSize,
                                                                                  maxCnt=maxCnt)

    rLen = len(reactants)
    rSize += 1

    # print("Same len:", rLen)

    # args = [[rKeys[index], pKeys[index], atomBags[index], bondR1, bondList[index], rSize, pSize] for index in range(rLen)]
    #
    # no_processors = adjustNoProcessors(config.noProcessors)
    #
    # mp = Pool(processes=no_processors)
    #
    # func = partial(divide_work, parallel_same_stoch)
    #
    # import time
    # start = time.clock()
    #
    # rxnList = mp.imap_unordered(func, args, chunksize=1)
    #
    # mp.close()
    # mp.join()
    #
    # print time.clock() - start, "s."
    #
    # rxnList = [rxn for rxn in rxnList if rxn]

    rxnList = []

    for index in range(rLen):
        curRKeys, curPKeys, curAtomList, curBonds = rKeys[index], pKeys[index], atomBags[index], bondList[index]

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
                # rxn = printResultSameParallel(reactants=curRKeys, products=curPKeys)
                rxnList.append(printResultSameParallel(reactants=curRKeys, products=curPKeys))

    return rxnList

def parallel_diff_stoch(curReactants, curProducts, curRKeys, curPKeys, curAtomList, bondR1, curBonds, rSize, pSize):

    coefList = BalanceEq().balanceEq(reactants=curReactants, products=curProducts, atomList=curAtomList)

    # print curRKeys, curPKeys

    rxns = []

    if coefList:

        for coefs in coefList:

            rBonds = multiplyCoef(bondR1, coefs[0])
            for rIndex in range(rSize - 1):
                rBonds = combine2Dicts(rBonds, multiplyCoef(curBonds[rIndex], coefs[rIndex + 1]))

            pBonds = {}
            for pIndex in range(pSize):
                pBonds = combine2Dicts(pBonds, multiplyCoef(curBonds[pIndex + rSize - 1], coefs[pIndex + rSize]))

            if compare2Dicts(rBonds, pBonds):
                rxns.append(printResultDiffParallel(reactants=curRKeys, products=curPKeys, coefs=coefs))

    return rxns

def generate_diff_rxn(keys, m1, r1, mols, atomR1, atomList, bondR1, bonds, rSize, pSize, maxCnt, config):

    reactants, products, rKeys, pKeys, atomBags, bondList = generate_permutations(keys=keys, m1=m1, r1=r1, mols=mols,
                                                                         atomR1=atomR1, atomList=atomList,
                                                                         bonds=bonds, config=config, rSize=rSize, pSize=pSize,
                                                                                  maxCnt=maxCnt)

    rLen = len(reactants)
    rSize += 1

    # print("Diff len:", rLen)

    if config.noSpecies > 4 or config.noTrials > 1000:

        args = [[reactants[index], products[index], rKeys[index], pKeys[index], atomBags[index], bondR1, bondList[index], rSize, pSize]
                for index in range(rLen)]

        no_processors = adjustNoProcessors(config.noProcessors)

        mp = Pool(processes=no_processors)

        func = partial(divide_work, parallel_diff_stoch)

        chunkSz = rLen // no_processors

        tmpList = mp.imap_unordered(func, args, chunksize=chunkSz)

        mp.close()
        mp.join()

        rxnList = []

        for rxns in tmpList:
            if rxns and len(rxns) > 0:
                rxnList += rxns

    else:

        rxnList = []

        for index in range(rLen):
            curReactants, curProducts, curRKeys, curPKeys, curAtomList, curBonds = reactants[index], products[index], rKeys[index], pKeys[index],\
                                                                                   atomBags[index], bondList[index]

            coefList = BalanceEq().balanceEq(reactants=curReactants, products=curProducts, atomList=curAtomList)

            # rxn = None

            if coefList:

                for coefs in coefList:

                    rBonds = multiplyCoef(bondR1, coefs[0])
                    for rIndex in range(rSize - 1):
                        rBonds = combine2Dicts(rBonds, multiplyCoef(curBonds[rIndex], coefs[rIndex + 1]))

                    pBonds = {}
                    for pIndex in range(pSize):
                        pBonds = combine2Dicts(pBonds, multiplyCoef(curBonds[pIndex + rSize - 1], coefs[pIndex + rSize]))

                    if compare2Dicts(rBonds, pBonds):
                        # rxn = printResultDiffParallel(reactants=curRKeys, products=curPKeys, coefs=coefs)
                        rxnList.append(printResultDiffParallel(reactants=curRKeys, products=curPKeys, coefs=coefs))

    return rxnList

def processSame_fast_3(r1, m1, keys, config):

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
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, config=config)

def processSame_fast_4(r1, m1, keys, config):

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
                             bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, config=config)

def processSame_fast_5(r1, m1, keys, config):

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
                             bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, config=config)

def processSame_fast_6(r1, m1, keys, config):

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
                             bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, config=config)

def processSame(r1, m1, keys, config):

    rxnList = []

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

    rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, config=config)

    # print("Current 2-1 len:", len(rxnList))

    if noSpecies > 3:

        # 2 - 2
        # start2 = time.clock()

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, config=config)

        # print("Current 2-2 len:", len(rxnList))

        # print "Time needed to generate rxns:", (time.clock() - start2), "s."


        # -----------------

        # 3 - 1
        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=1, maxCnt=maxCnt, config=config)

        # print("Current 3-1 len:", len(rxnList))

    if noSpecies > 4:

        # 2 - 3
        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=3, maxCnt=maxCnt, config=config)

        # -----------------

        # 3 - 2

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, config=config)

        # -----------------

        # 4 - 1

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=1, maxCnt=maxCnt, config=config)

    if noSpecies > 5:

        # 2 - 4
        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=4, maxCnt=maxCnt, config=config)

        # -----------------

        # 3 - 3

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, config=config)

        # -----------------

        # 4 - 2

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=2, maxCnt=maxCnt, config=config)

        # -----------------

        # 5 - 1

        rxnList += generate_same_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=4, pSize=1, maxCnt=maxCnt, config=config)

    return rxnList

def generate_permutations(keys, m1, r1, mols, atomR1, atomList, bonds, config, rSize=1, pSize=2, maxCnt=10):

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

    # print(totalPermutations)

    cur_seed = config.seed

    check_seed = True

    if cur_seed == None:
        check_seed = False

    for cnt in range(maxCnt):

        remainingArr = arr[:]
        curReactants = []
        reactantsKey = []

        for index in range(rSize):

            # Deterministic seed
            if check_seed:
                seed(cur_seed)
                cur_seed += 1

            chosenIndex = randrange(0, len(remainingArr))
            curReactants.append(mols[remainingArr[chosenIndex]])
            reactantsKey.append(remainingArr[chosenIndex])
            remainingArr.pop(chosenIndex)

        curProducts = []
        productsKey = []
        for index in range(pSize):

            # Deterministic seed
            if check_seed:
                seed(cur_seed)
                cur_seed += 1

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

    # print(len(dups), len(rKeys))

    # print("Max cnt", maxCnt, "seed", cur_seed)

    return rList, pList, rKeys, pKeys, atomBags, bondList

def processDiff(r1, m1, keys, config):

    rxnList = []

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

    rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, config=config)

    # print("Current len 2-1:", len(rxnList))


    if noSpecies > 3:

        # 2 - 2

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, config=config)

        # print("Current len 2-2:", len(rxnList))

        # -----------------

        # 3 - 1

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=1, maxCnt=maxCnt, config=config)

        # print("Current len 3-1:", len(rxnList))

    # cnt = 0

    if noSpecies > 4:

        # print("Starting with 2 - 3")

        # 2 - 3

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=3, maxCnt=maxCnt, config=config)

        # print("Current len 2-3:", len(rxnList))

        # 3 - 2

        # print("Starting with 3 - 2")

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, config=config)

        # print("Current len 3-2:", len(rxnList))

        # 4 - 1

        # print("Starting with 4 - 1")

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=1, maxCnt=maxCnt, config=config)

        # print("Current len 4-1:", len(rxnList))

    if noSpecies > 5:

        # print("Starting with 6 species")

        # 2 - 4

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=1, pSize=4, maxCnt=maxCnt, config=config)

        # 3 - 3

        # print("Starting with 3 - 3")

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, config=config)

        # -----------------

        # 4 - 2

        # print("Starting with 4 - 2")

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=3, pSize=2, maxCnt=maxCnt, config=config)

        # 5 - 1

        # print("Starting with 5 - 1")

        rxnList += generate_diff_rxn(keys=keys, m1=m1, r1=r1, mols=mols, atomR1=atomR1, atomList=atomList,
                                  bondR1=bondR1, bonds=bonds, rSize=4, pSize=1, maxCnt=maxCnt, config=config)

    return rxnList

def processDiff_fast_3(r1, m1, keys, config):

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
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=1, maxCnt=maxCnt, config=config)

def processDiff_fast_4(r1, m1, keys, config):

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
                              bondR1=bondR1, bonds=bonds, rSize=1, pSize=2, maxCnt=maxCnt, config=config)

def processDiff_fast_5(r1, m1, keys, config):

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
                              bondR1=bondR1, bonds=bonds, rSize=2, pSize=2, maxCnt=maxCnt, config=config)


def processDiff_fast_6(r1, m1, keys, config):

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
                              bondR1=bondR1, bonds=bonds, rSize=2, pSize=3, maxCnt=maxCnt, config=config)

def generate_rxns_stochastic(keys, s1, m1, fout, config):

    # The rxn generator parameters
    sameCoef, noSpecies, completeSet = config.sameCoef, config.noSpecies, config.completeSet

    print("Start generating reactions...")

    if sameCoef:
        if completeSet:
            rxnList = processSame(s1, m1, keys, config)
        else:

            if noSpecies > 6:
                noSpecies = 6

            fastSameProcessors = [processSame_fast_3, processSame_fast_4, processSame_fast_5, processSame_fast_6]

            index = noSpecies - 3

            rxnList = fastSameProcessors[index](s1, m1, keys, config)
    else:
        if completeSet:
            rxnList = processDiff(s1, m1, keys, config)
        else:
            if noSpecies > 6:
                noSpecies = 6

            fastDiffProcessors = [processDiff_fast_3, processDiff_fast_4, processDiff_fast_5, processDiff_fast_6]

            index = noSpecies - 3

            rxnList = fastDiffProcessors[index](s1, m1, keys, config)

    sCnt = len(rxnList)

    if not config.similarityOn:

        fout.write("SUMMARY\n\n")
        fout.write("Total number of reactions found: " + str(sCnt) + "\n\n")
        fout.write("=========================================================================\n")

        if sCnt > 0:

            maxWidth = 15 * config.noSpecies
            maxCnt = len(str(sCnt))

            fout.write("LIST OF REACTIONS FOUND\n\n")
            fout.write("%*s\t\t%*s\n\n" % (maxCnt, "No.", maxWidth, "Reaction"))

            for index in range(sCnt):

                line = rxnList[index]
                fout.write("%*d\t\t%*s\n" % (maxCnt, index + 1, maxWidth, line))

                fout.write("\n=========================================================================\n")

    print("End generating reactions...")
    print("Total isodesmic reactions found: " + str(sCnt))
    print("--------------------------------\n")

    return sCnt, rxnList
