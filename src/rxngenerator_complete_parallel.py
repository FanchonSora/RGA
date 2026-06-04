from utility_module import defaultdict, combine2Dicts, compare2Dicts, multiplyCoef, Chem, Pool, partial, adjustNoProcessors, divide_work, np
from balancing_module import BalanceEq
from bond_finder import findBonds
from rxn_writer import printResultSameParallel, printResultDiffParallel

def parallel_same_fast_3(r2, r1, keys, atomR1, atomList, bonds, bondR1):

    if keys[r2] == r1:
        return None

    rxns = []

    rList = combine2Dicts(atomR1, atomList[r2])

    bondR2 = bonds[r2]
    rBonds = combine2Dicts(bondR1, bondR2)

    sz = len(keys)

    for p1 in range(sz):
        if p1 == r2 or keys[p1] == r1:
            continue

        pList = atomList[p1]
        bondP1 = bonds[p1]

        if compare2Dicts(rList, pList) and compare2Dicts(rBonds, bondP1):
            rxns.append(printResultSameParallel(reactants=[r1, keys[r2]], products=[keys[p1]]))

    return rxns

def processSame_fast_3(r1, m1, keys, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_same_fast_3)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_same_fast_4(r2, r1, keys, atomR1, atomList, bonds, bondR1):

    if keys[r2] == r1:
        return None

    rList = combine2Dicts(atomR1, atomList[r2])

    bondR2 = bonds[r2]
    rBonds = combine2Dicts(bondR1, bondR2)

    rxns = []

    sz = len(keys)

    for p1 in range(sz - 1):
        if p1 == r2 or keys[p1] == r1:
            continue

        bondP1 = bonds[p1]

        for p2 in range(p1 + 1, sz):

            if p2 == r2 or keys[p2] == r1:
                continue

            pList = combine2Dicts(atomList[p1], atomList[p2])

            if compare2Dicts(rList, pList):

                bondP2 = bonds[p2]
                pBonds = combine2Dicts(bondP1, bondP2)

                if compare2Dicts(rBonds, pBonds):
                    rxns.append(printResultSameParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2]]))

    return rxns

def processSame_fast_4(r1, m1, keys, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_same_fast_4)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_same_fast_5(r2, r1, keys, atomR1, atomList, bonds, bondR1):

    if keys[r2] == r1:
        return None

    bondR2 = bonds[r2]

    sz = len(keys)

    rxns = []

    for r3 in range(r2 + 1, sz):
        if keys[r3] == r1:
            continue

        rList = combine2Dicts(atomR1, atomList[r2])
        rList = combine2Dicts(rList, atomList[r3])

        bondR3 = bonds[r3]
        rBonds = combine2Dicts(bondR1, bondR2)
        rBonds = combine2Dicts(rBonds, bondR3)

        for p1 in range(sz - 1):

            if p1 == r2 or p1 == r3 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz):

                if p2 == r2 or p2 == r3 or keys[p2] == r1:
                    continue

                pList = combine2Dicts(atomList[p1], atomList[p2])

                if compare2Dicts(rList, pList):

                    bondP2 = bonds[p2]
                    pBonds = combine2Dicts(bondP1, bondP2)

                    if compare2Dicts(rBonds, pBonds):
                        rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1], keys[p2]]))

    return rxns

def processSame_fast_5(r1, m1, keys, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_same_fast_5)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_same_fast_6(r2, r1, keys, atomR1, atomList, bonds, bondR1):

    if keys[r2] == r1:
        return None

    rxns = []

    sz = len(keys)

    bondR2 = bonds[r2]

    for r3 in range(r2 + 1, sz):
        if keys[r3] == r1:
            continue

        rList = combine2Dicts(atomR1, atomList[r2])
        rList = combine2Dicts(rList, atomList[r3])

        bondR3 = bonds[r3]
        rBonds = combine2Dicts(bondR1, bondR2)
        rBonds = combine2Dicts(rBonds, bondR3)

        for p1 in range(sz - 2):

            if p1 == r2 or p1 == r3 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 1):

                if p2 == r2 or p2 == r3 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz):

                    if p3 == r2 or p3 == r3 or keys[p3] == r1:
                        continue

                    pList = combine2Dicts(atomList[p1], atomList[p2])
                    pList = combine2Dicts(pList, atomList[p3])

                    if compare2Dicts(rList, pList):

                        bondP3 = bonds[p3]
                        pBonds = combine2Dicts(bondP1, bondP2)
                        pBonds = combine2Dicts(pBonds, bondP3)

                        if compare2Dicts(rBonds, pBonds):
                            rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3]],
                                            products=[keys[p1], keys[p2], keys[p3]]))

    return rxns

def processSame_fast_6(r1, m1, keys, config):

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_same_fast_6)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_same(r2, r1, keys, noSpecies, atomR1, atomList, bonds, bondR1):

    sz = len(keys)

    if keys[r2] == r1:
        return None

    bondR2 = bonds[r2]
    rList = combine2Dicts(atomR1, atomList[r2])
    rBonds = combine2Dicts(bondR1, bondR2)

    rxns = []

    # print(keys[r2], bondR2)

    # 2 - 1
    for p1 in range(sz):
        if p1 == r2 or keys[p1] == r1:
            continue

        if compare2Dicts(rList, atomList[p1]):

            bondP1 = bonds[p1]
            # cnt += 1
            if compare2Dicts(rBonds, bondP1):
                rxns.append(printResultSameParallel(reactants=[r1, keys[r2]], products=[keys[p1]]))

    if noSpecies > 3:

        # 2 - 2
        for p1 in range(sz - 1):
            if p1 == r2 or keys[p1] == r1:
                continue

            # Checking atom types
            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz):

                if p2 == r2 or keys[p2] == r1:
                    continue

                pList = combine2Dicts(atomList[p1], atomList[p2])

                if compare2Dicts(rList, pList):

                    bondP2 = bonds[p2]
                    pBonds = combine2Dicts(bondP1, bondP2)
                    # cnt += 1
                    if compare2Dicts(rBonds, pBonds):
                        rxns.append(printResultSameParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2]]))

        # -----------------

        # 3 - 1
        for r3 in range(r2 + 1, sz):
            if keys[r3] == r1:
                continue

            rList = combine2Dicts(atomR1, atomList[r2])
            rList = combine2Dicts(rList, atomList[r3])

            bondR3 = bonds[r3]
            rBonds = combine2Dicts(bondR1, bondR2)
            rBonds = combine2Dicts(rBonds, bondR3)

            for p1 in range(sz):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                if compare2Dicts(rList, atomList[p1]):

                    bondP1 = bonds[p1]
                    # cnt += 1
                    if compare2Dicts(rBonds, bondP1):
                        rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1]]))

    if noSpecies > 4:

        # 2 - 3
        rList = combine2Dicts(atomR1, atomList[r2])
        rBonds = combine2Dicts(bondR1, bondR2)

        for p1 in range(sz - 2):
            if p1 == r2 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 1):

                if p2 == r2 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz):

                    if p3 == r2 or keys[p3] == r1:
                        continue

                    pList = combine2Dicts(atomList[p1], atomList[p2])
                    pList = combine2Dicts(pList, atomList[p3])

                    if compare2Dicts(rList, pList):

                        bondP3 = bonds[p3]
                        pBonds = combine2Dicts(bondP1, bondP2)
                        pBonds = combine2Dicts(pBonds, bondP3)
                        # cnt += 1
                        if compare2Dicts(rBonds, pBonds):
                            rxns.append(printResultSameParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2], keys[p3]]))

        # -----------------

        # 3 - 2

        for r3 in range(r2 + 1, sz):

            if keys[r3] == r1:
                continue

            rList = combine2Dicts(atomR1, atomList[r2])
            rList = combine2Dicts(rList, atomList[r3])

            bondR3 = bonds[r3]
            rBonds = combine2Dicts(bondR1, bondR2)
            rBonds = combine2Dicts(rBonds, bondR3)

            for p1 in range(sz - 1):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                bondP1 = bonds[p1]

                for p2 in range(p1 + 1, sz):

                    if p2 == r2 or p2 == r3 or keys[p2] == r1:
                        continue

                    pList = combine2Dicts(atomList[p1], atomList[p2])

                    if compare2Dicts(rList, pList):

                        bondP2 = bonds[p2]
                        pBonds = combine2Dicts(bondP1, bondP2)
                        # cnt += 1
                        if compare2Dicts(rBonds, pBonds):
                            rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1], keys[p2]]))

        # -----------------

        # 4 - 1

        for r3 in range(r2 + 1, sz - 1):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz):

                if keys[r4] == r1:
                    continue

                rList = combine2Dicts(atomR1, atomList[r2])
                rList = combine2Dicts(rList, atomList[r3])
                rList = combine2Dicts(rList, atomList[r4])

                bondR4 = bonds[r4]
                rBonds = combine2Dicts(bondR1, bondR2)
                rBonds = combine2Dicts(rBonds, bondR3)
                rBonds = combine2Dicts(rBonds, bondR4)

                for p1 in range(sz):

                    if p1 == r2 or p1 == r3 or p1 == r4 or keys[p1] == r1:
                        continue

                    if compare2Dicts(rList, atomList[p1]):

                        bondP1 = bonds[p1]
                        # cnt += 1
                        if compare2Dicts(rBonds, bondP1):
                            rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3], keys[r4]], products=[keys[p1]]))

    if noSpecies > 5:

        # 2 - 4
        rList = combine2Dicts(atomR1, atomList[r2])
        rBonds = combine2Dicts(bondR1, bondR2)

        for p1 in range(sz - 3):

            if p1 == r2 or keys[p1] == r1:
                continue

            # Checking atom types
            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 2):

                if p2 == r2 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz - 1):

                    if p3 == r2 or keys[p3] == r1:
                        continue

                    bondP3 = bonds[p3]

                    for p4 in range(p3 + 1, sz):

                        if p4 == r2 or keys[p4] == r1:
                            continue

                        pList = combine2Dicts(atomList[p1], atomList[p2])
                        pList = combine2Dicts(pList, atomList[p3])
                        pList = combine2Dicts(pList, atomList[p4])

                        if compare2Dicts(rList, pList):

                            bondP4 = bonds[p4]

                            pBonds = combine2Dicts(bondP1, bondP2)
                            pBonds = combine2Dicts(pBonds, bondP3)
                            pBonds = combine2Dicts(pBonds, bondP4)

                            # cnt += 1

                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultSameParallel(reactants=[r1, keys[r2]],
                                                             products=[keys[p1], keys[p2], keys[p3], keys[p4]]))

        # -----------------

        # 3 - 3

        for r3 in range(r2 + 1, sz):

            if keys[r3] == r1:
                continue

            rList = combine2Dicts(atomR1, atomList[r2])
            rList = combine2Dicts(rList, atomList[r3])

            bondR3 = bonds[r3]
            rBonds = combine2Dicts(bondR1, bondR2)
            rBonds = combine2Dicts(rBonds, bondR3)

            for p1 in range(sz - 2):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                bondP1 = bonds[p1]

                for p2 in range(p1 + 1, sz - 1):

                    if p2 == r2 or p2 == r3 or keys[p2] == r1:
                        continue

                    bondP2 = bonds[p2]

                    for p3 in range(p2 + 1, sz):

                        if p3 == r2 or p3 == r3 or keys[p3] == r1:
                            continue

                        pList = combine2Dicts(atomList[p1], atomList[p2])
                        pList = combine2Dicts(pList, atomList[p3])

                        if compare2Dicts(rList, pList):

                            bondP3 = bonds[p3]

                            pBonds = combine2Dicts(bondP1, bondP2)
                            pBonds = combine2Dicts(pBonds, bondP3)
                            # cnt += 1
                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3]],
                                                             products=[keys[p1], keys[p2], keys[p3]]))

        # -----------------

        # 4 - 2

        for r3 in range(r2 + 1, sz - 1):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz):

                if keys[r4] == r1:
                    continue

                rList = combine2Dicts(atomR1, atomList[r2])
                rList = combine2Dicts(rList, atomList[r3])
                rList = combine2Dicts(rList, atomList[r4])

                bondR4 = bonds[r4]
                rBonds = combine2Dicts(bondR1, bondR2)
                rBonds = combine2Dicts(rBonds, bondR3)
                rBonds = combine2Dicts(rBonds, bondR4)

                for p1 in range(sz - 1):

                    if p1 == r2 or p1 == r3 or p1 == r4 or keys[p1] == r1:
                        continue

                    bondP1 = bonds[p1]

                    for p2 in range(p1 + 1, sz):

                        if p2 == r2 or p2 == r3 or p2 == r4 or keys[p2] == r1:
                            continue

                        pList = combine2Dicts(atomList[p1], atomList[p2])

                        if compare2Dicts(rList, pList):

                            bondP2 = bonds[p2]
                            pBonds = combine2Dicts(bondP1, bondP2)

                            # cnt += 1

                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3], keys[r4]],
                                                             products=[keys[p1], keys[p2]]))

        # -----------------

        # 5 - 1

        for r3 in range(r2 + 1, sz - 2):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz - 1):

                if keys[r4] == r1:
                    continue

                bondR4 = bonds[r4]

                for r5 in range(r4 + 1, sz):

                    if keys[r5] == r1:
                        continue

                    rList = combine2Dicts(atomR1, atomList[r2])
                    rList = combine2Dicts(rList, atomList[r3])
                    rList = combine2Dicts(rList, atomList[r4])
                    rList = combine2Dicts(rList, atomList[r5])

                    bondR5 = bonds[r5]
                    rBonds = combine2Dicts(bondR1, bondR2)
                    rBonds = combine2Dicts(rBonds, bondR3)
                    rBonds = combine2Dicts(rBonds, bondR4)
                    rBonds = combine2Dicts(rBonds, bondR5)

                    for p1 in range(sz):

                        if p1 == r2 or p1 == r3 or p1 == r4 or p1 == r5 or keys[p1] == r1:
                            continue

                        if compare2Dicts(rList, atomList[p1]):

                            bondP1 = bonds[p1]
                            # cnt += 1
                            if compare2Dicts(rBonds, bondP1):
                                rxns.append(printResultSameParallel(reactants=[r1, keys[r2], keys[r3], keys[r4], keys[r5]], products=[keys[p1]]))

    return rxns

def processSame(r1, m1, keys, config):

    noSpecies = config.noSpecies

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, noSpecies, atomR1, atomList, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_same)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_diff(r2, r1, keys, noSpecies, atomR1, atomList, m1, mols, bonds, bondR1):

    if keys[r2] == r1:
        return None

    sz = len(keys)

    # print(keys[r2])

    bondR2 = bonds[r2]

    rxns = []

    # print("Starting with 2 - 1")

    for p1 in range(sz):
        if p1 == r2 or keys[p1] == r1:
            continue

        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1]],
                                         atomList=[atomR1, atomList[r2], atomList[p1]])
        # coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1]])

        if coefList:

            for coefs in coefList:

                bondP1 = bonds[p1]

                bR1_coef, bR2_coef, bP1_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(bondR2,
                                                                                            coefs[1]), multiplyCoef(
                    bondP1, coefs[2])
                rBonds = combine2Dicts(bR1_coef, bR2_coef)

                if compare2Dicts(rBonds, bP1_coef):
                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]], products=[keys[p1]], coefs=coefs))

    if noSpecies > 3:

        # 2 - 2
        # for r2 in range(sz):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]
        # print("Starting with 2 - 2")

        for p1 in range(sz - 1):
            if p1 == r2 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz):

                if p2 == r2 or keys[p2] == r1:
                    continue

                coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1], mols[p2]],
                                                 atomList=[atomR1, atomList[r2], atomList[p1], atomList[p2]])

                # coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1], mols[p2]])

                if coefList:

                    for coefs in coefList:

                        bondP2 = bonds[p2]

                        bR1_coef, bR2_coef, bP1_coef, bP2_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(
                            bondR2, coefs[1]), \
                                                                 multiplyCoef(bondP1, coefs[2]), multiplyCoef(
                            bondP2, coefs[3])

                        rBonds = combine2Dicts(bR1_coef, bR2_coef)
                        pBonds = combine2Dicts(bP1_coef, bP2_coef)

                        if compare2Dicts(rBonds, pBonds):
                            rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2]], coefs=coefs))

        # -----------------

        # 3 - 1

        # for r2 in range(sz):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        # print("Starting with 3 - 1")

        for r3 in range(r2 + 1, sz):
            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for p1 in range(sz):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]], products=[mols[p1]],
                                                 atomList=[atomR1, atomList[r2], atomList[r3], atomList[p1]])

                # coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]], products=[mols[p1]])

                if coefList:

                    for coefs in coefList:

                        bondP1 = bonds[p1]

                        bR1_coef, bR2_coef, bR3_coef, bP1_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                 multiplyCoef(bondR3, coefs[2]), multiplyCoef(bondP1, coefs[3])

                        rBonds = combine2Dicts(bR1_coef, bR2_coef)
                        rBonds = combine2Dicts(rBonds, bR3_coef)

                        if compare2Dicts(rBonds, bP1_coef):
                            rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1]], coefs=coefs))

    # cnt = 0

    if noSpecies > 4:

        print("Starting with 2 - 3")

        # 2 - 3

        # for r2 in range(sz):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     print(keys[r2])
        #
        #     bondR2 = bonds[r2]

        for p1 in range(sz - 2):
            if p1 == r2 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 1):

                if p2 == r2 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz):

                    if p3 == r2 or keys[p3] == r1:
                        continue

                    # cnt += 1

                    coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]],
                                                     products=[mols[p1], mols[p2], mols[p3]],
                                                     atomList=[atomR1, atomList[r2], atomList[p1], atomList[p2],
                                                               atomList[p3]])

                    if coefList:

                        for coefs in coefList:

                            bondP3 = bonds[p3]

                            bR1_coef, bR2_coef, bP1_coef, bP2_coef, bP3_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                               multiplyCoef(bondP1,coefs[2]), multiplyCoef(bondP2, coefs[3]),\
                                                                               multiplyCoef(bondP3, coefs[4])

                            rBonds = combine2Dicts(bR1_coef, bR2_coef)
                            pBonds = combine2Dicts(bP1_coef, bP2_coef)
                            pBonds = combine2Dicts(pBonds, bP3_coef)

                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2], keys[p3]], coefs=coefs))

        # print(cnt)

        # cnt = 0

        # -----------------

        # 3 - 2

        print("Starting with 3 - 2")

        # for r2 in range(sz - 1):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for r3 in range(r2 + 1, sz):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for p1 in range(sz - 1):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                bondP1 = bonds[p1]

                for p2 in range(p1 + 1, sz):

                    if p2 == r2 or p2 == r3 or keys[p2] == r1:
                        continue

                    # cnt += 1

                    coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]],
                                                     products=[mols[p1], mols[p2]],
                                                     atomList=[atomR1, atomList[r2], atomList[r3], atomList[p1],
                                                               atomList[p2]])

                    if coefList:

                        for coefs in coefList:

                            bondP2 = bonds[p2]

                            bR1_coef, bR2_coef, bR3_coef, bP1_coef, bP2_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                               multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondP1, coefs[3]),\
                                                                               multiplyCoef(bondP2, coefs[4])

                            rBonds = combine2Dicts(bR1_coef, bR2_coef)
                            rBonds = combine2Dicts(rBonds, bR3_coef)
                            pBonds = combine2Dicts(bP1_coef, bP2_coef)

                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1], keys[p2]], coefs=coefs))

        # print(cnt)

        # if r2 == 0:
        #     break

        # -----------------

        # 4 - 1

        # cnt = 0

        print("Starting with 4 - 1")

        # for r2 in range(sz - 2):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for r3 in range(r2 + 1, sz - 1):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz):

                if keys[r4] == r1:
                    continue

                bondR4 = bonds[r4]

                for p1 in range(sz):

                    if p1 == r2 or p1 == r3 or p1 == r4 or keys[p1] == r1:
                        continue

                    coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3], mols[r4]],
                                                     products=[mols[p1]],
                                                     atomList=[atomR1, atomList[r2], atomList[r3], atomList[r4],
                                                               atomList[p1]])

                    # cnt += 1

                    if coefList:

                        for coefs in coefList:

                            bondP1 = bonds[p1]

                            bR1_coef, bR2_coef, bR3_coef, bR4_coef, bP1_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                               multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondR4, coefs[3]),\
                                                                               multiplyCoef(bondP1, coefs[4])

                            rBonds = combine2Dicts(bR1_coef, bR2_coef)
                            rBonds = combine2Dicts(rBonds, bR3_coef)
                            rBonds = combine2Dicts(rBonds, bR4_coef)

                            if compare2Dicts(rBonds, bP1_coef):
                                rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3], keys[r4]], products=[keys[p1]], coefs=coefs))

    # print(cnt)

    # cnt = 0

    if noSpecies > 5:

        print("Starting with 2 - 4")

        # 2 - 4

        # for r2 in range(sz):
        #
        #     print(keys[r2])
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for p1 in range(sz - 3):

            if p1 == r2 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 2):

                if p2 == r2 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz - 1):

                    if p3 == r2 or keys[p3] == r1:
                        continue

                    bondP3 = bonds[p3]

                    for p4 in range(p3 + 1, sz):

                        if p4 == r2 or keys[p4] == r1:
                            continue

                        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]],
                                                         products=[mols[p1], mols[p2], mols[p3], mols[p4]],
                                                         atomList=[atomR1, atomList[r2], atomList[p1], atomList[p2],
                                                                   atomList[p3], atomList[p4]])

                        # cnt += 1

                        if coefList:

                            for coefs in coefList:

                                bondP4 = bonds[p4]

                                bR1_coef, bR2_coef, bP1_coef, bP2_coef, bP3_coef, bP4_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                                             multiplyCoef(bondP1,coefs[2]), multiplyCoef(bondP2, coefs[3]),\
                                                                                             multiplyCoef(bondP3,coefs[4]), multiplyCoef(bondP4, coefs[5])

                                rBonds = combine2Dicts(bR1_coef, bR2_coef)
                                pBonds = combine2Dicts(bP1_coef, bP2_coef)
                                pBonds = combine2Dicts(pBonds, bP3_coef)
                                pBonds = combine2Dicts(pBonds, bP4_coef)

                                if compare2Dicts(rBonds, pBonds):
                                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]],
                                                    products=[keys[p1], keys[p2], keys[p3], keys[p4]], coefs=coefs))
        #
        # print(cnt)

        # if r2 == 0:
        #     break

        # -----------------

        # 3 - 3

        print("Starting with 3 - 3")

        # for r2 in range(sz - 1):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for r3 in range(r2 + 1, sz):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for p1 in range(sz - 2):

                if p1 == r2 or p1 == r3 or keys[p1] == r1:
                    continue

                bondP1 = bonds[p1]

                for p2 in range(p1 + 1, sz - 1):

                    if p2 == r2 or p2 == r3 or keys[p2] == r1:
                        continue

                    bondP2 = bonds[p2]

                    for p3 in range(p2 + 1, sz):

                        if p3 == r2 or p3 == r3 or keys[p3] == r1:
                            continue

                        # cnt += 1

                        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]],
                                                         products=[mols[p1], mols[p2], mols[p3]],
                                                         atomList=[atomR1, atomList[r2], atomList[r3], atomList[p1],
                                                                   atomList[p2], atomList[p3]])

                        if coefList:

                            for coefs in coefList:

                                bondP3 = bonds[p3]

                                bR1_coef, bR2_coef, bR3_coef, bP1_coef, bP2_coef, bP3_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                                             multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondP1, coefs[3]),\
                                                                                             multiplyCoef(bondP2,coefs[4]), multiplyCoef(bondP3, coefs[5])

                                rBonds = combine2Dicts(bR1_coef, bR2_coef)
                                rBonds = combine2Dicts(rBonds, bR3_coef)
                                pBonds = combine2Dicts(bP1_coef, bP2_coef)
                                pBonds = combine2Dicts(pBonds, bP3_coef)

                                if compare2Dicts(rBonds, pBonds):
                                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3]],
                                                    products=[keys[p1], keys[p2], keys[p3]], coefs=coefs))

        # -----------------

        # 4 - 2

        print("Starting with 4 - 2")

        # for r2 in range(sz - 2):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for r3 in range(r2 + 1, sz - 1):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz):

                if keys[r4] == r1:
                    continue

                bondR4 = bonds[r4]

                for p1 in range(sz - 1):

                    if p1 == r2 or p1 == r3 or p1 == r4 or keys[p1] == r1:
                        continue

                    bondP1 = bonds[p1]

                    for p2 in range(p1 + 1, sz):

                        if p2 == r2 or p2 == r3 or p2 == r4 or keys[p2] == r1:
                            continue

                        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3], mols[r4]],
                                                         products=[mols[p1], mols[p2]],
                                                         atomList=[atomR1, atomList[r2], atomList[r3], atomList[r4],
                                                                   atomList[p1], atomList[p2]])

                        if coefList:

                            for coefs in coefList:

                                bondP2 = bonds[p2]

                                bR1_coef, bR2_coef, bR3_coef, bR4_coef, bP1_coef, bP2_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                                             multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondR4, coefs[3]),\
                                                                                             multiplyCoef(bondP1,coefs[4]), multiplyCoef(bondP2, coefs[5])

                                rBonds = combine2Dicts(bR1_coef, bR2_coef)
                                rBonds = combine2Dicts(rBonds, bR3_coef)
                                rBonds = combine2Dicts(rBonds, bR4_coef)
                                pBonds = combine2Dicts(bP1_coef, bP2_coef)

                                if compare2Dicts(rBonds, pBonds):
                                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3], keys[r4]],
                                                    products=[keys[p1], keys[p2]], coefs=coefs))
        #
        # # -----------------
        #
        # # 5 - 1
        #
        print("Starting with 5 - 1")

        # for r2 in range(sz - 3):
        #
        #     if keys[r2] == r1:
        #         continue
        #
        #     bondR2 = bonds[r2]

        for r3 in range(r2 + 1, sz - 2):

            if keys[r3] == r1:
                continue

            bondR3 = bonds[r3]

            for r4 in range(r3 + 1, sz - 1):

                if keys[r4] == r1:
                    continue

                bondR4 = bonds[r4]

                for r5 in range(r4 + 1, sz):

                    if keys[r5] == r1:
                        continue

                    bondR5 = bonds[r5]

                    for p1 in range(sz):

                        if p1 == r2 or p1 == r3 or p1 == r4 or p1 == r5 or keys[p1] == r1:
                            continue

                        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3], mols[r4], mols[r5]],
                                                         products=[mols[p1]],
                                                         atomList=[atomR1, atomList[r2], atomList[r3], atomList[r4],
                                                                   atomList[r5], atomList[p1]])

                        if coefList:

                            for coefs in coefList:

                                bondP1 = bonds[p1]

                                bR1_coef, bR2_coef, bR3_coef, bR4_coef, bR5_coef, bP1_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                                             multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondR4, coefs[3]),\
                                                                                             multiplyCoef(bondR5,coefs[4]), multiplyCoef(bondP1, coefs[5])

                                rBonds = combine2Dicts(bR1_coef, bR2_coef)
                                rBonds = combine2Dicts(rBonds, bR3_coef)
                                rBonds = combine2Dicts(rBonds, bR4_coef)
                                rBonds = combine2Dicts(rBonds, bR5_coef)

                                if compare2Dicts(rBonds, bP1_coef):
                                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3], keys[r4], keys[r5]],
                                                    products=[keys[p1]],coefs=coefs))


    return rxns

def processDiff(r1, m1, keys, config):

    noSpecies = config.noSpecies

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, noSpecies, atomR1, atomList, m1, mols, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_diff)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_diff_fast_3(r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1):

    if keys[r2] == r1:
        return None

    sz = len(keys)

    rxns = []

    bondR2 = bonds[r2]

    for p1 in range(sz):
        if p1 == r2 or keys[p1] == r1:
            continue

        coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1]],
                                         atomList=[atomR1, atomList[r2], atomList[p1]])

        if coefList:

            for coefs in coefList:

                bondP1 = bonds[p1]

                bR1_coef, bR2_coef, bP1_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(bondR2,coefs[1]), multiplyCoef(bondP1, coefs[2])
                rBonds = combine2Dicts(bR1_coef, bR2_coef)

                if compare2Dicts(rBonds, bP1_coef):
                    rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]], products=[keys[p1]], coefs=coefs))
    return rxns

def processDiff_fast_3(r1, m1, keys, config):

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_diff_fast_3)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    return rxnList

def parallel_diff_fast_4(r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1):

    if keys[r2] == r1:
        return None

    sz = len(keys)

    rxns = []

    bondR2 = bonds[r2]

    # 2 - 2
    for p1 in range(sz - 1):
        if p1 == r2 or keys[p1] == r1:
            continue

        bondP1 = bonds[p1]

        for p2 in range(p1 + 1, sz):

            if p2 == r2 or keys[p2] == r1:
                continue

            coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2]], products=[mols[p1], mols[p2]],
                                             atomList=[atomR1, atomList[r2], atomList[p1], atomList[p2]])

            if coefList:

                for coefs in coefList:

                    bondP2 = bonds[p2]

                    bR1_coef, bR2_coef, bP1_coef, bP2_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(bondR2,coefs[1]),\
                                                             multiplyCoef(bondP1, coefs[2]), multiplyCoef(bondP2,coefs[3])

                    rBonds = combine2Dicts(bR1_coef, bR2_coef)
                    pBonds = combine2Dicts(bP1_coef, bP2_coef)

                    if compare2Dicts(rBonds, pBonds):
                        rxns.append(printResultDiffParallel(reactants=[r1, keys[r2]], products=[keys[p1], keys[p2]], coefs=coefs))

    return rxns

def processDiff_fast_4(r1, m1, keys, config):

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_diff_fast_4)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_diff_fast_5(r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1):

    if keys[r2] == r1:
        return None

    sz = len(keys)

    rxns = []

    bondR2 = bonds[r2]

    for r3 in range(r2 + 1, sz):

        if keys[r3] == r1:
            continue

        bondR3 = bonds[r3]

        for p1 in range(sz - 1):

            if p1 == r2 or p1 == r3 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz):

                if p2 == r2 or p2 == r3 or keys[p2] == r1:
                    continue

                coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]], products=[mols[p1], mols[p2]],
                                                 atomList=[atomR1, atomList[r2], atomList[r3], atomList[p1],
                                                           atomList[p2]])

                if coefList:

                    for coefs in coefList:

                        bondP2 = bonds[p2]

                        bR1_coef, bR2_coef, bR3_coef, bP1_coef, bP2_coef = multiplyCoef(bondR1,coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                           multiplyCoef(bondR3,coefs[2]), multiplyCoef(bondP1, coefs[3]),\
                                                                           multiplyCoef(bondP2, coefs[4])

                        rBonds = combine2Dicts(bR1_coef, bR2_coef)
                        rBonds = combine2Dicts(rBonds, bR3_coef)
                        pBonds = combine2Dicts(bP1_coef, bP2_coef)

                        if compare2Dicts(rBonds, pBonds):
                            rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3]], products=[keys[p1], keys[p2]], coefs=coefs))
    return rxns

def processDiff_fast_5(r1, m1, keys, config):

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_diff_fast_5)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def parallel_diff_fast_6(r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1):

    if keys[r2] == r1:
        return None

    sz = len(keys)

    rxns = []

    bondR2 = bonds[r2]

    for r3 in range(r2 + 1, sz):

        if keys[r3] == r1:
            continue

        bondR3 = bonds[r3]

        for p1 in range(sz - 2):

            if p1 == r2 or p1 == r3 or keys[p1] == r1:
                continue

            bondP1 = bonds[p1]

            for p2 in range(p1 + 1, sz - 1):

                if p2 == r2 or p2 == r3 or keys[p2] == r1:
                    continue

                bondP2 = bonds[p2]

                for p3 in range(p2 + 1, sz):

                    if p3 == r2 or p3 == r3 or keys[p3] == r1:
                        continue

                    coefList = BalanceEq().balanceEq(reactants=[m1, mols[r2], mols[r3]],
                                                     products=[mols[p1], mols[p2], mols[p3]],
                                                     atomList=[atomR1, atomList[r2], atomList[r3], atomList[p1],
                                                               atomList[p2], atomList[p3]])

                    if coefList:

                        for coefs in coefList:

                            bondP3 = bonds[p3]

                            bR1_coef, bR2_coef, bR3_coef, bP1_coef, bP2_coef, bP3_coef = multiplyCoef(bondR1, coefs[0]), multiplyCoef(bondR2, coefs[1]),\
                                                                                         multiplyCoef(bondR3, coefs[2]), multiplyCoef(bondP1, coefs[3]),\
                                                                                         multiplyCoef(bondP2, coefs[4]), multiplyCoef(bondP3, coefs[5])

                            rBonds = combine2Dicts(bR1_coef, bR2_coef)
                            rBonds = combine2Dicts(rBonds, bR3_coef)
                            pBonds = combine2Dicts(bP1_coef, bP2_coef)
                            pBonds = combine2Dicts(pBonds, bP3_coef)

                            if compare2Dicts(rBonds, pBonds):
                                rxns.append(printResultDiffParallel(reactants=[r1, keys[r2], keys[r3]],
                                                products=[keys[p1], keys[p2], keys[p3]], coefs=coefs))

    return rxns

def processDiff_fast_6(r1, m1, keys, config):

    sz = len(keys)

    mols = [Chem.MolFromSmiles(keys[s]) for s in range(sz)]
    bonds = [findBonds(m, coef=1, config=config) for m in mols]

    atomR1 = defaultdict(int)

    for atom in m1.GetAtoms():
        atomR1[atom.GetSymbol()] += 1
        atomR1["H"] += atom.GetTotalNumHs()

    atomList = []

    for s in range(sz):
        m = mols[s]

        atoms = defaultdict(int)

        for atom in m.GetAtoms():
            atoms[atom.GetSymbol()] += 1
            atoms["H"] += atom.GetTotalNumHs()

        atomList.append(atoms)

    bondR1 = findBonds(m1, coef=1, config=config)

    args = [[r2, r1, keys, atomR1, atomList, m1, mols, bonds, bondR1] for r2 in range(sz) if keys[r2] != r1]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_diff_fast_6)

    rxnList = mp.imap_unordered(func, args, chunksize=1)

    mp.close()
    mp.join()

    return rxnList

def generate_rxns_complete(keys, s1, m1, fout, config):

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

            print("\nWarning !!! The reaction generation process may take long !!!\nPlease wait or choose another faster method such as stochastic or sameCoef=True !!!")
            rxnList = processDiff(s1, m1, keys, config)
        else:
            if noSpecies > 6:
                noSpecies = 6

            fastDiffProcessors = [processDiff_fast_3, processDiff_fast_4, processDiff_fast_5, processDiff_fast_6]

            index = noSpecies - 3

            rxnList = fastDiffProcessors[index](s1, m1, keys, config)

    resList = []

    for rxns in rxnList:
        if rxns and len(rxns) > 0:
            resList += rxns

    sCnt = len(resList)

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

                line = resList[index]
                fout.write("%*d\t\t%*s\n" % (maxCnt, index + 1, maxWidth, line))

            fout.write("\n=========================================================================\n")

    print("End generating reactions...")
    print("Total isodesmic reactions found: " + str(sCnt))
    print("--------------------------------\n")

    return sCnt, resList
