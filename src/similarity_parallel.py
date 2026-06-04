from utility_module import np, combine2Dicts, math, defaultdict, Chem, AllChem, multiplyCoef, Pairs, Torsions, Pool, partial, adjustNoProcessors, divide_work
from coulomb_matrix import compute_coulomb_matrix, coulomb_matrix_eig

def cosineSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    rxnSim = 0

    commonSet = set(rDict.keys()) & set(pDict.keys())

    for elem in commonSet:
        rxnSim += rDict[elem] * pDict[elem]

    return rxnSim * 1.0 / (1e-3 + np.linalg.norm(list(rDict.values())) * np.linalg.norm(list(pDict.values())))

def cosineSimilarity_2(rDict1, rDict2, pDict1, pDict2):

    rxnSim = 0

    commonSet = set(rDict1.keys()) & set(pDict1.keys())

    for elem in commonSet:
        rxnSim += rDict1[elem] * pDict1[elem]

    commonSet = set(rDict2.keys()) & set(pDict2.keys())

    for elem in commonSet:
        rxnSim += rDict2[elem] * pDict2[elem]

    rVals = list(rDict1.values()) + list(rDict2.values())
    pVals = list(pDict1.values()) + list(pDict2.values())

    return rxnSim * 1.0 / (1e-3 + np.linalg.norm(rVals) * np.linalg.norm(pVals))

def tanimotoSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    numerator = 0
    denominator = 0

    commonSet = set(rDict.keys()) & set(pDict.keys())

    for elem in commonSet:
        numerator += min(rDict[elem], pDict[elem])
        denominator += max(rDict[elem], pDict[elem])

    return numerator * 1.0 / (1e-3 + denominator)

def diceSimilarity(rDict, pDict):

    tanimoto = tanimotoSimilarity(rDict, pDict)

    return 2.0 * tanimoto / (1 + tanimoto)

def euclideanSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    unionSet = set(rDict.keys()) | set(pDict.keys())

    # rNorm = np.linalg.norm(list(rDict.values()))
    # pNorm = np.linalg.norm(list(pDict.values()))

    # if rNorm == 0 or pNorm == 0:
    #     return 0

    distance = 0

    for elem in unionSet:
        # rElem = rDict.get(elem, 0) / (rNorm * 1.0)
        # pElem = pDict.get(elem, 0) / (pNorm * 1.0)
        rElem = rDict.get(elem, 0)
        pElem = pDict.get(elem, 0)
        # if rElem == 0:
        #     print("R:", elem)
        # if pElem == 0:
        #     print("P:", elem)
        distance += (rElem - pElem) ** 2

    return 1.0 / (1e-3 + math.sqrt(distance))

def manhattanSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    unionSet = set(rDict.keys()) | set(pDict.keys())

    # rNorm = np.linalg.norm(list(rDict.values()))
    # pNorm = np.linalg.norm(list(pDict.values()))
    #
    # if rNorm == 0 or pNorm == 0:
    #     return 0

    distance = 0

    for elem in unionSet:
        # rElem = rDict.get(elem, 0) / (rNorm * 1.0)
        # pElem = pDict.get(elem, 0) / (pNorm * 1.0)
        rElem = rDict.get(elem, 0)
        pElem = pDict.get(elem, 0)
        distance += abs(rElem - pElem)

    return 1.0 / (1e-3 + distance)

def baryCurtisSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    unionSet = set(rDict.keys()) | set(pDict.keys())

    distance = 0

    # rNorm = np.linalg.norm(list(rDict.values()))
    # pNorm = np.linalg.norm(list(pDict.values()))
    #
    # if rNorm == 0 or pNorm == 0:
    #     return 0

    for elem in unionSet:
        # rElem = rDict.get(elem, 0) / (rNorm * 1.0)
        # pElem = pDict.get(elem, 0) / (pNorm * 1.0)
        rElem = rDict.get(elem, 0)
        pElem = pDict.get(elem, 0)
        distance += abs(rElem - pElem) * 1.0 / (abs(rElem) + abs(pElem))

    return 1.0 / (1e-3 + distance)

def canberraSimilarity(rDict, pDict):

    if len(rDict.keys()) == 0 and len(pDict.keys()) == 0:
        return 0

    unionSet = set(rDict.keys()) | set(pDict.keys())

    distance = 0

    # rNorm = np.linalg.norm(list(rDict.values()))
    # pNorm = np.linalg.norm(list(pDict.values()))

    for elem in unionSet:
        # rElem = rDict.get(elem, 0) / (rNorm * 1.0)
        # pElem = pDict.get(elem, 0) / (pNorm * 1.0)
        rElem = rDict.get(elem, 0)
        pElem = pDict.get(elem, 0)
        distance += abs(rElem - pElem) * 1.0 / abs(rElem + pElem)

    return 1.0 / (1e-3 + distance)

def describeAtom(atom):

    curFrag = ""

    symbol = atom.GetSymbol()
    hybridization = atom.GetHybridization()
    nRadicalHs = atom.GetNumRadicalElectrons()
    nHs = atom.GetTotalNumHs()

    curFrag += symbol

    # curFrag += "n" + str(len(atom.GetNeighbors()))

    # if symbol == "C":

    if nRadicalHs > 0:
        curFrag += "[H" + str(nHs) + "]"
    else:
        curFrag += "H" + str(nHs)

    if atom.GetIsAromatic():
        curFrag += "{A}"
    else:
        if atom.IsInRing():
            curFrag += "{R}"

    if symbol == "C":

        curFrag += "{" + str(hybridization) + "}"

        neighbors = atom.GetNeighbors()

        carbonCount = 0

        for neighbor in neighbors:
            if neighbor.GetSymbol() == "C":
                carbonCount += 1

        curFrag += "{d" + str(carbonCount) + "}"
    else:
        curFrag += "{d" + str(atom.GetDegree()) + "}"


    return curFrag

def circularfrag(species, radius=2, useChirality=True):

    # print(species)

    bitInfo = {}
    fragmol = defaultdict(int)
    fragList = defaultdict(list)
    mol = Chem.MolFromSmiles(species)
    AllChem.GetMorganFingerprint(mol, radius=radius, useFeatures=False, bitInfo=bitInfo, useChirality=useChirality)

    for bit, info in bitInfo.items():

        for atmidx, rad in info:

            # print "Cur index:", atmidx

            env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, atmidx)

            submol = Chem.PathToSubmol(mol, env)
            smi = Chem.MolToSmiles(submol, isomericSmiles=True)
            # smi = Chem.MolToSmiles(submol)

            if smi != "":

                # print smi
                # print sorted(env), len(mol.GetBonds())

                curFrags = ""

                found = False

                if smi in fragList:
                    for curMol in fragList[smi]:
                        if curMol[0].HasSubstructMatch(submol):
                            found = True
                            curFrags = curMol[1]
                            break

                # print found

                if not found:

                    if smi in species:
                        sorted_bonds = sorted(env)
                    else:
                        sorted_bonds = sorted(env, reverse=True)
                    # print sorted_bonds

                    atomList = []

                    startAtom = mol.GetAtomWithIdx(mol.GetBondWithIdx(sorted_bonds[0]).GetBeginAtomIdx())
                    atomList.append(describeAtom(startAtom))

                    for bIdx in sorted_bonds:

                        bond = mol.GetBondWithIdx(bIdx)
                        atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                        atomList.append(describeAtom(atom))

                    # print(atomList)

                    index, fragIndex = 0, 0
                    curFrags = smi
                    # curFrags = curFrags.replace("H", "")

                    # print(curFrags)

                    while index < len(curFrags):
                        c = curFrags[index].lower()
                        if c != "(" and c != ")" and c != "[" and c != "]" and c != "/" and c != "\\"\
                                    and c != "=" and c != "#" and c != "+" and c != "-"\
                                    and c != "@" and c != "h" and (c <= "0" or c >= "9"):

                            # print(index, len(curFrags), c)

                            if index < len(curFrags) - 1:
                                curFrags = curFrags[:index] + atomList[fragIndex] + curFrags[index + 1:]
                            else:
                                curFrags = curFrags[:index] + atomList[fragIndex]
                            index += len(atomList[fragIndex])
                            # print curFrags, index
                            fragIndex += 1
                        else:
                            index += 1

                    fragList[smi].append((submol, curFrags))

                # print curFrags

                # non_H_atom_cnt = 0
                # for c in smi:
                #     c = c.lower()
                #     if c == "c" or c == "o" or c == "n" or c == "s":
                #         non_H_atom_cnt += 1

                # fragmol[curFrags] += 1.0 / (non_H_atom_cnt * 1.0)
                fragmol[curFrags] += 1

            else:
                atom = mol.GetAtomWithIdx(atmidx)

                # print curFrags
                fragmol[describeAtom(atom)] += 1

    # print(dict(fragmol))

    return fragmol

def circularfrag_v2(species, radius=2, useChirality=True):

    bitInfo = {}
    fragmol = defaultdict(int)
    fragList = defaultdict(list)
    mol = Chem.MolFromSmiles(species)

    atoms = mol.GetAtoms()
    for atom in atoms:
        fragmol[describeAtom(atom)] += 1

    bondList = {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": "+"}

    bonds = mol.GetBonds()

    for bond in bonds:
        startAtom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        endAtom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        bondType = bondList[str(bond.GetBondType())]

        sStr = describeAtom(startAtom)
        eStr = describeAtom(endAtom)

        if sStr < eStr:
            sStr, eStr = eStr, sStr

        conn = sStr + bondType + eStr
        fragmol[conn] += 1

    AllChem.GetMorganFingerprint(mol, radius=radius, useFeatures=False, bitInfo=bitInfo, useChirality=useChirality)

    for bit, info in bitInfo.items():

        for atmidx, rad in info:

            # print "Cur index:", atmidx

            env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, atmidx)

            submol = Chem.PathToSubmol(mol, env)
            smi = Chem.MolToSmiles(submol, isomericSmiles=True)
            # smi = Chem.MolToSmiles(submol)

            if smi != "":

                CCnt = smi.lower().count('c')
                OCnt = smi.lower().count('o')

                if CCnt + OCnt < 3:
                    continue

                # print smi
                # print sorted(env), len(mol.GetBonds())

                curFrags = ""

                found = False

                if smi in fragList:
                    for curMol in fragList[smi]:
                        if curMol[0].HasSubstructMatch(submol):
                            found = True
                            curFrags = curMol[1]
                            break

                # print found

                if not found:

                    if smi in species:
                        sorted_bonds = sorted(env)
                    else:
                        sorted_bonds = sorted(env, reverse=True)
                    # print sorted_bonds

                    atomList = []

                    startAtom = mol.GetAtomWithIdx(mol.GetBondWithIdx(sorted_bonds[0]).GetBeginAtomIdx())
                    atomList.append(describeAtom(startAtom))

                    for bIdx in sorted_bonds:

                        bond = mol.GetBondWithIdx(bIdx)
                        atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                        atomList.append(describeAtom(atom))

                    # print(atomList)

                    index, fragIndex = 0, 0
                    curFrags = smi
                    # curFrags = curFrags.replace("H", "")

                    # print(curFrags)

                    while index < len(curFrags):
                        c = curFrags[index].lower()
                        if c != "(" and c != ")" and c != "[" and c != "]" and c != "/" and c != "\\"\
                                    and c != "=" and c != "#" and c != "+" and c != "-"\
                                    and c != "@" and c != "h" and (c <= "0" or c >= "9"):

                            # print(index, len(curFrags), c)

                            if index < len(curFrags) - 1:
                                curFrags = curFrags[:index] + atomList[fragIndex] + curFrags[index + 1:]
                            else:
                                curFrags = curFrags[:index] + atomList[fragIndex]
                            index += len(atomList[fragIndex])
                            # print curFrags, index
                            fragIndex += 1
                        else:
                            index += 1

                    fragList[smi].append((submol, curFrags))

                # print curFrags

                # non_H_atom_cnt = 0
                # for c in smi:
                #     c = c.lower()
                #     if c == "c" or c == "o" or c == "n" or c == "s":
                #         non_H_atom_cnt += 1

                # fragmol[curFrags] += 1.0 / (non_H_atom_cnt * 1.0)
                fragmol[curFrags] += 1

    # print(dict(fragmol))

    return fragmol

def filter_frags(species, radius=2, useChirality=True):

    bitInfo = {}
    fragmol = defaultdict(int)
    mol = Chem.MolFromSmiles(species)
    AllChem.GetMorganFingerprint(mol, radius=radius, useFeatures=True, bitInfo=bitInfo, useChirality=useChirality)

    for bit, info in bitInfo.items():

        for atmidx, rad in info:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, atmidx)
            submol = Chem.PathToSubmol(mol, env)
            smi = Chem.MolToSmiles(submol)

            if smi != '':
                fragmol[smi] += 1
            else:
                fragmol[mol.GetAtomWithIdx(atmidx).GetSymbol()] += 1

    return fragmol

def bag_of_bonds(species):

    bondList = {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": "+"}

    fragmol = defaultdict(int)
    m = Chem.MolFromSmiles(species)
    bonds = m.GetBonds()

    for bond in bonds:
        startAtom = m.GetAtomWithIdx(bond.GetBeginAtomIdx())
        endAtom = m.GetAtomWithIdx(bond.GetEndAtomIdx())
        bondType = bondList[str(bond.GetBondType())]

        sStr = startAtom.GetSymbol()
        eStr = endAtom.GetSymbol()

        if sStr < eStr:
            sStr, eStr = eStr, sStr

        conn = sStr + bondType + eStr

        fragmol[conn] += 1

    # print(dict(fragmol))

    return fragmol

def hr_bag_of_bonds(species):

    bondList = {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": "+"}

    fragmol = defaultdict(int)
    m = Chem.MolFromSmiles(species)
    bonds = m.GetBonds()

    for bond in bonds:
        startAtom = m.GetAtomWithIdx(bond.GetBeginAtomIdx())
        endAtom = m.GetAtomWithIdx(bond.GetEndAtomIdx())
        bondType = bondList[str(bond.GetBondType())]

        # sStr = startAtom.GetSymbol()
        # eStr = endAtom.GetSymbol()

        sStr = describeAtom(startAtom)
        eStr = describeAtom(endAtom)

        if sStr < eStr:
            sStr, eStr = eStr, sStr

        conn = sStr + bondType + eStr

        fragmol[conn] += 1

    # print(dict(fragmol))

    return fragmol

def coulomb_matrices_sim(reactants, products):
    rMatrices = [coulomb_matrix_eig(Chem.MolFromSmiles(r)) for r in reactants]
    pMatrices = [coulomb_matrix_eig(Chem.MolFromSmiles(p)) for p in products]

    rSize = len(rMatrices)
    pSize = len(pMatrices)

    # if len(rMatrices) > 1:
    #     rSim = 0
    #     for i in range(rSize - 1):
    #         for j in range(i + 1, rSize):
    #             rSim += np.linalg.norm(rMatrices[j] - rMatrices[i]) / (np.linalg.norm(rMatrices[j]) * np.linalg.norm(rMatrices[i]))
    #     rSim /= (rSize * 1.0)
    # else:
    #     rSim = 1
    #
    # if len(pMatrices) > 1:
    #     pSim = 0
    #     for i in range(pSize - 1):
    #         for j in range(i + 1, pSize):
    #             pSim += np.linalg.norm(pMatrices[j] - pMatrices[i])
    #     pSim /= (pSize * 1.0)
    # else:
    #     pSim = 1

    rSim = 1
    pSim = 1

    norms = 0
    for rMatrix in rMatrices:
        for pMatrix in pMatrices:
            norms += np.linalg.norm(pMatrix - rMatrix) * rSim * pSim / (np.linalg.norm(pMatrix) * np.linalg.norm(rMatrix))

    return norms / (rSize * pSize * 1.0)

def similarity_dict(dList):

    sz = len(dList)

    sim = 0

    for i in range(sz):
        for j in range(i + 1, sz):
            sim += cosineSimilarity(dList[j], dList[i])

    return sim / (sz * 1.0)

def parallel_similarity(line, featureMethod, radius, useFeatures, useChirality):

    toks = line.split("--->")
    # toks = [t.replace(" ", "") for t in toks]

    lSide, rSide = toks[0].split(" + "), toks[1].split(" + ")

    coefs = []

    reactants = []
    products = []

    for l in lSide:
        toks = l.split()
        if len(toks) == 2:
            coefs.append(int(toks[0]))
            reactants.append(toks[1])
        else:
            coefs.append(1)
            reactants.append(toks[0])

    for r in rSide:
        toks = r.split()
        if len(toks) == 2:
            coefs.append(int(toks[0]))
            products.append(toks[1])
        else:
            coefs.append(1)
            products.append(toks[0])

    if featureMethod == "cm":
        rxnSim = coulomb_matrices_sim(reactants, products)
    else:
        if featureMethod == "hrcf":
            rDicts = [circularfrag(r, radius=radius, useChirality=useChirality) for r in reactants]
            pDicts = [circularfrag(p, radius=radius, useChirality=useChirality) for p in products]
            # rDicts = [circularfrag_v2(r, radius=radius, useChirality=useChirality) for r in reactants]
            # pDicts = [circularfrag_v2(p, radius=radius, useChirality=useChirality) for p in products]
        elif featureMethod == "cf":
            rDicts = [filter_frags(r, radius=radius, useChirality=useChirality) for r in reactants]
            pDicts = [filter_frags(p, radius=radius, useChirality=useChirality) for p in products]
        elif featureMethod == "bob":
            rDicts = [bag_of_bonds(r) for r in reactants]
            pDicts = [bag_of_bonds(p) for p in products]
        elif featureMethod == "hr_bob":
            rDicts = [hr_bag_of_bonds(r) for r in reactants]
            pDicts = [hr_bag_of_bonds(p) for p in products]
        else:

            # rDicts = [dict(DataStructs.UIntSparseIntVect.GetNonzeroElements(
            #     AllChem.GetMorganFingerprint(Chem.MolFromSmiles(r), radius, useFeatures=useFeatures,
            #                                  useChirality=useChirality))) for r in reactants]
            #
            # pDicts = [dict(DataStructs.UIntSparseIntVect.GetNonzeroElements(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(p), radius, useFeatures=useFeatures,
            #                                  useChirality=useChirality))) for p in products]
            if featureMethod == "fcfp" or featureMethod == "ecfp":
                rDicts = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(r), radius, useFeatures=useFeatures,
                                                       useChirality=useChirality).GetNonzeroElements() for r in
                          reactants]

                pDicts = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(p), radius, useFeatures=useFeatures,
                                                       useChirality=useChirality).GetNonzeroElements() for p in
                          products]
            elif featureMethod == "ap":
                rDicts = [Pairs.GetAtomPairFingerprint(Chem.MolFromSmiles(r),
                                                       includeChirality=useChirality).GetNonzeroElements()
                          for r in reactants]

                pDicts = [Pairs.GetAtomPairFingerprint(Chem.MolFromSmiles(p),
                                                       includeChirality=useChirality).GetNonzeroElements()
                          for p in products]
            elif featureMethod == "tt":
                rDicts = [Torsions.GetTopologicalTorsionFingerprint(Chem.MolFromSmiles(r),
                                                                    includeChirality=useChirality).GetNonzeroElements()
                          for r in reactants]

                pDicts = [Torsions.GetTopologicalTorsionFingerprint(Chem.MolFromSmiles(p),
                                                                    includeChirality=useChirality).GetNonzeroElements()
                          for p in products]

        # if len(rDicts) > 1:
        #     rSim = similarity_dict(rDicts)
        # else:
        #     rSim =  1
        #
        # if len(pDicts) > 1:
        #     pSim = similarity_dict(pDicts)
        # else:
        #     pSim = 1

        # print(rSim, pSim)

        for i in range(len(rDicts)):
            rDicts[i] = multiplyCoef(rDicts[i], coefs[i])

        for i in range(len(pDicts)):
            pDicts[i] = multiplyCoef(pDicts[i], coefs[i + len(reactants)])

        rDict = {}
        for d in rDicts:
            rDict = combine2Dicts(rDict, d)

        pDict = {}
        for d in pDicts:
            pDict = combine2Dicts(pDict, d)

        # print(rDict, pDict)

        # rxnSim = cosineSimilarity(rDict, pDict) * rSim * pSim

        # print(rxnSim)
        # rxnSim = cosineSimilarity_2(rDicts[0], rDicts[1], pDicts[0], pDicts[1])
        rxnSim = cosineSimilarity(rDict, pDict)
        # rxnSim = tanimotoSimilarity(rDict, pDict)
        # rxnSim = diceSimilarity(rDict, pDict)
        # rxnSim = euclideanSimilarity(rDict, pDict)
        # rxnSim = manhattanSimilarity(rDict, pDict)
        # rxnSim = baryCurtisSimilarity(rDict, pDict)
        # rxnSim = canberraSimilarity(rDict, pDict)

    return line, rxnSim

def checkSimilarity(rxnList, config, fout):

    # Similarity checking parameters
    radius, useChirality = config.radius, config.chirality

    featureMethod = config.featureMethod.lower()

    if featureMethod == "fcfp":
        useFeatures = True
    else:
        useFeatures = False

    reactions = {}
    # maxWidth = 15 * config.noSpecies
    # maxCnt = len(str(rxnCnt))

    print("Start checking similarity...")

    args = [[rxn, featureMethod, radius, useFeatures, useChirality] for rxn in rxnList]

    no_processors = adjustNoProcessors(config.noProcessors)

    mp = Pool(processes=no_processors)

    func = partial(divide_work, parallel_similarity)

    if no_processors > len(rxnList):
        no_processors = len(rxnList)

    chunkSz = max(1, len(rxnList) // no_processors)

    # print(chunkSz)

    # chunkSz = 1
    simList = mp.imap_unordered(func, args, chunksize=chunkSz)

    mp.close()
    mp.join()

    # with open(simFile, "w") as fout:
    #     fout.write("%*s\t\t%*s\t\t%s\n" % (maxCnt, "No.", maxWidth, "Reaction", "Similarity"))

    sCnt = len(rxnList)
    maxWidth = 15 * config.noSpecies
    maxCnt = len(str(sCnt))

    if not config.analysisOn:

        fout.write("SUMMARY\n\n")
        fout.write("Total number of reactions found: " + str(sCnt) + "\n\n")
        fout.write("=========================================================================\n")

        if sCnt > 0:
            fout.write("LIST OF REACTIONS FOUND\n\n")
            fout.write("%*s\t\t%*s\t\t%s\n\n" % (maxCnt, "No.", maxWidth, "Reaction", "Similarity"))

    maxSim = -1
    maxRxn = ""
    cnt = 1

    for sim in simList:

        line = sim[0]
        rxnSim = sim[1]

        # fout.write("%*d\t\t%*s\t\t%.5f\n" % (maxCnt, cnt, maxWidth, line, rxnSim))

        if not config.analysisOn:
            fout.write("%*d\t\t%*s\t\t%.5f\n" % (maxCnt, cnt, maxWidth, line, rxnSim))
            cnt += 1

        if rxnSim > maxSim:
            maxRxn = line
            maxSim = rxnSim

        reactions[line] = rxnSim

    if not config.analysisOn:
        fout.write("\n=========================================================================\n")

    print("End checking similarity...")

    print("Reaction with max similarity: " + maxRxn)
    print("Max Similarity: " + str(maxSim))

    print("--------------------------------\n")

    return reactions
