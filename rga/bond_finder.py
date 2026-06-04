from .utils import defaultdict, Chem, NumValenceElectrons

bondList = {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": "+"}

def atomDegree(atom, symbol):
    neighbors = atom.GetNeighbors()
    carbon_cnt = 0

    for neighbor in neighbors:
        if neighbor.GetSymbol() == symbol:
            carbon_cnt += 1

    return carbon_cnt

def findBonds(m, coef, config):

    lElectron, ePair, ePairType, atomWithHydro = config.lElectron, config.ePair, config.ePairType, config.atomWithHydro
    atomHybridization, atomHybridWithHydro, atomWithDegree, radical = config.atomHybridization, config.atomHybridWithHydro, config.atomDegree, config.radical
    hydroBond, normalBond, bondHybrid, bondHydro = config.hydroBond, config.normalBond, config.bondHybrid, config.bondHydro
    bondDegree, bondRadical, ringIncluded = config.bondDegree, config.bondRadical, config.ring

    atomTypes = config.atomType.split(",")
    bondAtomTypes = config.bondAtomType.split(",")

    global bondList

    bonds = defaultdict(int)

    if ePair or lElectron or ePairType:

        numValenceElectrons = NumValenceElectrons(m)

        if lElectron:
            bonds["LE"] = numValenceElectrons % 2

        if ePair:
            bonds['EP'] = int(numValenceElectrons / 2)

        if ePairType:
            bonds["EP_shared"] = m.GetNumBonds(onlyHeavy=False)
            bonds["EP_unshared"] = int(numValenceElectrons / 2)

    if ePair or hydroBond or atomWithHydro or atomHybridization or radical:
        for atom in m.GetAtoms():

            symbol = atom.GetSymbol()
            bond = "a" + symbol

            if hydroBond:
                for bondAtomType in bondAtomTypes:
                    if symbol.lower() == bondAtomType.lower():
                        nHs = atom.GetTotalNumHs()

                        bond += "-H"
                        bonds[bond] += nHs * coef

            for atomType in atomTypes:
                if symbol.lower() == atomType.lower():
                    if atomWithHydro:

                        bond = symbol

                        nHs = atom.GetTotalNumHs()
                        bond += "-" + str(nHs) + "H"

                        bonds[bond] += coef

                    if atomHybridization or atomHybridWithHydro:

                        bond = symbol

                        hybrid = str(atom.GetHybridization())
                        bond += "(" + hybrid + ")"

                        if atomHybridWithHydro:
                            nHs = atom.GetTotalNumHs()
                            bond += "-" + str(nHs) + "H"

                        bonds[bond] += coef

                    if atomWithDegree:

                        bond = symbol

                        degree = str(atomDegree(atom, "C"))
                        bond += "(d" + degree + ")"

                        bonds[bond] += coef

            if radical:
                radicalElectron = atom.GetNumRadicalElectrons()
                if radicalElectron > 0:
                    nExplicitHs = atom.GetNumExplicitHs()
                    rad_bond = "[" + symbol

                    if nExplicitHs > 0:
                        rad_bond += "H"

                        if nExplicitHs > 1:
                            rad_bond += str(nExplicitHs)

                    rad_bond += "]"
                    bonds[rad_bond] += coef

    if normalBond or bondHybrid or bondHydro or bondDegree or bondRadical or ePairType:

        for bond in m.GetBonds():

            startAtom = m.GetAtomWithIdx(bond.GetBeginAtomIdx())
            endAtom = m.GetAtomWithIdx(bond.GetEndAtomIdx())
            bondType = bondList[str(bond.GetBondType())]

            if ePairType:
                if bondType == "=" or bondType == "+":
                    bonds['EP_shared'] += 1
                elif bondType == "#":
                    bonds['EP_shared'] += 2

            if normalBond:
                sStr = startAtom.GetSymbol()
                eStr = endAtom.GetSymbol()

                if sStr < eStr:
                    sStr, eStr = eStr, sStr

                conn = sStr + bondType + eStr

                bonds[conn] += coef

            if bondHybrid or bondHydro or bondDegree or bondRadical:

                sStr = startAtom.GetSymbol()
                eStr = endAtom.GetSymbol()

                if bondRadical:

                    sRadElectron = startAtom.GetNumRadicalElectrons()
                    eRadElectron = endAtom.GetNumRadicalElectrons()

                    if sRadElectron > 0:
                        sExplicitHs = startAtom.GetNumExplicitHs()

                        if sExplicitHs > 1:
                            sStr = "[" + sStr + "H" + str(sExplicitHs) + "]"
                        elif sExplicitHs == 1:
                            sStr = "[" + sStr + "H]"
                        else:
                            sStr = "[" + sStr + "]"

                    if eRadElectron > 0:
                        eExplicitHs = endAtom.GetNumExplicitHs()

                        if eExplicitHs > 1:
                            eStr = "[" + eStr + "H" + str(eExplicitHs) + "]"
                        elif eExplicitHs == 1:
                            eStr = "[" + eStr + "H]"
                        else:
                            eStr = "[" + eStr + "]"

                for bondAtomType in bondAtomTypes:

                    if sStr.lower() == bondAtomType.lower():

                        if bondHybrid:
                            startHybrid = str(startAtom.GetHybridization())
                            sStr += "(" + startHybrid + ")"

                        if bondHydro:
                            startHydro = str(startAtom.GetTotalNumHs())
                            sStr += "(" + startHydro + "H)"

                        if bondDegree:
                            startDegree = str(atomDegree(startAtom, "C"))
                            sStr += "(d" + startDegree + ")"

                    if eStr.lower() == bondAtomType.lower():
                        if bondHybrid:
                            endHybrid = str(endAtom.GetHybridization())
                            eStr += "(" + endHybrid + ")"

                        if bondHydro:
                            endHydro = str(endAtom.GetTotalNumHs())
                            eStr += "(" + endHydro + "H)"

                        if bondDegree:
                            endDegree = str(atomDegree(endAtom, "C"))
                            eStr += "(d" + endDegree + ")"

                if sStr < eStr:
                    sStr, eStr = eStr, sStr

                conn = sStr + bondType + eStr

                bonds[conn] += coef

    if ePairType:
        bonds["EP_unshared"] -= bonds["EP_shared"]

    if ringIncluded:

        ssrs = Chem.GetSymmSSSR(m)

        for ssr in ssrs:

            ssr_sorted = sorted(ssr)

            if ssr != ssr_sorted:

                rems = ssr_sorted[1:]
                res = [ssr_sorted[0]]

                i = 0
                while len(rems) > 0:
                    curAtom = res[i]

                    j = 0
                    while j < len(rems):
                        remAtom = rems[j]
                        if m.GetBondBetweenAtoms(curAtom, remAtom):
                            rems.remove(remAtom)
                            res.append(remAtom)
                            i += 1
                            break
                        else:
                            j += 1

                ssr = res[:]
            else:
                ssr = ssr_sorted[:]

            bond = m.GetAtomWithIdx(ssr[0]).GetSymbol()

            for i in range(len(ssr) - 1):
                curIndex = ssr[i]
                nextIndex = ssr[i + 1]
                bond += bondList[str(m.GetBondBetweenAtoms(curIndex, nextIndex).GetBondType())] + m.GetAtomWithIdx(nextIndex).GetSymbol()

            bonds[bond] += coef

    return bonds
