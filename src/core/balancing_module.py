from utility_module import Chem, defaultdict, np, permutations, combine2Dicts, math, matrix_rank, re, sympy, Matrix, LpProblem, LpMinimize, LpVariable, LpInteger, LpStatus, PulpSolverError, GLPK_CMD

class BalanceEq:

    def __init__(self):
        self.found = False
        self.results = {}

    def processSpecies(self, rList, pList, atomList):

        rSize = len(rList)
        pSize = len(pList)

        atoms = atomList[0]

        totalLen = rSize + pSize

        for i in range(1, totalLen):
            atoms = combine2Dicts(atoms, atomList[i])

        resultArr = []

        lineCnt = 0
        for atom in atoms:

            curArr = []
            for index in range(rSize):
                curArr.append(atomList[index].get(atom, 0))

            for index in range(pSize):
                curArr.append(-atomList[index + rSize].get(atom, 0))

            resultArr.append(curArr)
            lineCnt += 1

        return resultArr

        # atomList["H"] = [0] * (rSize + pSize)
        #
        # for index in range(rSize):
        #
        #     m = rList[index]
        #
        #     nHs = 0
        #
        #     atoms = m.GetAtoms()
        #
        #     for atom in atoms:
        #
        #         nHs += atom.GetTotalNumHs()
        #
        #         symbol = atom.GetSymbol()
        #
        #         if not symbol in atomList:
        #             atomList[symbol] = [0] * (rSize + pSize)
        #
        #         atomList[symbol][index] += 1
        #
        #     if nHs > 0:
        #         atomList["H"][index] += nHs
        #
        # for index in range(pSize):
        #
        #     m = pList[index]
        #
        #     nHs = 0
        #
        #     atoms = m.GetAtoms()
        #
        #     for atom in atoms:
        #
        #         nHs += atom.GetTotalNumHs()
        #
        #         symbol = atom.GetSymbol()
        #
        #         if not symbol in atomList:
        #             atomList[symbol] = [0] * (rSize + pSize)
        #
        #         atomList[symbol][index + rSize] -= 1
        #
        #     if nHs > 0:
        #         atomList["H"][index + rSize] -= nHs
        #
        # return atomList

    def checkInteger(self, val):
        return math.fabs(val - int(val)) < 1e-9 or math.fabs(val - math.ceil(val)) < 1e-9

    def checkCoefs(self, coefs):

        for elem in coefs:
            if elem <= 0:
                return False

        return True

    def hasSameElems(self, coefs):

        startElem = coefs[0]
        sz = len(coefs)

        for i in range(1, sz):
            if coefs[i] != startElem:
                return False

        return True

    def findOptimalValue(self, curIndex, keys, chars, sols, maxCnt, step, valueFound):

        # global results
        # global found

        curKey = keys[curIndex]

        for i in range(1, maxCnt + 1):

            if not valueFound and i > maxCnt - step:
                valueFound = True

            oldVal = chars[curKey]
            chars[curKey] = i

            if valueFound:

                # print curKey, i, chars, results, minVal

                isInteger = True
                smallerThanOne = True

                GCD = 1

                for elem in sols:

                    for c in chars.keys():
                        elem = elem.replace(c, str(chars[c]))

                    # Avoid numerical error
                    val = round(eval(elem), 9)

                    if val < 1.0:
                        smallerThanOne = False
                        break

                    if not self.checkInteger(val):
                        isInteger = False
                        break

                    GCD = self.gcd(GCD, int(val))

                if smallerThanOne and isInteger:

                    if GCD == 1:
                        self.found = True
                        self.results = chars.copy()
                        return

            if curIndex + 1 < len(chars.keys()):
                self.findOptimalValue(curIndex + 1, keys, chars, sols, maxCnt, step, valueFound)

            if self.found:
                return

            chars[curKey] = oldVal

    def findSol(self, M, totalLen):

        sols = ["" for _ in range(totalLen)]
        noFreeIndices = [index for index in range(totalLen)]
        chars = {}
        cIndex = -1

        for r in M:
            firstOne = False
            for index in range(len(r)):

                val = float(r[index])

                if not firstOne and val == 1:
                    firstOne = True
                    cIndex = index
                else:

                    if val == 0:
                        continue

                    val = -val

                    char = chr(ord('a') + index)

                    if sols[cIndex] == "":
                        sols[cIndex] = str(val) + "*" + char
                    else:
                        sols[cIndex] += " + " + str(val) + "*" + char
                    if char not in chars:
                        chars[char] = 1

        for index in range(totalLen):
            if sols[index] == "":
                noFreeIndices.remove(index)
                sols[index] = chr(ord('a') + index)

        selectedSols = [sols[val] for val in noFreeIndices]

        return {'sols': sols, 'symbols': chars, "selectedSols": selectedSols}

    def gcd(self, a, b):
        """Compute the greatest common divisor of a and b"""
        while b > 0:
            a, b = b, a % b
        return a

    def lcm(self, a, b):
        """Compute the lowest common multiple of a and b"""
        return a * b // self.gcd(a, b)

    def balanceEq(self, reactants, products, atomList):

        # global results
        # global found

        rLen = len(reactants)
        pLen = len(products)
        totalLen = rLen + pLen

        # atomList = []
        #
        # for s in range(rLen):
        #     m = reactants[s]
        #
        #     atoms = defaultdict(int)
        #
        #     for atom in m.GetAtoms():
        #         atoms[atom.GetSymbol()] += 1
        #         atoms["H"] += atom.GetTotalNumHs()
        #
        #     atomList.append(atoms)
        #
        # for s in range(pLen):
        #     m = products[s]
        #
        #     atoms = defaultdict(int)
        #
        #     for atom in m.GetAtoms():
        #         atoms[atom.GetSymbol()] += 1
        #         atoms["H"] += atom.GetTotalNumHs()
        #
        #     atomList.append(atoms)
        #
        # atomList = self.processSpecies(rList=reactants, pList=products, atomList={})
        #
        # A = list(atomList.values())[:]

        A = self.processSpecies(rList=reactants, pList=products, atomList=atomList)
        # print("Execution time:", (time.clock() - start), "s.")

        # Checking whether both sides contains same atom types
        for row in A:
            rCnt = 0
            for col_index in range(rLen):
                if row[col_index] != 0:
                    rCnt += 1
            if rCnt == 0:
                return None

            pCnt = 0
            for col_index in range(pLen):
                if row[col_index + rLen] != 0:
                    pCnt += 1
            if pCnt == 0:
                return None

        # print(atomList)

        mRank = matrix_rank(A)

        # print(mRank)

        # There're many solutions for the coefficients of the equation (multiples of the least-value set)
        # Therefore, the rank of the matrix must be at most N - 1 (N is the total no. of species of both sides of the equations)
        if mRank >= totalLen:
            return None

        if mRank == totalLen - 1:

            # find first basis vector == primary solution

            coefs = sympy.Matrix(A).nullspace()[0]

            if self.checkCoefs(coefs):
                if not self.hasSameElems(coefs):

                    p = re.compile("\/(\d+)")

                    LCM = 1

                    for elem in coefs:
                        tok = re.findall(pattern=p, string=str(elem))

                        if tok:
                            LCM = self.lcm(LCM, int(tok[0]))

                    return [[int(elem * LCM) for elem in coefs]]

                return [[elem for elem in coefs]]

            return None

        else:

            checkedA = np.asarray(Matrix(A).rref()[0])

            # print(checkedA)

            for rIndex in range(mRank):
                row = checkedA[rIndex]
                negFound = False
                for col in row:
                    if col < 0:
                        negFound = True
                        break

                if not negFound:
                    return None

            mSols = self.findSol(checkedA, totalLen)

            # print(mSols)

            sols = mSols['sols']
            chars = mSols['symbols']
            selectedSols = mSols['selectedSols']

            # print(sols)
            # print(chars)

            if totalLen > 4:
                keyList = [list(chars.keys())]
            else:
                keyList = permutations(list(chars.keys()))

            coefs = []

            if totalLen > 4:
                step = 1
                cutOffMaxCnt = 1
            else:
                step = 10
                cutOffMaxCnt = 30
                # cutOffMaxCnt = 100

            # step = 10
            # cutOffMaxCnt = 30

            # step = 5
            # cutOffMaxCnt = 10

            # step = 5
            # cutOffMaxCnt = 20

            # print(step, cutOffMaxCnt)

            for keys in keyList:

                self.results = {}
                self.found = False
                maxCnt = step

                # Trying ascending way
                while True:
                    self.findOptimalValue(curIndex=0, keys=keys, chars=chars.copy(), sols=selectedSols, maxCnt=maxCnt, step=step,
                                          valueFound=False)

                    if self.results:
                        break

                    maxCnt += step
                    # print(maxCnt)
                    if maxCnt > cutOffMaxCnt:

                        # return None

                        # print(len(reactants), len(products))
                        # print(Chem.MolToSmiles(reactants[0]), Chem.MolToSmiles(reactants[1]),
                        #       Chem.MolToSmiles(reactants[2]))
                        # print(Chem.MolToSmiles(products[0]))
                        # print(checkedA)

                        varList = []
                        objective = None
                        # print(checkedA)

                        for index in range(totalLen):
                            varList.append(LpVariable(chr(ord('t') + index), 1, 50, LpInteger))
                            # varList.append(LpVariable("x" + str(index), 1, None, LpInteger))

                            if index == 0:
                                objective = varList[-1]
                            else:
                                objective += varList[-1]

                        prob = LpProblem("Balancingrxn", LpMinimize)

                        prob += objective

                        for rIndex in range(mRank):
                            row = checkedA[rIndex]
                            curCol = None
                            firstFound = False

                            for cIndex in range(len(row)):

                                col = row[cIndex]

                                if col == 0:
                                    continue

                                if not firstFound:
                                    firstFound = True
                                    curCol = col * 1.0 * varList[cIndex]
                                else:
                                    curCol += col * 1.0 * varList[cIndex]

                            prob += curCol == 0

                        # print(prob)

                        try:
                            prob.solve(solver=GLPK_CMD(msg=0))
                            # prob.solve()

                            if LpStatus[prob.status] == 'Optimal':
                                return [[int(round(v.varValue, 9)) for v in prob.variables()]]

                            return None

                        except PulpSolverError:
                            # print(A)
                            # print("This equation is infeasible")
                            return None
                        # if maxCnt > 100:
                        #     print(Chem.MolToSmiles(reactants[0]), Chem.MolToSmiles(reactants[1]), Chem.MolToSmiles(reactants[2]))
                        #     print(Chem.MolToSmiles(products[0]))
                        #     print(maxCnt)

                sols_copy = ",".join(sols)

                for c in self.results:
                    sols_copy = sols_copy.replace(c, str(self.results[c]))

                sols_copy = sols_copy.split(",")
                # print(sols_copy)

                coefs.append([int(round(eval(elem))) for elem in sols_copy])

            return np.unique(coefs, axis=0).tolist()


# times = []
# for i in range(100):
#     start = time.clock()
#
    # balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C#C[CH2]")], products=[Chem.MolFromSmiles("C#CC"), Chem.MolFromSmiles("[CH3]")])
    # balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("CC(C)(C)C")],
    #           products=[Chem.MolFromSmiles("CC(C)C"), Chem.MolFromSmiles("CC(C)CC(C)(C)C")])
    #
    # times.append(float(time.clock() - start))

# C=CC(=O)OC + [CH2]C + [CH3] ---> C=[C]C + COC + C[C]=O

# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C1CC1")], products=[Chem.MolFromSmiles("C1CCC1"), Chem.MolFromSmiles("CC")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C=CC(=O)OC"), Chem.MolFromSmiles("[CH2]C"), Chem.MolFromSmiles("[CH3]")], products=[Chem.MolFromSmiles("C=[C]C"), Chem.MolFromSmiles("COC"), Chem.MolFromSmiles("C[C]=O")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C1=CC1")], products=[Chem.MolFromSmiles("C=CC"), Chem.MolFromSmiles("C1CC1")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("[CH3]"), Chem.MolFromSmiles("C#C")], products=[Chem.MolFromSmiles("[C]#C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("[CH3]"), Chem.MolFromSmiles("[CH]1CC1")], products=[Chem.MolFromSmiles("[CH2]C")]))

# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C=CC")], products=[Chem.MolFromSmiles("C=C"), Chem.MolFromSmiles("CC")]))
# start = time.clock()
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C=C=C")], products=[Chem.MolFromSmiles("C=C")], atomList=[]))
# end = time.clock()
# print("Execution time:", (end - start), "s.")

# ['CC', '[CH3]']
# ['C=[C]C', '[C]1=CC1', 'C1CCCC1', 'C']

# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("COC(=O)C=C"), Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CO")], products=[Chem.MolFromSmiles("CC=C"), Chem.MolFromSmiles("COC"), Chem.MolFromSmiles("CC(=O)O")]))

# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("[CH3]")], products=[Chem.MolFromSmiles("C=[C]C"), Chem.MolFromSmiles("[C]1=CC1"), Chem.MolFromSmiles("C1CCCC1"), Chem.MolFromSmiles("C")], atomList=[]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C=CC(=O)OC"), Chem.MolFromSmiles("[CH2]C=C"), Chem.MolFromSmiles("CC=O"), Chem.MolFromSmiles("C=[C]C")], products=[Chem.MolFromSmiles("C#CC"), Chem.MolFromSmiles("CC")], atomList=[]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C=CC(=O)OC"), Chem.MolFromSmiles("[CH3]")], products=[Chem.MolFromSmiles("[CH]1CC1"), Chem.MolFromSmiles("COC"), Chem.MolFromSmiles("C=C[O]"), Chem.MolFromSmiles("C[C]=O")], atomList=[]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C=CC(=O)OC"), Chem.MolFromSmiles("C"), Chem.MolFromSmiles("[CH]1CO1")], products=[Chem.MolFromSmiles("C#CC"), Chem.MolFromSmiles("C[C]=O"), Chem.MolFromSmiles("CCO")], atomList=[]))

# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C#CC")], products=[Chem.MolFromSmiles("C#C"), Chem.MolFromSmiles("CC(C)C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C#CC")], products=[Chem.MolFromSmiles("C#C"), Chem.MolFromSmiles("CCCC")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C#CC")], products=[Chem.MolFromSmiles("C#C"), Chem.MolFromSmiles("CC(C)(C)C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C#C[CH2]")], products=[Chem.MolFromSmiles("CCCCCC"), Chem.MolFromSmiles("[C]#C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C1CCCC1"), Chem.MolFromSmiles("C=CC")], products=[Chem.MolFromSmiles("C=C"), Chem.MolFromSmiles("CC=C(C)C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C1CCCC1"), Chem.MolFromSmiles("C=CC")], products=[Chem.MolFromSmiles("C=CCC"), Chem.MolFromSmiles("CC=C(C)C")]))
# print(BalanceEq().balanceEq(reactants=[Chem.MolFromSmiles("C"), Chem.MolFromSmiles("C1=CC1")], products=[Chem.MolFromSmiles("C=CCC")]))
# C + 3 C#CC ---> 3 C#C + CC(C)C
# C + 3 C#CC ---> 3 C#C + CCCC
# C + 4 C#CC ---> 4 C#C + CC(C)(C)C
# C + 5 C#C[CH2] ---> CCCCCC + 5 [C]#C
# C + C1 = CC1 - --> C = CCC
# C1CCCC1 + 4 C=CC ---> C=C + 3 CC=C(C)C
# C1CCCC1 + 3 C=CC ---> C=CCC + 2 CC=C(C)C
# print "Execution time:", np.mean(np.asarray(times)), "s."
# print("Execution time:", (time.clock() - start), "s.")
# raw_input("Waiting")