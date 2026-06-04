from .utils import Chem, defaultdict, np, permutations, combine2Dicts, math, matrix_rank, re, sympy, Matrix, LpProblem, LpMinimize, LpVariable, LpInteger, LpStatus, PulpSolverError, PULP_CBC_CMD

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

        curKey = keys[curIndex]

        for i in range(1, maxCnt + 1):

            if not valueFound and i > maxCnt - step:
                valueFound = True

            oldVal = chars[curKey]
            chars[curKey] = i

            if valueFound:

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

        rLen = len(reactants)
        pLen = len(products)
        totalLen = rLen + pLen

        A = self.processSpecies(rList=reactants, pList=products, atomList=atomList)

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

        mRank = matrix_rank(A)

        if mRank >= totalLen:
            return None

        if mRank == totalLen - 1:

            coefs = sympy.Matrix(A).nullspace()[0]

            if self.checkCoefs(coefs):
                if not self.hasSameElems(coefs):

                    p = re.compile("\\/(\\d+)")

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

            sols = mSols['sols']
            chars = mSols['symbols']
            selectedSols = mSols['selectedSols']

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
                    if maxCnt > cutOffMaxCnt:

                        varList = []
                        objective = None

                        for index in range(totalLen):
                            varList.append(LpVariable(chr(ord('t') + index), 1, 50, LpInteger))

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

                        try:
                            prob.solve(solver=PULP_CBC_CMD(msg=0))

                            if LpStatus[prob.status] == 'Optimal':
                                return [[int(round(v.varValue, 9)) for v in prob.variables()]]

                            return None

                        except PulpSolverError:
                            return None

                sols_copy = ",".join(sols)

                for c in self.results:
                    sols_copy = sols_copy.replace(c, str(self.results[c]))

                sols_copy = sols_copy.split(",")

                coefs.append([int(round(eval(elem))) for elem in sols_copy])

            return np.unique(coefs, axis=0).tolist()
