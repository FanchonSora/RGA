def printResultSame(reactants, products, fout):

    rSize = len(reactants)
    pSize = len(products)

    for i in range(rSize):

        if i > 0:
            fout.write(" + ")

        fout.write(reactants[i])

    fout.write(" ---> ")

    for i in range(pSize):

        if i > 0:
            fout.write(" + ")

        fout.write(products[i])

    fout.write("\n")


def printResultSameParallel(reactants, products):
    result = ""

    rSize = len(reactants)
    pSize = len(products)

    for i in range(rSize):

        if i > 0:
            result += " + "

        result += reactants[i]

    result += " ---> "

    for i in range(pSize):

        if i > 0:
            result += " + "

        result += products[i]

    # result += "\n"
    return result

def printResultDiff(reactants, products, coefs, fout):

    rSize = len(reactants)
    pSize = len(products)

    for i in range(rSize):

        if i > 0:
            fout.write(" + ")

        curR = reactants[i]
        curCoef = coefs[i]
        if curCoef > 1:
            fout.write(str(curCoef) + " ")

        fout.write(curR)

    fout.write(" ---> ")

    for i in range(pSize):

        if i > 0:
            fout.write(" + ")

        curP = products[i]
        curCoef = coefs[i + rSize]

        if curCoef > 1:
            fout.write(str(curCoef) + " ")

        fout.write(curP)
    fout.write("\n")

    # for i in range(pSize):
    #
    #     if i > 0:
    #         fout.write(" + ")
    #
    #     curP = products[i]
    #     curCoef = coefs[i + rSize]
    #     if curCoef > 1:
    #         fout.write(str(curCoef) + " ")
    #
    #     fout.write(curP)
    #
    # fout.write(" ---> ")
    #
    # for i in range(rSize):
    #
    #     if i > 0:
    #         fout.write(" + ")
    #
    #     curR = reactants[i]
    #     curCoef = coefs[i]
    #     if curCoef > 1:
    #         fout.write(str(curCoef) + " ")
    #
    #     fout.write(curR)
    # fout.write("\n")

def printResultDiffParallel(reactants, products, coefs):

    rSize = len(reactants)
    pSize = len(products)

    results = ""

    for i in range(rSize):

        if i > 0:
            results += " + "

        curR = reactants[i]
        curCoef = coefs[i]
        if curCoef > 1:
            results += str(curCoef) + " "

        results += curR

    results += " ---> "

    for i in range(pSize):

        if i > 0:
            results += " + "

        curP = products[i]
        curCoef = coefs[i + rSize]

        if curCoef > 1:
            results += str(curCoef) + " "

        results += curP

    # results += "\n"

    return results
