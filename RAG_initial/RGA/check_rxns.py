from itertools import permutations

rxnFile = "rxn_test.out"

with open(rxnFile) as f:

    dups = {}

    for line in f:

        line = line.rstrip("\n").strip()

        toks = line.split("--->")
        lSide, rSide = toks[0].split(" + "), toks[1].split(" + ")

        reactantsKey = []
        productsKey = []

        for r in rSide:
            toks = r.split()
            if len(toks) == 2:
                s = toks[1]
            else:
                s = toks[0]

            productsKey.append(s)

        for l in lSide:
            toks = l.split()
            if len(toks) == 2:
                s = toks[1]
            else:
                s = toks[0]

            reactantsKey.append(s)

        print(len(reactantsKey), len(productsKey))

        rPermutations = list(permutations(reactantsKey[1:]))
        pPermutations = list(permutations(productsKey))

        for rPermutation in rPermutations:
            for pPermutation in pPermutations:
                curCombinedKeys = ','.join(rPermutation + pPermutation)
                dups[curCombinedKeys] = 1

print(len(dups.keys()))