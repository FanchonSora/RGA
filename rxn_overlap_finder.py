rxnFile1 = "res_test.out"
rxnFile2 = "res_test2.out"

rxnList = {}

with open(rxnFile1) as f1, open(rxnFile2) as f2:

    dups = 0
    isodesmicCnt = 0

    for line in f2:
        line = line.rstrip("\n").split("\t")[0]
        rxnList[line] = 1
        isodesmicCnt += 1

    cnt = 0

    for line in f1:

        if cnt == isodesmicCnt:
            break

        line = line.rstrip("\n").split("\t")[0]

        cnt += 1

        if line in rxnList:
            dups += 1

    print(isodesmicCnt)
    print(dups, dups / (isodesmicCnt * 1.0))