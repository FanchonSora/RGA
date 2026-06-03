from utility_module import np, math, operator

def analyzeReactions(resFile, r1, reactions, inputCalcData, uncerts, expts, calcs, config):

    topReactions = config.top
    # print reactions

    print("Start analyzing reactions...")

    sorted_reactions = sorted(reactions.items(), key=operator.itemgetter(1), reverse=True)

    # sims = [sorted_reactions[i][1] for i in range(len(sorted_reactions))]

    # topReactions = len(sorted_reactions)

    if topReactions == -1:
        topReactions = len(sorted_reactions)
    elif topReactions < 1.0:
        topReactions = int(math.ceil(topReactions * len(sorted_reactions)))
    elif topReactions > len(sorted_reactions):
        topReactions = len(sorted_reactions)
    else:
        topReactions = int(topReactions)

    exptData = []
    uncertData = []
    # sumUncertData = []
    # sims = []

    # rxnSims = []
    # reactantSims = []
    # productSims = []
    # diffs = []

    rxnList = []

    # with open(resFile, "w") as fout:

    # fout.write("The top " + str(topReactions) + " isodesmic reactions with highest similarity for input " + r1 + " :\n")

    for i in range(topReactions):

        # if sorted_reactions[i][1] < 0.7:
        #     continue

        # sims.append(sorted_reactions[i][1])

        line = sorted_reactions[i][0]

        # fout.write(sorted_reactions[i][0] + "\t" + str(sorted_reactions[i][1]) + "\n")
        # fout.write(sorted_reactions[i][0] + "\t" + str(sorted_reactions[i][1]) + "\t")
        # fout.write("-------------\n")

        # print(sorted_reactions[i][0], str(sorted_reactions[i][1]),)
        # print "-------------"

        line = line.rstrip("\n").strip()

        toks = line.split("--->")
        lSide, rSide = toks[0].split(" + "), toks[1].split(" + ")

        exptHof = 0
        totalUncert = 0.0
        # sumTotalUncert = 0.0

        # print line

        for r in rSide:
            toks = r.split()
            if len(toks) == 2:
                coef = int(toks[0])
                val = toks[1]
                exptHof += coef * (expts[val] - calcs[val])
                # totalUncert += coef * uncerts[val] * 1.0
                totalUncert += (coef * uncerts[val] * 1.0) ** 2
                # sumTotalUncert += coef * uncerts[val] * 1.0
            else:
                val = toks[0]
                exptHof += expts[val] - calcs[val]
                # totalUncert += uncerts[val]
                totalUncert += uncerts[val] ** 2
                # sumTotalUncert += uncerts[val]

        inputCoef = 1

        for index in range(len(lSide)):
            l = lSide[index]
            toks = l.split()
            if len(toks) == 2:
                coef = int(toks[0])
                if index == 0:
                    inputCoef = coef
                    exptHof += coef * inputCalcData
                else:
                    val = toks[1]
                    exptHof += coef * (calcs[val] - expts[val])
                    # totalUncert += coef * uncerts[val] * 1.0
                    # sumTotalUncert += coef * uncerts[val] * 1.0

                    # totalUncert -= coef * uncerts[val] * 1.0
                    totalUncert += (coef * uncerts[val] * 1.0) ** 2
            else:
                if index == 0:
                    exptHof += inputCalcData
                else:
                    val = toks[0]
                    exptHof += calcs[val] - expts[val]
                    # totalUncert += uncerts[val]
                    # sumTotalUncert += uncerts[val]

                    # totalUncert -= uncerts[val]
                    totalUncert += uncerts[val] ** 2

        exptHof /= inputCoef * 1.0
        totalUncert /= inputCoef * 1.0

        exptData.append(exptHof)

        # fout.write(sorted_reactions[i][0] + "\t")

        # simList = sorted_reactions[i][1]
        # # print(simList)
        #
        # for index in range(len(simList)):
        #     pass
            # print(simList[index])
            # fout.write(str(simList[index]) + "\t")

        # fout.write(str(abs(exptHof - expts[r1])) + "\n")
        # fout.write(str(exptHof) + "\n")
        # fout.write("-------------\n")
        # fout.write("\n")

        totalUncert = math.sqrt(totalUncert)
        rxnList.append([sorted_reactions[i][0], sorted_reactions[i][1], exptHof, totalUncert])

        # rxnSims.append(simList[0])
        # reactantSims.append(simList[1])
        # productSims.append(simList[2])

        # diffs.append(abs(exptHof - expts[r1]))

        # print(exptHof)

        uncertData.append(totalUncert)
        # sumUncertData.append(sumTotalUncert)

    print("End analyzing reactions...\n")

    # if np.mean(diffs) > 1:
    #     print("Species with caution:", r1)

    # plt.plot(rxnSims, diffs, color='r', label='Rxn sim')
    # plt.plot(reactantSims, diffs, color='g', label='Reactant Sim')
    # plt.plot(productSims, diffs, color='b', label='Product Sim')
    # # plt.plot(sims, diffs, color='b')
    # plt.legend()
    # plt.show()

    # print "Values of experimental HoF of", r1, ":"
    # fout.write("Values of experimental HoF of " + r1 + " :\n")
    # for val in exptData:
    #     fout.write("%.3f\t" % val)
    #     print "%.3f\t" % val,
    # fout.write("\n\n")
    # print
    #
    # print "Values of uncertainty of", r1, ":"
    # fout.write("Values of uncertainty of " + r1 + " :\n")
    # for val in uncertData:
    #     fout.write("%.3f\t" % val)
    #     print "%.3f\t" % val,
    # fout.write("\n\n")
    # print

    rxnList = np.asarray(rxnList)
    exptData = np.asarray(exptData)
    uncertData = np.asarray(uncertData)
    sims = np.asarray(sims)

    # while True:
    #
    #     prevLen = len(exptData)

        # print(mean_expt, std_expt)

    mean_expt = np.mean(exptData)
    std_expt = np.std(exptData)

    if config.filter:

        cut_off = 2

        selectedIndices = np.where(np.logical_and(exptData <= mean_expt + cut_off * std_expt,
                                                           exptData >= mean_expt - cut_off * std_expt))


        # sims = sims[selectedIndices]
        exptData = exptData[selectedIndices]
        rxnList = rxnList[selectedIndices]
        uncertData = uncertData[selectedIndices]

        mean_expt = np.mean(exptData)
        std_expt = np.std(exptData)

    with open(resFile, "w") as fout:
        fout.write("The top " + str(len(rxnList)) + " isodesmic reactions with highest similarity for input " + r1 + " :\n")

        maxWidth = 15 * config.noSpecies

        for rxn in rxnList:
            fout.write("%*s\t\t%.5f\t\t%.5f\t|\n" % (maxWidth, rxn[0], float(rxn[1]), float(rxn[2])))
            fout.write("-"*(maxWidth + 33) + "\n")
    #
    # # print(len(exptData))
    # # print(exptData)
    #
    # # if math.fabs(mean_expt - expts[r1]) > 1:
    # #     print("Alert:", r1)
    #
    # # print("Real expt value of", r1, ":")
    # # fout.write("Real experimental HoF of " + r1 + " :\n")
    # # print(expts[r1])
    # # fout.write("%.3f" % expts[r1] + "\n\n")
    #
    # # if math.fabs(mean_expt - expts[r1]) > 1:
    # #     print("Alert", r1)
    # #
    # # print("Mean value of experimental HoF of", r1, ":")
    #     fout.write("\nMean value of experimental HoF of " + r1 + " :\n")
    #     # print("%.3f" % np.mean(np.asarray(exptData)))
    #     fout.write("%.3f\n\n" % mean_expt)
    #     #
    #     # print("Standard deviation of experimental HoF of", r1, ":")
    #     # print("%.3f" % np.std(np.asarray(exptData)))
    #
    #     fout.write("Standard deviation of experimental HoF of " + r1 + " :\n")
    #     fout.write("%.3f\n\n" % std_expt)
    #     #
    #     # print("Mean value of uncertainty of", r1, ":")
    #     fout.write("Mean value of uncertainty of " + r1 + " :\n")
    #     # print("%.3f" % np.mean(np.asarray(uncertData)))
    #     # print("--------------------------------\n")
    #     fout.write("%.3f\n\n" % np.mean(np.asarray(uncertData)))


    # diffs = np.abs(np.asarray(exptData) - expts[r1])
    # diffs = np.asarray(exptData) - expts[r1]

    # plt.hist(exptData, bins=200)
    # plt.title('Distribution of HoF value of ' + r1)
    # plt.show()

    # mean = np.mean(np.asarray(exptData))
    # std = np.std(np.asarray(exptData))
    #
    # std = 1
    #
    # coef = 1
    # upBound = mean + coef * std
    # lowBound = mean - coef * std
    #
    # exptData = np.asarray(exptData)
    #
    # exptData = exptData[exptData <= upBound]
    # exptData = exptData[exptData >= lowBound]
    #
    # print("Len of exptData:", len(exptData))
    #
    # print("New mean", np.mean(exptData))

    # diffs = np.abs(np.asarray(exptData) - expts[r1])
    # diffs = np.asarray(exptData)
    # # diffs = np.asarray(uncertData)
    #
    # import matplotlib.pyplot as plt
    # sims = sims[:topReactions]
    #
    # print sims
    # #
    # plt.plot(sims, exptData, color='r', label='HoF value')
    # # plt.plot(sims, diffs, color='r')
    # # # plt.plot(sims, diffs, color='b')
    # # # plt.plot(sims, uncertData, color='b', label='Uncertainty')
    # plt.axhline(y=expts[r1], color='b', linestyle='--', label="Experimental data")
    # # plt.axhline(y=mean_expt, color='b', linestyle='-.', label="Mean predicted value")
    # # # plt.axhline(y=uncerts[r1], color='g', linestyle='-')
    # plt.xlabel('Reaction similarity')
    # plt.ylabel('Value (kcal/mol)')
    # # plt.title('Relationship between HoF and similarity for isodesmic rxn of ' + r1)
    # # # plt.title('Relationship between uncertainty and similarity for isodesmic rxn of ' + r1)
    # plt.title('Similarity vs. HoF value (all rxns) for ' + r1)
    # plt.xlim([0, 1])
    # plt.legend()
    # plt.show()

    # meanVal = np.mean(np.asarray(exptData))

    # fout.write("%d\t%.3f\n" % (len(reactions), math.fabs(mean_expt - expts[r1])))

    # print("Number of selected reactions:", len(selectedRxnList))

    print("Mean value of experimental HoF of", r1, ":")
    print("%.3f" % mean_expt)
    # # print("Real expt value of", r1, ":")
    # # print("%.3f" % expts[r1])
    # # print("Difference between predicted and true value of HoF of", r1, ":")
    # # print("%.3f" % math.fabs(mean_expt - expts[r1]))
    print("Standard deviation of experimental HoF of", r1, ":")
    print("%.3f" % std_expt)
    print("Mean value of uncertainty of", r1, ":")
    print("%.3f" % math.fabs(np.mean(np.asarray(uncertData))))
    # print("Mean value of sum uncertainty of", r1, ":")
    # print("%.3f" % np.mean(np.asarray(sumUncertData)))
    # print("Real value of uncertainty of", r1, ":")
    # print("%.3f" % uncerts[r1])

    # return [r1, topReactions, expts[r1], mean_expt, math.fabs(mean_expt - expts[r1]), std_expt, uncerts[r1], np.mean(np.asarray(uncertData))]
    return [r1, len(exptData), expts[r1], mean_expt, math.fabs(mean_expt - expts[r1]), std_expt, uncerts[r1], math.fabs(np.mean(np.asarray(uncertData)))]
    #
    # if math.fabs(meanVal - expts[r1]) > 1:
    #     print "Mean value of experimental HoF of", r1, ":"
    #     print meanVal
    #     print "Real expt value of", r1, ":"
    #     print expts[r1]
    #     print "Mean value of uncertainty of", r1, ":"
    #     print "%.3f" % np.mean(np.asarray(uncertData))
    #     print "Real value of uncertainty of", r1, ":"
    #     print uncerts[r1]
