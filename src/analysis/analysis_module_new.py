from utility_module import np, math, operator

def analyzeReactions(fout, r1, reactions, inputCalcData, uncerts, expts, calcs, config):
    topReactions = config.top
    print("Start analyzing reactions...")
    sorted_reactions = sorted(reactions.items(), key=operator.itemgetter(1), reverse=True)
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
    rxnList = []
    sz = len(sorted_reactions)
    for i in range(sz):
        line = sorted_reactions[i][0]
        line = line.rstrip("\n").strip()
        toks = line.split("--->")
        lSide, rSide = toks[0].split(" + "), toks[1].split(" + ")
        exptHof = 0
        totalUncert = 0.0
        for r in rSide:
            toks = r.split()
            if len(toks) == 2:
                coef = int(toks[0])
                val = toks[1]
                exptHof += coef * (expts[val] - calcs[val])
                totalUncert += (coef * uncerts[val] * 1.0) ** 2
            else:
                val = toks[0]
                exptHof += expts[val] - calcs[val]
                totalUncert += uncerts[val] ** 2
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
                    totalUncert += (coef * uncerts[val] * 1.0) ** 2
            else:
                if index == 0:
                    exptHof += inputCalcData
                else:
                    val = toks[0]
                    exptHof += calcs[val] - expts[val]
                    totalUncert += uncerts[val] ** 2
        exptHof /= inputCoef * 1.0
        totalUncert = math.sqrt(totalUncert)
        totalUncert /= inputCoef * 1.0
        exptData.append(exptHof)
        rxnList.append([sorted_reactions[i][0], sorted_reactions[i][1], exptHof, totalUncert])
        uncertData.append(totalUncert)
    print("End analyzing reactions...\n")
    rxnList = np.asarray(rxnList)
    exptData = np.asarray(exptData)
    uncertData = np.asarray(uncertData)
    leftWidth = 50
    print("Species".ljust(leftWidth) + ": " + r1)
    print("Calculated Heat of Formation with all reactions found:")
    print("\t+ Number of reactions".ljust(leftWidth) + ": " + str(len(rxnList)))
    print("\t+ Average HoF at 298 K (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.mean(exptData)))
    print("\t+ Standard deviation (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.std(exptData)))
    print("\t+ Average uncertainty (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.mean(uncertData)) + "\n")
    fout.write("SUMMARY\n\n")
    fout.write("Calculated Heat of Formation with all reactions found:\n")
    fout.write("\t+ Number of reactions".ljust(leftWidth) + ": " + str(len(rxnList)) + "\n")
    fout.write("\t+ Average HoF at 298 K (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.mean(exptData)) + "\n")
    fout.write("\t+ Standard deviation (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.std(exptData)) + "\n")
    fout.write("\t+ Average uncertainty (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.mean(uncertData)) + "\n\n")
    maxWidth = 15 * config.noSpecies
    sz = len(sorted_reactions)
    maxCnt = len(str(sz))
    selected_mean_expt = np.mean(exptData)
    selected_std_expt = np.std(exptData)
    selectedUncertData = uncertData
    selectedExptData = exptData
    if config.filter or config.top != -1:
        selectedRxnList = rxnList[:topReactions]
        selectedExptData = exptData[:topReactions]
        selectedUncertData = uncertData[:topReactions]
        selected_mean_expt = np.mean(selectedExptData)
        selected_std_expt = np.std(selectedExptData)
        cut_off = 2
        selectedIndices = np.where(np.logical_and(selectedExptData <= selected_mean_expt + cut_off * selected_std_expt,
                                                                   selectedExptData >= selected_mean_expt - cut_off * selected_std_expt))
        selectedExptData = selectedExptData[selectedIndices]
        selectedRxnList = selectedRxnList[selectedIndices]
        selectedUncertData = selectedUncertData[selectedIndices]
        selected_mean_expt = np.mean(selectedExptData)
        selected_std_expt = np.std(selectedExptData)
        print("Calculated Heat of Formation with top-similarity reactions:")
        print("\t+ Number of reactions".ljust(leftWidth) + ": " + str(len(selectedRxnList)))
        print("\t+ Average HoF at 298 K (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % selected_mean_expt))
        print("\t+ Standard deviation (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % selected_std_expt))
        print("\t+ Average uncertainty (kcal/mol)".ljust(leftWidth) + ": " + str.format(
            "%.3f" % np.mean(selectedUncertData)) + "\n")
        fout.write("Calculated Heat of Formation with top-similarity reactions:\n")
        fout.write("\t+ Number of reactions".ljust(leftWidth) + ": " + str(len(selectedRxnList)) + "\n")
        fout.write("\t+ Average HoF at 298 K (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % selected_mean_expt) + "\n")
        fout.write("\t+ Standard deviation (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % selected_std_expt) + "\n")
        fout.write("\t+ Average uncertainty (kcal/mol)".ljust(leftWidth) + ": " + str.format("%.3f" % np.mean(selectedUncertData)) + "\n\n")
        fout.write("=========================================================================\n")
        fout.write("LIST OF REACTIONS FOUND\n\n")
        fout.write("%*s\t\t%*s\t\t%s\t\t%s\t\t%s\n\n" % (maxCnt, "No.", maxWidth, "Reaction", "Similarity", "HoF (298 K)", "Uncertainty"))
        selectedSz = len(selectedRxnList)
        for index in range(sz):
            curRxn = rxnList[index]
            line = curRxn[0]
            rxnSim = float(curRxn[1])
            hof = float(curRxn[2])
            uncert = float(curRxn[3])
            fout.write("%*d\t\t%*s\t\t%.5f\t\t\t%.5f\t\t\t%.5f\n" % (maxCnt, index + 1, maxWidth, line, rxnSim, hof, uncert))
            if index == selectedSz - 1:
                fout.write("\n----------------------------------------------- End of top " + str(len(selectedRxnList)) + " reactions ----------------------------------------\n\n")
    else:
        fout.write("=========================================================================\n")
        fout.write("LIST OF REACTIONS FOUND\n\n")
        fout.write("%*s\t\t%*s\t\t%s\t\t%s\t\t%s\n\n" % (maxCnt, "No.", maxWidth, "Reaction", "Similarity", "HoF (298 K)", "Uncertainty"))
        for index in range(sz):
            curRxn = rxnList[index]
            line = curRxn[0]
            rxnSim = float(curRxn[1])
            hof = float(curRxn[2])
            uncert = float(curRxn[3])
            fout.write("%*d\t\t%*s\t\t%.5f\t\t\t%.5f\t\t\t%.5f\n" % (maxCnt, index + 1, maxWidth, line, rxnSim, hof, uncert))
    fout.write("\n=========================================================================\n")
    print("Mean value of experimental HoF of " + r1 + ":")
    print("%.3f" % selected_mean_expt)
    print("Standard deviation of experimental HoF of " + r1 + ":")
    print("%.3f" % selected_std_expt)
    print("Mean value of uncertainty of " + r1 + ":")
    print("%.3f" % math.fabs(np.mean(selectedUncertData)))