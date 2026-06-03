from utility_module import time, is_number, Chem, os, freeze_support, printHeader, printExecutionDetails, traceback, \
	urlopen, HTTPError, URLError, exit
# from utility_module import *
from rxngenconfig import Config
# from rxngenerator_complete import generate_rxns_complete
from rxngenerator_complete_parallel import generate_rxns_complete
# from rxngenerator_stochastic import generate_rxns_stochastic
from rxngenerator_stochastic_parallel import generate_rxns_stochastic
# from similarity_module import checkSimilarity
from similarity_parallel import checkSimilarity
# from analysis_module import analyzeReactions
from analysis_module_new import analyzeReactions

def main():

	# res = urlopen('http://just-the-time.appspot.com/')
	# res = urlopen('https://nvd.nist.gov/feeds/json/cve/1.1/nvdcve-1.1-modified.meta')
	# time_str = res.read().strip().splitlines()[0]

	# start_index, end_index = time_str.index(':'), time_str.index('T')
	# time_str = time_str[start_index + 1 : end_index]

	# date = time_str.split()[0]
	# year, month, day = [int(s) for s in date.split("-")]

	# print("\n*** ------------------------------------- ***")
	# print("*** RGA: Reaction Generator and Analyzer  ***")
	# print("*** Copyright @ Triet Le & Lam Huynh " + str(year) + " ***")
	# print("*** ------------------------------------- ***\n")

	# check_expiry = False

	# if check_expiry:

	# 	# Setting the expiry date: Dec 31, 2022

	# 	expiry_day, expiry_month, expiry_year = 31, 12, 2022
	# 	expiry_date = str(expiry_year) + '/' + str(expiry_month) + '/' + str(expiry_day)

	# 	if year > expiry_year or (year == expiry_year and month > expiry_month) or \
	# 			(year == expiry_year and month == expiry_month and day > expiry_day):
	# 		print("Your version has been expired on " + expiry_date + ". Please contact the authors to get an updated version. Thank you !!!")
	# 		raw_input("Press Enter to continue...")
	# 		exit(0)

	# Input species parameters
	config = Config()

	if config.parseConfig() == False:
		print("There is something wrong with the config file (config.yaml or config.json). Please check again !!!")
		raw_input("Press Enter to continue...")
		exit(0)

	speciesFile = config.speciesFile
	# rxnFile = config.rxnFile
	# simFile = config.simFile
	resFile = config.resFile

	species = {}
	expts = {}
	uncerts = {}
	calcs = {}
	inputCalcData = None

	if not os.path.exists(speciesFile):
		print("The species file " + speciesFile + " does not exist. Please check again !!!")
		raw_input("Press Enter to continue...")
		exit(0)

	# speciesFile = "literature_data.csv"
	# speciesFile = "qm9.csv"

	with open(speciesFile) as f:
		for line in f:

			line = line.rstrip("\n")

			if line == "":
				continue

			# toks = line.split()
			toks = line.split()

			if len(toks) == 3:
				uncert = 1
				calc = float(toks[2])
			else:
				uncert = float(toks[2])
				calc = float(toks[3])

			s = toks[0]
			expt = float(toks[1])

			# s = toks[0]
			# expt = float(toks[1])
			# uncert = float(toks[2])
			# calc = float(toks[3])

			# s = toks[-1].replace("\"", "")
			# expt = float(toks[1].replace("\"", ""))
			# calc = float(toks[2].replace("\"", ""))
			# uncert = float(toks[3].replace("\"", ""))

			# s = toks[1]
			# expt = calc = uncert = float(toks[2])

			# species[line] = 1
			species[s] = 1
			expts[s] = expt
			uncerts[s] = uncert
			calcs[s] = calc

	# keys = sorted(species.keys())

	keys = list(species.keys())

	print("Number of species in the database: " + str(len(keys)) + "\n")

	# ------------------------------------

	executor = config.executor

	if config.executor == "":
		executor = raw_input("Please enter the name of the executor: ")

	if config.speciesSmiles == None:
		input_species = [""]
		calc_vals = [""]
	else:
		input_species = str(config.speciesSmiles).split()
		calc_vals = str(config.calcValue).split()

	parentResFile = resFile

	for index in range(len(input_species)):

		s1 = input_species[index]

		if s1 == "":
			s1 = raw_input("Input the SMILES of the species: ")

			m1 = Chem.MolFromSmiles(s1)

			while s1 == "" or not m1:
				print("The input SMILES is invalid. Please input again or enter \"exit\" to close the program")
				s1 = raw_input("")
				if s1.lower() == "exit":
					raw_input("Press Enter to continue...")
					exit(0)
				m1 = Chem.MolFromSmiles(s1)

			s1_standard = Chem.MolToSmiles(m1, isomericSmiles=True)

			# if s1 != s1_standard:
			#     print("The standardized form of", s1, "is", s1_standard)
			#     print(s1_standard, "will be used instead to ensure the consistency\n")
			#     s1 = s1_standard

			if s1 not in species and s1_standard not in species:
				print("This species is not in the database")
				inputCalcData = raw_input("Please enter the calculated value for the input species: ")
				while inputCalcData == "" or not is_number(inputCalcData):
					print("The input calculated value is invalid (not a number). Please input again or enter \"exit\" to close the program")
					inputCalcData = raw_input("")
					if inputCalcData.lower() == "exit":
						raw_input("Press Enter to continue...")
						exit(0)
				inputCalcData = float(inputCalcData)

			else:
				if s1 in calcs:
					inputCalcData = calcs[s1]
				else:
					inputCalcData = calcs[s1_standard]
		else:
			m1 = Chem.MolFromSmiles(s1)

			print("\n########################################\n")
			print("Processing " + s1)

			if not m1:
				print("The input SMILES is invalid. Please double-check the config file.")
				raw_input("Press Enter to continue...")
				exit(0)

			inputCalcData = float(calc_vals[index])

		resFile = s1.replace("/", "'").replace("\\", "`") + "_" + parentResFile

		# ------------------------------------

		with open(resFile, "w") as fout:

			fout = printHeader(fout, s1, executor, config, year)

			start = time.clock()

			if config.stochastic:
				# rxnCnt = generate_rxns_stochastic(keys, s1, m1, rxnFile, config)
				rxnCnt, rxnList = generate_rxns_stochastic(keys, s1, m1, fout, config)
			else:
				# rxnCnt = generate_rxns_complete(keys, s1, m1, rxnFile, config)
				rxnCnt, rxnList = generate_rxns_complete(keys, s1, m1, fout, config)

			# -------------------------------------------------

			if config.similarityOn:
				if rxnCnt > 0:
					# reactions = checkSimilarity(rxnCnt=rxnCnt, rxnFile=rxnFile, simFile=simFile, config=config)
					reactions = checkSimilarity(rxnList=rxnList, config=config, fout=fout)

				# -------------------------------------------------

					if config.analysisOn:
						# analyzeReactions(resFile=resFile, r1=s1, reactions=reactions, inputCalcData=inputCalcData, uncerts=uncerts, expts=expts, calcs=calcs, config=config)
						analyzeReactions(fout=fout, r1=s1, reactions=reactions, inputCalcData=inputCalcData,
										 uncerts=uncerts, expts=expts, calcs=calcs, config=config)

				# else:
				#     with open(simFile, "w") as fout1:
				#         fout1.write("There is no reaction found for species " + s1 + "\n")
				#
				#     if config.analysisOn:
				#         with open(resFile, "w") as fout2:
				#             fout2.write("There is no reaction found for species " + s1 + "\n")

			# --------------------------------------------------

			exe_time = time.clock() - start
			print("Execution time: " + str(exe_time) + " s.")

			printExecutionDetails(fout, exe_time)

	raw_input("Press Enter to continue...")

if __name__ == '__main__':

	freeze_support()

	try:
		main()
	except (HTTPError, URLError) as e:

		print("\n*** ------------------------------------- ***")
		print("*** RGA: Reaction Generator and Analyzer  ***")
		print("***   Copyright @ Triet Le & Lam Huynh    ***")
		print("*** ------------------------------------- ***\n")
		print("Please check your Internet connection or contact the authors for more information. Thank you !!!\n")

		with open("error.out", "w") as fout:
			fout.write("Please check your Internet connection or contact the authors for more information. Thank you !!!\n")

		raw_input("Press Enter to continue...")
		exit(0)

	except Exception as e:

		print("There is something wrong with the program. Please check the details in error.out file. Thank you !!!\n")

		with open("error.out", "w") as fout:
			fout.write(str(traceback.format_exc()) + "\n")

		raw_input("Press Enter to continue...")
		exit(0)
