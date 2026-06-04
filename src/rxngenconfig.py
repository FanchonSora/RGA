from utility_module import json, yaml, os

class Config:
    def __init__(self):
        self.filesKW = "files"
        self.speciesFileKW = "species_file"
        self.rxnFileKW = "rxn_file"
        self.simFileKW = "sim_file"
        self.resFileKW = "res_file"

        # For testing
        self.executorKW = "executor"
        self.speciesSmilesKW = "species_smiles"
        self.calcValueKW = "calc_value"

        self.rxnContraintsKW = "rxn_constraints"
        self.noProcessorsKW = "no_of_cores"
        self.sameCoefKW = "same_coef"
        self.noSpeciesKW = "no_species"
        self.completeSetKW = "complete_set"
        self.stochasticKW = "stochastic"
        self.noTrialsKW = "no_trials"

        # For testing
        self.seedKW = "seed"

        self.lElectronKW = "lone_electron"
        self.ePairKW = "electron_pair"
        self.ePairTypeKW = "electron_pair_type"
        self.atomWithHydroKW = "atom_with_hydro"
        self.atomHybridizationKW = "atom_hybridization"
        self.atomHybridWithHydroKW = "atom_hybrid_with_hydro"
        self.atomDegreeKW = "atom_degree"
        self.atomTypeKW = "atom_type"
        self.radicalKW = "radical"

        self.bondtypeKW = "bond_type"
        self.hydroBondKW = "hydro_bond"
        self.normalBondKW = "normal_bond"
        self.hybridizationKW = "bond_hybridization"
        self.cHydroKW = "bond_with_hydro"
        self.cDegreeKW = "bond_degree"
        self.bondRadicalKW = "bond_radical"
        self.bondAtomTypeKW = "bond_atom_type"
        self.ringKW = "ring"

        self.similarityKW = "similarity"
        self.similarityOnKW = "similarity_on"
        self.radiusKW = "radius"
        self.chiralityKW = "chirality"
        self.featureMethodKW = "feature_method"

        self.analysisKW = "analysis"
        self.analysisOnKW = "analysis_on"
        self.filterKW = "filter"
        self.topKW = "top"

    def decodeBool(self, s):
        if s.lower() == "true":
            return True

        return False

    def parseConfig(self):

        print("Start parsing config file...")

        parseResult = False

        if os.path.exists("config.yaml"):
            parseResult = self.parseYAML()
        elif os.path.exists("config.json"):
            parseResult = self.parseJSON()

        print("End parsing config file...")
        print("--------------------------------\n")

        return parseResult

    def parseJSON(self, fileName="config.json"):

        with open(fileName) as f:
            try:
                self.jsonfile = json.load(f)
            except (ValueError, KeyError) as e:
                return False

        self.parseFiles(self.jsonfile[self.filesKW])
        self.parseContraints(self.jsonfile[self.rxnContraintsKW])
        self.parseSimInfo(self.jsonfile[self.similarityKW])
        self.parseAnalysisInfo(self.jsonfile[self.analysisKW])

        return True


    def parseYAML(self, fileName="config.yaml"):

        with open(fileName) as f:
            try:
                self.configFile = yaml.safe_load(f)
            except yaml.YAMLError as e:
                return False

        self.parseFiles(self.configFile[self.filesKW])
        self.parseContraints(self.configFile[self.rxnContraintsKW])
        self.parseSimInfo(self.configFile[self.similarityKW])
        self.parseAnalysisInfo(self.configFile[self.analysisKW])

        return True

    def parseFiles(self, filesInfo):
        self.speciesFile = filesInfo[self.speciesFileKW]
        # self.rxnFile = filesInfo[self.rxnFileKW]
        # self.simFile = filesInfo[self.simFileKW]
        self.resFile = filesInfo[self.resFileKW]

        # Used for testing
        self.executor = filesInfo[self.executorKW]
        self.speciesSmiles = filesInfo[self.speciesSmilesKW]
        self.calcValue = filesInfo[self.calcValueKW]

    def parseContraints(self, contraints):
        # self.sameCoef = self.decodeBool(contraints[self.sameCoefKW])
        self.noProcessors = int(contraints[self.noProcessorsKW])
        self.sameCoef = contraints[self.sameCoefKW]
        self.noSpecies = int(contraints[self.noSpeciesKW])
        # self.completeSet = self.decodeBool(contraints[self.completeSetKW])
        self.completeSet = contraints[self.completeSetKW]
        # self.stochastic = self.decodeBool(contraints[self.stochasticKW])
        self.stochastic = contraints[self.stochasticKW]
        self.noTrials = int(contraints[self.noTrialsKW])

        self.seed = contraints[self.seedKW] # Empty or int

        # self.lElectron = self.decodeBool(contraints[self.lElectronKW])
        self.lElectron = contraints[self.lElectronKW]
        # self.ePair = self.decodeBool(contraints[self.ePairKW])
        self.ePair = contraints[self.ePairKW]
        # self.ePairType = self.decodeBool(contraints[self.ePairTypeKW])
        self.ePairType = contraints[self.ePairTypeKW]
        # self.atomWithHydro = self.decodeBool(contraints[self.atomWithHydroKW])
        self.atomWithHydro = contraints[self.atomWithHydroKW]
        # self.atomHybridization = self.decodeBool(contraints[self.atomHybridizationKW])
        self.atomHybridization = contraints[self.atomHybridizationKW]
        # self.atomHybridWithHydro = self.decodeBool(contraints[self.atomHybridWithHydroKW])
        self.atomHybridWithHydro = contraints[self.atomHybridWithHydroKW]
        # self.atomDegree = self.decodeBool(contraints[self.atomDegreeKW])
        self.atomDegree = contraints[self.atomDegreeKW]
        self.atomType = str(contraints[self.atomTypeKW])
        # self.radical = self.decodeBool(contraints[self.radicalKW])
        self.radical = contraints[self.radicalKW]

        # self.hydroBond = self.decodeBool(contraints[self.bondtypeKW][self.hydroBondKW])
        # self.normalBond = self.decodeBool(contraints[self.bondtypeKW][self.normalBondKW])
        # self.bondHybrid = self.decodeBool(contraints[self.bondtypeKW][self.hybridizationKW])
        # self.bondHydro = self.decodeBool(contraints[self.bondtypeKW][self.cHydroKW])
        # self.bondDegree = self.decodeBool(contraints[self.bondtypeKW][self.cDegreeKW])
        # self.bondRadical = self.decodeBool(contraints[self.bondtypeKW][self.bondRadicalKW])

        self.hydroBond = contraints[self.bondtypeKW][self.hydroBondKW]
        self.normalBond = contraints[self.bondtypeKW][self.normalBondKW]
        self.bondHybrid = contraints[self.bondtypeKW][self.hybridizationKW]
        self.bondHydro = contraints[self.bondtypeKW][self.cHydroKW]
        self.bondDegree = contraints[self.bondtypeKW][self.cDegreeKW]
        self.bondRadical = contraints[self.bondtypeKW][self.bondRadicalKW]

        self.bondAtomType = str(contraints[self.bondtypeKW][self.bondAtomTypeKW])
        # self.ring = self.decodeBool(contraints[self.bondtypeKW][self.ringKW])
        self.ring = contraints[self.bondtypeKW][self.ringKW]

    def parseSimInfo(self, simInfo):
        # self.similarityOn = self.decodeBool(simInfo[self.similarityOnKW])
        self.similarityOn = simInfo[self.similarityOnKW]
        self.radius = int(simInfo[self.radiusKW])
        # self.chirality = self.decodeBool(simInfo[self.chiralityKW])
        self.chirality = simInfo[self.chiralityKW]
        self.featureMethod = str(simInfo[self.featureMethodKW])

    def parseAnalysisInfo(self, analysisInfo):
        # self.analysisOn = self.decodeBool(analysisInfo[self.analysisOnKW])
        self.analysisOn = analysisInfo[self.analysisOnKW]
        # self.filter = self.decodeBool(analysisInfo[self.filterKW])
        self.filter = analysisInfo[self.filterKW]
        self.top = float(analysisInfo[self.topKW])
