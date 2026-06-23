import argparse
from species_io import load_species
from models import *
from generation import *
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--species_file', type=str, default='CBSQB3_2019.10.25.txt')
    opt = parser.parse_args()
    
    species_dict = load_species(opt.species_file)
    target = "C1=CC=C2C=CC=CC2=C1"
    opts = FeatureOptions(electron_pair = True, hydro_bond = True, normal_bond = False)
    print(gen_reactions(target, species_dict, opts))
    
