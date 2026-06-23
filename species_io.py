from models import Species
from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict

def load_species(path) -> dict[str, Species]:
    species_dict = {}

    with open(path, "r") as file:
        lines = file.readlines()
        for line in lines:
            smiles, experimental_hof, uncertainty, calculated_hof = line.split()
            species = Species(smiles, float(experimental_hof), float(uncertainty), float(calculated_hof))

            if smiles in species_dict:
                raise ValueError(f"Duplicate species: {smiles}")
            
            species_dict[smiles] = species
    return species_dict

def parse_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        raise ValueError(f"Invalid smiles: {smiles}")
    return mol
