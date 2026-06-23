from models import *
from balancing import *
from features import atom_counts, species_features, dict_add_dict
from reaction_io import format_reaction_list
from itertools import combinations

def gen_2_to_1(target, species_list: dict[str, Species], opts: FeatureOptions, atoms_list: dict, features_list: dict) -> list[Reaction]:
    reactions = []

    target_atoms = atoms_list[target]
    target_features = features_list[target]

    candidates = [smiles for smiles in species_list.keys() if smiles != target]

    for reactant in candidates:
        reactant_atoms = dict_add_dict(target_atoms, atoms_list[reactant])

        for product in candidates:
            if product == reactant:
                continue

            product_atoms = atoms_list[product]

            if reactant_atoms != product_atoms:
                continue

            reactant_features = dict_add_dict(target_features, features_list[reactant])
            product_features = features_list[product]

            if reactant_features == product_features:
                reactions.append(Reaction(
                    reactants=((1, target), (1, reactant)),
                    products=((1, product),)
                ))

    return reactions


def gen_2_to_2(target, species_list: dict[str, Species], opts: FeatureOptions, atoms_list: dict, features_list: dict) -> list[Reaction]:
    reactions = []

    target_atoms = atoms_list[target]
    target_features = features_list[target]

    candidates = [smiles for smiles in species_list.keys() if smiles != target]

    for product1, product2 in combinations(candidates, 2):
        product_atoms = dict_add_dict(atoms_list[product1], atoms_list[product2])

        for reactant in candidates:
            if reactant == product1 or reactant == product2:
                continue

            reactant_atoms = dict_add_dict(target_atoms, atoms_list[reactant])

            if reactant_atoms != product_atoms:
                continue

            reactant_features = dict_add_dict(target_features, features_list[reactant])
            product_features = dict_add_dict(features_list[product1], features_list[product2])

            if reactant_features == product_features:
                reactions.append(Reaction(
                    reactants=((1,  target), (1, reactant)),
                    products=((1, product1), (1, product2))
                ))

    return reactions

def gen_3_to_1(target, species_list: dict[str, Species], opts: FeatureOptions, atoms_list: dict, features_list: dict) -> list[Reaction]:
    reactions = []

    target_atoms = atoms_list[target]
    target_features = features_list[target]
    
    candidates = [smiles for smiles in species_list.keys() if smiles != target]

    for reactant1, reactant2 in combinations(candidates, 2):
        reactant_atoms = dict_add_dict(atoms_list[reactant1], atoms_list[reactant2])
        reactant_atoms = dict_add_dict(reactant_atoms, target_atoms)

        for product in candidates:
            if product == reactant1 or product == reactant2:
                continue

            product_atoms = atoms_list[product]

            if reactant_atoms != product_atoms:
                continue

            reactant_features = dict_add_dict(features_list[reactant1], features_list[reactant2])
            reactant_features = dict_add_dict(reactant_features, target_features)
            product_features = features_list[product]

            if reactant_features == product_features:
                reactions.append(Reaction(
                    reactants=((1, target), (1, reactant1), (1, reactant2)),
                    products=((1, product),)
                ))

    return reactions

def gen_reactions(target, species_list: dict[str, Species], opts: FeatureOptions, max_species: int = 4) -> str:
    atoms_list = {smiles: atom_counts(smiles) for smiles in species_list}
    features_list = {smiles: species_features(smiles, opts) for smiles in species_list}
    reactions = []

    if target not in atoms_list:
        atoms_list[target] = atom_counts(target)
    if target not in features_list:
        features_list[target] = species_features(target, opts)

    reactions.extend(gen_2_to_1(target, species_list, opts, atoms_list, features_list))

    if (max_species >= 4):
        reactions.extend(gen_2_to_2(target, species_list, opts, atoms_list, features_list))
        reactions.extend(gen_3_to_1(target, species_list, opts, atoms_list, features_list))


    return format_reaction_list(reactions)
