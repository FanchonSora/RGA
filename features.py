from rdkit import Chem
from rdkit.Chem.Descriptors import NumValenceElectrons
from collections import defaultdict
from models import FeatureOptions

BOND = {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": "+"}

def parse_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol

def atom_counts(smiles) -> dict[str, int]:
    mol = parse_molecule(smiles)
    mol = Chem.AddHs(mol)

    atom_freq = defaultdict(int)
    for atom in mol.GetAtoms():
        atom_freq[atom.GetSymbol()] += 1
    
    return dict(atom_freq)

def dict_add_dict(a: dict[str, int], b: dict[str, int]) -> dict[str, int]:
    c = a.copy()
    for key, value in b.items():
        if key not in c:
            c[key] = value
        else:
            c[key] += value
    return c

def dict_multiply_scalar(a: dict[str, int], coeff: int) -> dict[str, int]:
    c = {}
    for atom_a, cnt_a in a.items():
        c[atom_a] = cnt_a * coeff
    return c

def bond_counts(smiles) -> dict[str, int]:
    bonds = {}
    mol = parse_molecule(smiles)

    for bond in mol.GetBonds():
        startAtom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        endAtom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        key = str(startAtom.GetSymbol()) + BOND[str(bond.GetBondType())] + str(endAtom.GetSymbol())
        if key not in bonds:
            bonds[key] = 1
        else:
            bonds[key] += 1
    
    return bonds

def hydrogen_bound_counts(smiles, symbol) -> dict[str, int]:
    mol = parse_molecule(smiles)
    cnt = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol().lower() == symbol.lower():
            cnt += atom.GetTotalNumHs()
    
    key = f"a{symbol}-H"
    return {key: cnt}

def electron_counts(smiles: str, opts: FeatureOptions) -> dict[str, int]:
    electron_counts = {}
    mol = parse_molecule(smiles)
    numValenceElectrons = NumValenceElectrons(mol)
    
    if opts.electron_pair:
        electron_counts["EP"] = numValenceElectrons / 2
    if opts.lone_electron:
        electron_counts["LE"] = numValenceElectrons % 2

    return electron_counts

def species_features(smiles, opts: FeatureOptions) -> dict[str, int]:
    features = {}
    if opts.electron_pair:
        features = dict_add_dict(features, electron_counts(smiles))
    if opts.hydro_bond:
        features = dict_add_dict(features, hydrogen_bound_counts(smiles, "C"))
    if opts.normal_bond:
        features = dict_add_dict(features, bond_counts(smiles))

    return features

    
