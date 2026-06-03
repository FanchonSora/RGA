from rdkit import Chem
from rdkit.Chem import AllChem, Mol
import numpy as np

def compute_coulomb_matrix(mol, max_atoms=50, addHs=True, sorted_matrix=False):

    if addHs:
        mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=1, params=AllChem.ETKDG())

    n_atoms = mol.GetNumAtoms()
    z = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    rval = []

    for conf in mol.GetConformers():
        d = get_interatomic_distances(conf)
        m = np.zeros((max_atoms, max_atoms))
        for i in range(n_atoms):
            for j in range(n_atoms):
                if i == j:
                    m[i, j] = 0.5 * z[i] ** 2.4
                elif i < j:
                    m[i, j] = (z[i] * z[j]) / d[i, j]
                    m[j, i] = m[i, j]
                else:
                    continue

        if sorted_matrix:
            sorted_indices = np.argsort(np.linalg.norm(m, axis=1))[::-1]
            m = m[sorted_indices]

        rval.append(m)
    rval = np.asarray(rval)
    return rval

def coulomb_matrix_eig(mol, max_atoms=50):

    cmat = compute_coulomb_matrix(mol, max_atoms=max_atoms)

    w, v = np.linalg.eig(cmat)
    w_abs = np.abs(w)
    sortidx = np.argsort(w_abs[0])[::-1]

    return w[0][sortidx]

def randomize_coulomb_matrix(m, n_samples):

    rval = []
    cmat = compute_coulomb_matrix(m)
    row_norms = np.asarray([np.linalg.norm(row) for row in cmat], dtype=float)
    rng = np.random.RandomState(42)
    for i in range(n_samples):
        e = rng.normal(size=row_norms.size)
        p = np.argsort(row_norms + e)
        new = cmat[p][:, p]  # permute rows first, then columns
        rval.append(new)
    return rval


def get_interatomic_distances(conf):

    n_atoms = conf.GetNumAtoms()
    coords = [conf.GetAtomPosition(i).__idiv__(0.52917721092) for i in range(n_atoms)]
    d = np.zeros((n_atoms, n_atoms), dtype=float)
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i < j:
                d[i, j] = coords[i].Distance(coords[j])
                d[j, i] = d[i, j]
            else:
                continue
    return d