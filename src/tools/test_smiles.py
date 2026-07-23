from rdkit import Chem
from rdkit.Chem import AllChem


m1 = Chem.MolFromSmiles('C1=CC(O)=CC=C1C#N')
m2 = Chem.MolFromSmiles('C1=CC(C#N)=CC=C1O')
# m2 = Chem.MolFromSmiles('C1=CC(O)=CC=C1C#N')

print(Chem.MolToSmiles(m2))
print(Chem.MolToSmiles(m1))

# print(Chem.MolToSmiles(m1) == Chem.MolToSmiles(m2))
Morgan_Fingerprint = AllChem.GetMorganFingerprintAsBitVect(m1, radius=2, nBits=128)
MorganBinary1 = AllChem.DataStructs.BitVectToText(Morgan_Fingerprint)
print(MorganBinary1)

Morgan_Fingerprint = AllChem.GetMorganFingerprintAsBitVect(m2, radius=2, nBits=128)
MorganBinary2 = AllChem.DataStructs.BitVectToText(Morgan_Fingerprint)
print(MorganBinary1 == MorganBinary2)
