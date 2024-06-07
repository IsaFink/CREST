from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import os

SMILES=input("Enter the SMILES of the drug: ")
name=input("Enter the name of the drug: ")
SMILES=	'Cn1cnc2c1c(=O)n(C)c(=O)n2C'
drug=Chem.MolFromSmiles(SMILES)

drug_H = Chem.AddHs(drug)
# Generate 3D coordinates
AllChem.EmbedMolecule(drug_H)
drug_H_xyz=AllChem.rdmolfiles.MolToXYZBlock(drug_H)

try:
    os.mkdir(f'../{name}')
except:
    pass

#Save as file using os
with open(f'../{name}/{name}.xyz', 'w') as f:
    f.write(drug_H_xyz)

f.close()

#Add a check-in that makes sure the charge of the mmolecule is != 0.
#zero --> error
#Else, just print charge.


#explicit hydrogens
