from rdkit import Chem
from rdkit.Chem import rdmolfiles
import pandas as pd
import os

os.chdir("C:\\Users\\jeppe\\OneDrive\\Documenten\\Bioinformatics\\Tweede master\\Master Thesis\\Assembly_Theory\\GEMs")

aa = pd.read_csv("data/amino_acids.csv")

os.chdir("bin/assembly_go/molfiles")

for i in range(len(aa)):
    for j in range(i, len(aa)):
        
        id1 = aa["id"][i]
        id2 = aa["id"][j]
        sym1 = aa["symbol"][i]
        sym2 = aa["symbol"][j]

        mol1 = Chem.MolFromMolFile("{}.mol".format(id1), sanitize=False)
        mol2 = Chem.MolFromMolFile("{}.mol".format(id2), sanitize=False)

        combined = Chem.CombineMols(mol1, mol2)
        combined_block = Chem.MolToMolBlock(combined)

        with open("combined_aa/{}_{}.mol".format(sym1, sym2), "w") as f:
            f.write(combined_block)
