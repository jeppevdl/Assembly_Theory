from rdkit import Chem
from rdkit.Chem import rdmolfiles
import os
os.chdir("GEMs/bin/assembly_go/molfiles")

mol1 = Chem.MolFromMolFile("C00078.mol", sanitize=False)
mol2 = Chem.MolFromMolFile("C00064.mol", sanitize=False)

combined = Chem.CombineMols(mol1, mol2)
combined_block = Chem.MolToMolBlock(combined)

print(combined_block)
with open("acombined.mol", "w") as f:
    f.write(combined_block)
