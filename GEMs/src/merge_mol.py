from rdkit import Chem
from rdkit.Chem import rdmolfiles
import pandas as pd
import os

os.chdir("GEMs")

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

# attempt to combine all the cpds in the lookup file

# os.chdir("..\\..\\kegg-small\\data\\lookup")

# my_file = open("oyw.lookup_latest.txt", "r")
# data = my_file.read()
# data_into_list = data.replace('\n', ' ').split(" ")
# my_file.close()

# cpds = [cpd[4:] for cpd in data_into_list if cpd[0:4] == "cpd:"]

# os.chdir("..\\..\\..\\assembly_go\\molfiles")

# mol1 = Chem.MolFromMolFile("{}.mol".format(cpds[0]), sanitize=False)
# mol2 = Chem.MolFromMolFile("{}.mol".format(cpds[1]), sanitize=False)

# combined = Chem.CombineMols(mol1, mol2)
# editable = Chem.EditableMol(Chem.Mol(combined))

# for i in range(2, len(cpds)):
#     print(i)
#     mol_path = "{}.mol".format(cpds[i])
#     if os.path.isfile(mol_path):
#         print(mol_path)
#         mol = Chem.MolFromMolFile(mol_path, sanitize=False)
#         if mol is None:
#             continue
#         combined = Chem.CombineMols(editable.GetMol(), mol)
#         editable = Chem.EditableMol(Chem.Mol(combined))

# final_combined = editable.GetMol()
# Chem.SanitizeMol(final_combined)
# combined_block = Chem.MolToMolBlock(final_combined, forceV3000=False)

# with open("combined_org/{}.mol".format("oyw"), "w") as f:
#     f.write(combined_block)