using MolecularGraph

cd("../bin/assembly_go/molfiles")

function join_molfiles(file1::String, file2::String, output::String)
    # Read the two molecules
    mol1 = sdftomol(file1)
    mol2 = sdftomol(file2)

    # Get atoms and bonds from both molecules
    atoms1 = atom_symbol(mol1)
    bonds1 = 

    atoms2 = mol2.atoms
    bonds2 = mol2.bonds

    # Shift atom indices in mol2 so they follow mol1's atoms
    offset = length(atoms1)

    # Combine atoms
    combined_atoms = vcat(atoms1, atoms2)

    # Adjust bonds in mol2 to account for new indices
    shifted_bonds2 = map(b -> MolecularGraph.Bond(b.src + offset, b.dst + offset, b.order), bonds2)

    # Combine bonds
    combined_bonds = vcat(bonds1, shifted_bonds2)

    # Create a new molecule
    combined_mol = MolecularGraph.GraphMol(combined_atoms, combined_bonds)

    # Write the combined molecule to a file
    writemol(combined_mol, output)

    println("✅ Combined molecule written to: $output")
end

join_molfiles("C00001.mol", "C00002.mol", "test/C00001_C00002.mol")