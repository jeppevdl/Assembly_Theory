using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using MolecularGraph
id = "C00001"
graph = MolecularGraph.sdftomol("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/molfiles/$id.mol")

svgstring = drawsvg(graph)

# display("image/svg+xml",svgstring)

open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/$id.svg", "w") do file
    write(file, svgstring)
end