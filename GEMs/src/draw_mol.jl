using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using MolecularGraph

# VISUALIZE MOLECULE USING MOLFILE

id = "C00010"
graph = MolecularGraph.sdftomol("../data/molfiles/$id.mol")

svgstring = drawsvg(graph)

display("image/svg+xml",svgstring)

open("../figures/$id.svg", "w") do file
    write(file, svgstring)
end