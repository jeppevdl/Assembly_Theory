using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using DataFrames, CSV, CairoMakie

df = CSV.read("../data/pathway_complexities/complete_descriptors.csv", DataFrame)

CairoMakie.scatter(df.mass, df.ma)