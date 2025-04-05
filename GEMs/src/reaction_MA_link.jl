using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using JSON3, KEGGAPI, ProgressMeter, DataFrames, CSV, Plots, StatsPlots

MAs = CSV.read("../data/bash_MA_output/pathway_MA_values.tsv", DataFrame; delim="\t")

dens = density(MAs.ma, size=(800, 600), xlabel="Molecular Assembly", ylabel="Density", title="Density plot of MA values", legend=false, lw=2, alpha = 0.75, color=:blue, dpi=600)
savefig(dens, "../data/bash_MA_output/density_plot_MA.png")

hist = histogram(MAs.ma, size=(800, 600), xlabel="Molecular Assembly", ylabel="Frequency", title="Histogram of MA values", legend=false, alpha=0.75, lw=2, color=:blue)
savefig(hist, "../data/bash_MA_output/histogram_MA.png")

pathway = "map00010"
descriptors = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/pathway_descriptors.csv", DataFrame; delim=",", header=true)


reactions = KEGGAPI.link("rn", pathway).data[2]

reaction_df = DataFrame(id=String[], equation=Vector{String}[], cpd_in=Vector{String}[], cpd_out=Vector{String}[])
 
for reaction in reactions
    response = KEGGAPI.kegg_get([reaction])[2][1]
    equation = [split(line, "EQUATION")[2][5:end] for line in split(response, "\n") if startswith(line, "EQUATION")]
    ins = [compound.match for compound in eachmatch(r"C\d{5}", split(equation[1], "<=>")[1])]
    outs = [compound.match for compound in eachmatch(r"C\d{5}", split(equation[1], "<=>")[2])]
    push!(reaction_df, (id=reaction[4:end], equation=equation[1], cpd_in=ins, cpd_out=outs), promote=true)
end

for reaction in eachrow(reaction_df)
    println("Reaction ID: $(reaction.id)")
    println("Equation: $(reaction.equation)")
    println("IN:")
    for cpd in reaction.cpd_in
        ma = descriptors.ma[descriptors.id .== cpd]
        println("$cpd -> MA:$ma")
    end
    println("OUT:")
    for cpd in reaction.cpd_out
        ma = descriptors.ma[descriptors.id .== cpd]
        println("$cpd -> MA:$ma")
    end
    println("\n")
end
