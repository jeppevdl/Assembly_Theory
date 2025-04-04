using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using JSON3, KEGGAPI, ProgressMeter, DataFrames, CSV, Plots, StatsPlots

MAs = CSV.read("../data/bash_MA_output/intermediate_MA_pathway_no_MA.tsv", DataFrame; delim="\t")
MAs2 = CSV.read("../data/bash_MA_output/intermediate_MA_v2_pathway_no_MA.tsv", DataFrame, delim="\t")

MAs = append!(MAs, MAs2)

#make density plot of the MAs
density(MAs.ma, xlabel="Molecular assembly", ylabel="Density", title="Density plot of molecular complexities", legend=false, alpha=0.5, lw=2, color=:blue)


pathway = "map00010"
json_data = JSON3.read(open("data/complexities_$pathway.json", "r"))
columns = json_data[:columns] 
complexities = DataFrame(columns, json_data[:colindex][:names])

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
        ma = complexities.ma[complexities.id .== cpd]
        println("$cpd -> MA:$ma")
    end
    println("OUT:")
    for cpd in reaction.cpd_out
        ma = complexities.ma[complexities.id .== cpd]
        println("$cpd -> MA:$ma")
    end
    println("\n")
end
