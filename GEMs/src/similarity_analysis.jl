using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using CSV, DataFrames, CairoMakie, CategoricalArrays

aa = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/amino_acids.csv", DataFrame)
# order aa based on this symbol order: A, V, I, L, M, F, Y, W, C, U, G, P, S, T, N, Q, R, H, K, D, E
custom_order = ["G", "A", "V", "I", "L", "M", "P", "C",
                "F", "Y", "W", "S", "T", "N", "Q", 
                "D", "E", "R", "H", "K"]

aa.symbol = CategoricalArray(aa.symbol; levels=custom_order, ordered=true)
sort!(aa, :symbol)

combined_values = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/(combined_aa)_MA.tsv", DataFrame)
aa_values = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/(aa)_MA.tsv", DataFrame)

distance_matrix = fill(NaN, size(aa)[1], size(aa)[1])

for i in 1:size(aa)[1]
    for j in i:size(aa)[1]

        sym1 = aa[i, :symbol] |> string
        sym2 = aa[j, :symbol] |> string

        id1 = aa[i, :id]
        id2 = aa[j, :id]

        if sym1*"_"*sym2 in combined_values[!, :cpd]
            comb_ma = combined_values[combined_values[!, :cpd] .== sym1*"_"*sym2, :ma][1]
        else
            comb_ma = combined_values[combined_values[!, :cpd] .== sym2*"_"*sym1, :ma][1]
        end

        ma1 = aa_values[aa_values[!, :cpd] .== id1, :ma][1]
        ma2 = aa_values[aa_values[!, :cpd] .== id2, :ma][1]

        ncd = (comb_ma - min(ma1, ma2)) / max(ma1, ma2)
        distance_matrix[i, j] = ncd
        distance_matrix[j, i] = ncd

    end
end

hm = Figure();
ax = Axis(
    hm[1, 1], 
    title = "Normalized Compression Distance Matrix",
    xticks = (1:size(aa)[1], aa.symbol .|> string),
    yticks = (1:size(aa)[1], aa.symbol .|> string)
)
hmplot = CairoMakie.heatmap!(ax, distance_matrix, colormap = :viridis)
Colorbar(hm[1, 2], hmplot, label = "NCD", width = 15)
hm