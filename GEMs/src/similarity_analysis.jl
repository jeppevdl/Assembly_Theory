using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using CSV, DataFrames, CairoMakie, CategoricalArrays

aa = CSV.read("../data/amino_acids.csv", DataFrame)

custom_order = ["G", "A", "V", "L", "M", "I", # nonpolar aliphatic
                "S", "T", "C", "P", "N", "Q", # polar uncharged
                "D", "E", # negatively charged
                "K", "R", "H", # positively charged
                "F", "Y", "W"] # aromatic

group_labels = [
    "nonpolar aliphatic", "nonpolar aliphatic", "nonpolar aliphatic",
    "nonpolar aliphatic", "nonpolar aliphatic", "nonpolar aliphatic",
    "polar uncharged", "polar uncharged", "polar uncharged", "polar uncharged", "polar uncharged", "polar uncharged",
    "negatively charged", "negatively charged",
    "positively charged", "positively charged", "positively charged",
    "aromatic", "aromatic", "aromatic"
    ]


aa.symbol = CategoricalArray(aa.symbol; levels=custom_order, ordered=true)
sort!(aa, :symbol)
symbol2group = Dict(custom_order .=> group_labels)
aa.group = [symbol2group[s] for s in aa.symbol]

combined_values = CSV.read("../data/bash_MA_output/(combined_aa)_MA.tsv", DataFrame)
aa_values = CSV.read("../data/bash_MA_output/(aa)_MA.tsv", DataFrame)

distance_matrix = fill(NaN, size(aa)[1], size(aa)[1])

# method can be "JAO", "NCD", or "Absolute difference"
method = "NCD"

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

        if method == "JAO"
            jao = 1 - ((ma1 + ma2) / comb_ma - 1)
            distance_matrix[i, j] = jao
            distance_matrix[j, i] = jao
        elseif method == "NCD"
            ncd = (comb_ma - min(ma1, ma2)) / max(ma1, ma2)
            distance_matrix[i, j] = ncd
            distance_matrix[j, i] = ncd
        elseif method == "Absolute_difference"
            absol = abs(ma1 - ma2)
            distance_matrix[i, j] = absol
            distance_matrix[j, i] = absol
        end
    end
end

CSV.write("../data/phylogeny/aa_distance_matrix_$method.csv", Tables.table(distance_matrix))

hm = Figure();
ax = Axis(
    hm[1, 1], 
    title = "Amino Acid Distance Matrix: $method",
    xticks = (1:size(aa)[1], aa.symbol .|> string),
    yticks = (1:size(aa)[1], aa.symbol .|> string)
)
hmplot = CairoMakie.heatmap!(ax, distance_matrix, colormap = :viridis)
Colorbar(hm[1, 2], hmplot, label = method, width = 15)
hm

CairoMakie.save("../figures/similarity_analysis/aa_similarity_$method.png", hm, px_per_unit=5)

D = distance_matrix
tiplabels = aa.symbol .|> String

using Clustering, NewickTreeTools

# Perform WPGMA (average linkage)
hc = hclust(D, linkage=:average)

open("../data/phylogeny/aa_$method.nw", "w") do io
    tree = NewickTreeTools.newick(hc, names=tiplabels)
    write(io, nwstr(tree) * "\n")
end


# using NeighborJoining
# njclusts = regNJ(D)
# nwstring = newickstring(njclusts, labels; labelinternalnodes=true)
# write("../data/phylogeny/aa_tree_nj.nw", nwstring)

using BasicTreePlots
aa_tree = readnw(String(read("../data/phylogeny/aa_$method.nw")))

# Build mapping dictionaries
symbol_to_name = Dict(row.symbol => row.name for row in eachrow(aa))
symbol_to_group = Dict(row.symbol => row.group for row in eachrow(aa))

layoutstyle = :cladogram

# Compute layout
f = Figure(size = (800, 600));
ax = Axis(f[1,1:3], xautolimitmargin = (0.1, 0.7), xtrimspine = true)
Label(f[0, 1:4], "Amino Acid WPGMA Tree Based on $method", fontsize=20, tellwidth=false)

leaves_xy, leaves_names = BasicTreePlots.nodepositions(aa_tree, layoutstyle=layoutstyle) |> BasicTreePlots.tipannotations

# Map group to categorical values and colors
groups = [symbol_to_group[s] for s in leaves_names]
group_ids = to_category(groups)
colors = generate_colors(groups);

# Plot tree
treeplot!(ax, aa_tree, layoutstyle=layoutstyle, tipannotationsvisible=false, openangle=deg2rad(10))

# Plot tip points
CairoMakie.scatter!(ax, leaves_xy, color=group_ids, colormap=colors)

# Add labels (use full names)
max_depth = maximum(first.(leaves_xy))
for (i, (xy, symbol)) in enumerate(zip(leaves_xy, leaves_names))
    name = symbol_to_name[symbol]
    text!(
        ax,
        max_depth,
        xy[2],
        text = symbol*": "*name,
        color = colors[group_ids[i]],
        fontsize = 13,
        align = (:left, :center),
        offset = (50.0f0, 0.0f0)
    )
end

# Add legend
unique_groups = unique(groups)
group_cats = to_category(unique_groups)
cmap = generate_colors(group_cats);

Legend(
    f[1,4],
    [MarkerElement(color=cmap[group_cats[i]], marker=:rect, markersize=20) for i in eachindex(unique_groups)],
    unique_groups,
    "Amino Acid Groups"
)

# Final touches
hideydecorations!(ax)
hidespines!(ax, :r, :l, :t)
ax.xticks = range(0, max_depth, length=3)
CairoMakie.ylims!(ax, (-1, last(argmax(last, leaves_xy)) + 1))

f

CairoMakie.save("../figures/similarity_analysis/aa_tree_$method.png", f, px_per_unit=5)