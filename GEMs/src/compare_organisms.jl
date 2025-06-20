using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI, StatsPlots, 
ExactOptimalTransport, Distributions, StatsBase, NeighborJoining, Phylo, MultivariateStats, CairoMakie, UMAP, TSne

function read_lines(file_path::String)
    open(file_path) do file
        return readlines(file)
    end
end
function KEGGAPI.kegg_get(query::Vector{String}, option::String, retries::Int)
	i = 0
	while i < retries
		try
			return KEGGAPI.kegg_get(query, option)
		catch e
			if occursin("404", string(e))
				error(e)
			elseif occursin("403", string(e))
				i += 1
				sleep(10)
			end
		end
	end
end

# Read the KEGG compound - MA value ----------------------------------------------------------------------------------------------------------------------
all_MAs = CSV.read("../data/bash_MA_output/complete_MAs.csv", DataFrame; header=true)
all_MAs = all_MAs[all_MAs.ma .!= "na", :]
all_MAs.ma = parse.(Int64, all_MAs.ma)

# Make histogram of all MA values----------------------------------------------------------------------------------------------------------------------
allhist = Figure();
allhist[1,1] = Axis(allhist, title="MA values for all gathered compounds", xlabel="MA", ylabel="Frequency", xticks=0:10:ceil(maximum(all_MAs.ma)/10)*10)
hist!(allhist[1,1], all_MAs.ma, bins=75, color = palette(:seaborn_colorblind)[1])
allhist;
CairoMakie.save("../figures/all_MA_hist/all_MA_histogram_makie.png", allhist, px_per_unit=5)

# Make density plot of all MA values----------------------------------------------------------------------------------------------------------------------
alldens = Figure();
alldens[1,1] = Axis(alldens, title="MA values for all gathered compounds", xlabel="MA", ylabel="Density")
CairoMakie.density!(alldens[1,1], all_MAs.ma, color = (:lightseagreen, 0.5), strokewidth=2, strokecolor = :lightseagreen)
alldens;
CairoMakie.save("../figures/all_MA_hist/all_MA_density_makie.png", alldens, px_per_unit=5)

# Study second peak in histogram----------------------------------------------------------------------------------------------------------------------
second_peak = all_MAs[(all_MAs.ma .>= 35 .&& all_MAs.ma .<= 50), :]
br_dict = Dict{String, Int}()
no_br = 0
n_coa = 0
@showprogress for cpd in second_peak.cpd
    sleep(0.4)
    response = kegg_get(["cpd:"*cpd])[2][1]
    name_line = match(r"NAME\s+([^\n]+)", response)
    if !isnothing(name_line) && occursin("CoA", name_line.captures[1])
        n_coa += 1
        println("CoA found for $cpd")
        continue
    end
    matches = collect(eachmatch(r"br\d{5}", response))
    br_codes = unique([m.match for m in matches])
    if br_codes == []
        println("no br codes found for $cpd")
        no_br += 1
    else
        for br_code in br_codes
            if haskey(br_dict, br_code)
                br_dict[br_code] += 1
            else
                br_dict[br_code] = 1
            end
        end
    end
end
sorted_br_dict = sort(collect(br_dict), by=x->x[2], rev=true)

# Read the list of organisms----------------------------------------------------------------------------------------------------------------------
organisms = read_lines("../bin/kegg-small/data/kegg-small.lst")
organism_df = CSV.read("../bin/kegg-small/data/kegg-small.tsv", DataFrame; header=true)

# Link MA values to organism compounds and produce histograms----------------------------------------------------------------------------------------------------------------------

missing_ma = []
@showprogress for organism in organisms
    name = organism_df.Species[organism_df.Organism .== organism][1]
    lines = read_lines("../bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
    cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

    cpd_df = DataFrame(cpd=cpds, ma=[NaN for _ in 1:length(cpds)])

    for i in 1:length(cpd_df.cpd)
        cpd = cpd_df.cpd[i]
        if cpd in all_MAs.cpd && !isnan(all_MAs[all_MAs.cpd .== cpd, :].ma[1])
            cpd_df[i, "ma"] = all_MAs[all_MAs.cpd .== cpd, :].ma[1]
        elseif !(cpd in all_MAs.cpd)
            # println("cpd $cpd not found in all_MAs")
            missing_ma = push!(missing_ma, cpd)
        end
    end
    
    histo = Figure();
    histo[1,1] = Axis(histo, title="MA values for $name", xlabel="MA", ylabel="Frequency", xticks=0:10:ceil(maximum([ma for ma in cpd_df.ma if !isnan(ma)])/10)*10,
    titlesize=20, xlabelsize=17, ylabelsize=17, xticklabelsize=14, yticklabelsize=14)
    hist!(histo[1,1], [ma for ma in cpd_df.ma if !isnan(ma)], bins=40, color = palette(:seaborn_colorblind)[2])
    histo
    CairoMakie.save("../figures/organisms_MA_hist/makie_MA_histogram_$organism.png", histo, px_per_unit=5)

    CSV.write("../data/pathway_complexities/organism_MA_values/MA_$organism.csv", cpd_df; header=true)
end
missings = unique(missing_ma)

# Make histogram of all unique MA values per kingdom----------------------------------------------------------------------------------------------------------------------
# Dictionary to collect unique compound -> MA per kingdom
kingdom_cpds = Dict{String, Set{String}}()
kingdom_ma_values = Dict{String, Vector{Float64}}()

@showprogress for organism in unique(organism_df.Organism)
    name = organism_df.Species[organism_df.Organism .== organism][1]
    kingdom = organism_df.Kingdom[organism_df.Organism .== organism][1]

    lines = read_lines("../bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
    cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

    if !haskey(kingdom_cpds, kingdom)
        kingdom_cpds[kingdom] = Set{String}()
    end
    union!(kingdom_cpds[kingdom], cpds)
end

# Now look up MA values per unique compound per kingdom
for (kingdom, cpds) in kingdom_cpds
    kingdom_ma_values[kingdom] = Float64[]
    for cpd in cpds
        if cpd in all_MAs.cpd
            ma_val = all_MAs[all_MAs.cpd .== cpd, :].ma[1]
            if !isnan(ma_val)
                push!(kingdom_ma_values[kingdom], ma_val)
                if ma_val > 100 && kingdom == "Bacteria"
                    println("$cpd: $ma_val")
                end
            end
        end
    end
end

# Plotting all histograms together
fig = Figure(size=(800,600));
ax = Axis(fig[1, 1], title="MA Value Distributions by Kingdom (Unique Compounds)", xlabel="Kingdom", ylabel="MA", width=600, xticklabelsvisible=false)

base_colors = palette(:seaborn_colorblind);
i = 1
medians = []
for (kingdom, ma_values) in sort(collect(kingdom_ma_values))
    CairoMakie.boxplot!(ax, i*ones(length(ma_values)), ma_values; color = base_colors[i], label=kingdom, width=0.5)
    push!(medians, median(ma_values))
    i += 1
end

CairoMakie.Legend(fig[1,2], ax, "Kingdoms")
fig

CairoMakie.save("../figures/organism_comparison/kingdom_MA_boxplot.png", fig, px_per_unit=5)
CSV.write("../figures/organism_comparison/kingdom_MA_medians.csv", DataFrame(kingdom=sort(collect(keys(kingdom_ma_values))), median=medians); header=true)

# Make histogram of all unique MA values per superkingdom----------------------------------------------------------------------------------------------------------------------
# Dictionary to collect unique compound -> MA per superkingdom
superkingdom_cpds = Dict{String, Set{String}}()
superkingdom_ma_values = Dict{String, Vector{Float64}}()

@showprogress for organism in unique(organism_df.Organism)
    name = organism_df.Species[organism_df.Organism .== organism][1]
    superkingdom = organism_df.Superkingdom[organism_df.Organism .== organism][1]

    lines = read_lines("../bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
    cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

    if !haskey(superkingdom_cpds, superkingdom)
        superkingdom_cpds[superkingdom] = Set{String}()
    end
    union!(superkingdom_cpds[superkingdom], cpds)
end

# Now look up MA values per unique compound per kingdom
for (superkingdom, cpds) in superkingdom_cpds
    superkingdom_ma_values[superkingdom] = Float64[]
    for cpd in cpds
        if cpd in all_MAs.cpd
            ma_val = all_MAs[all_MAs.cpd .== cpd, :].ma[1]
            if !isnan(ma_val)
                push!(superkingdom_ma_values[superkingdom], ma_val)
                if ma_val > 100 && superkingdom == "Prokaryotes"
                    println("$cpd: $ma_val")
                end
            end
        end
    end
end

# Plotting all histograms together
fig = Figure(size=(800,600));
ax = Axis(fig[1, 1], title="MA Value Distributions by Superkingdom (Unique Compounds)", xlabel="Superkingdom", ylabel="MA", width=600, xticklabelsvisible=false)

base_colors = palette(:seaborn_colorblind)[9:10];
i = 1
medians = []
for (superkingdom, ma_values) in sort(collect(superkingdom_ma_values))
    CairoMakie.boxplot!(ax, i*ones(length(ma_values)), ma_values; color = base_colors[i], label=superkingdom, width=0.5)
    push!(medians, median(ma_values))
    i += 1
end

CairoMakie.Legend(fig[1,2], ax, "Superkingdoms")
fig

CairoMakie.save("../figures/organism_comparison/superkingdom_MA_boxplot.png", fig, px_per_unit=5)
CSV.write("../figures/organism_comparison/superkingdom_MA_medians.csv", DataFrame(superkingdom=sort(collect(keys(superkingdom_ma_values))), median=medians); header=true)

# Preallocate the distance matrix----------------------------------------------------------------------------------------------------------------------
d = zeros(length(organisms), length(organisms))

# Set to desired quantiles (e.g. lower percentile 90 = Q90 => lower_quantile = 0.9)
lower_quantile = 0
upper_quantile = 1

# Calculate the Wasserstein distance between the MA distributions between the quantiles of each pair of organisms
for i in 1:length(organisms)
    organism = organisms[i]
    for j in i+1:length(organisms)
        organism2 = organisms[j]
        if organism != organism2
            cpd_df_1 = CSV.read("../data/pathway_complexities/organism_MA_values/MA_$organism.csv", DataFrame; header=true)
            cpd_df_2 = CSV.read("../data/pathway_complexities/organism_MA_values/MA_$organism2.csv", DataFrame; header=true)

            q1_1 = quantile(cpd_df_1.ma[isnan.(cpd_df_1.ma) .|> !], lower_quantile)
            q1_2 = quantile(cpd_df_1.ma[isnan.(cpd_df_1.ma) .|> !], upper_quantile)
            q2_1 = quantile(cpd_df_2.ma[isnan.(cpd_df_2.ma) .|> !], lower_quantile)
            q2_2 = quantile(cpd_df_2.ma[isnan.(cpd_df_2.ma) .|> !], upper_quantile)

            # Filter out the values outside the quantiles
            cpd_df_1 = cpd_df_1[cpd_df_1.ma .>= q1_1 .&& cpd_df_1.ma .<= q1_2, :]
            cpd_df_2 = cpd_df_2[cpd_df_2.ma .>= q2_1 .&& cpd_df_2.ma .<= q2_2, :]

            # Clean and aggregate x
            x_clean = filter(!isnan, cpd_df_1.ma)
            x_counts = countmap(x_clean)
            x_vals = collect(keys(x_counts))
            # x_probs = collect(values(x_counts)) ./ length(x_clean)

            # Clean and aggregate y
            y_clean = filter(!isnan, cpd_df_2.ma)
            y_counts = countmap(y_clean)
            y_vals = collect(keys(y_counts))
            # y_probs = collect(values(y_counts)) ./ length(y_clean)

            all_vals = union(keys(x_counts), keys(y_counts))
            x_probs = [get(x_counts, v, 0) / length(x_clean) for v in all_vals]
            y_probs = [get(y_counts, v, 0) / length(y_clean) for v in all_vals]

            μ = DiscreteNonParametric(collect(all_vals), x_probs)
            ν = DiscreteNonParametric(collect(all_vals), y_probs)

            w = wasserstein(μ, ν; p=Val(1))
            d[i, j] = d[j, i] = w
        end
    end
end

# Plot the distance matrix as a heatmap----------------------------------------------------------------------------------------------------------------------
Q1 = lower_quantile * 100 |> Int
Q2 = upper_quantile * 100 |> Int

hm = Figure(size=(800,600));
ax = Axis(hm[1,1], title="Heatmap of MA distance matrix: MA distribution between percentile $Q1 and $Q2", xlabel="Organisms", ylabel="Organisms",
yticks=(1:length(organisms), organisms), xticklabelsvisible=false, xticksvisible=false, yticklabelsize=7, titlesize=17)
heatmap_plot = CairoMakie.heatmap!(ax, d, colormap = :viridis)
Colorbar(hm[1, 2], heatmap_plot, label = "Wasserstein Distance", width = 15)
hm
CairoMakie.save("../figures/organism_comparison/makie_MA_heatmap_Q$(Q1)_Q$(Q2).png",
hm, px_per_unit=5)

CSV.write("../data/phylogeny/MA_distance_matrix_Q$(Q1)_Q$(Q2).csv", Tables.table(d))

# Scale the distance matrix-------------------------------------------------------------------------------------------------------------
# d_scaled = (d .- mean(d)) ./ std(d)
d_minmax = (d .- minimum(d)) ./ (maximum(d) - minimum(d))  # scale to [0, 1]

# Make an MDS plot----------------------------------------------------------------------------------------------------------------------
mds = fit(MDS, d_minmax; distances=true, maxoutdim=2)
Y = predict(mds)

animals = Y[:, organism_df.Kingdom .== "Animals"]
plants = Y[:, organism_df.Kingdom .== "Plants"]
fungi = Y[:, organism_df.Kingdom .== "Fungi"]
bacteria = Y[:, organism_df.Kingdom .== "Bacteria"]
protists = Y[:, organism_df.Kingdom .== "Protists"]
archaea = Y[:, organism_df.Kingdom .== "Archaea"]

mdsplot = Figure(size=(800,600));
ax = Axis(mdsplot[1, 1], title = "MDS Plot: MA distribution between Q$Q1 and Q$Q2", xlabel = "Dimension 1", ylabel = "Dimension 2", width = 600)
palet = palette(:seaborn_colorblind);
CairoMakie.scatter!(ax, animals[1, :], animals[2, :], color = palet[1], label = "Animals")
CairoMakie.scatter!(ax, plants[1, :], plants[2, :], color = palet[2], label = "Plants")
CairoMakie.scatter!(ax, fungi[1, :], fungi[2, :], color = palet[3], label = "Fungi")
CairoMakie.scatter!(ax, bacteria[1, :], bacteria[2, :], color = palet[4], label = "Bacteria")
CairoMakie.scatter!(ax, protists[1, :], protists[2, :], color = palet[5], label = "Protists")
CairoMakie.scatter!(ax, archaea[1, :], archaea[2, :], color = palet[9], label = "Archaea")
Legend(mdsplot[1, 2], ax, "Kingdoms", tellwidth = false)
mdsplot
CairoMakie.save("../figures/organism_comparison/mds_plot_Q$(Q1)_Q$(Q2).png", mdsplot, px_per_unit=5)

# Make a UMAP plot---------------------------------------------------------------------------------------------------------------------- 
using Random
Random.seed!(666)
embedding = umap(d_minmax, 2; metric=:precomputed)
animals = embedding[:, organism_df.Kingdom .== "Animals"]
plants = embedding[:, organism_df.Kingdom .== "Plants"]
fungi = embedding[:, organism_df.Kingdom .== "Fungi"]
bacteria = embedding[:, organism_df.Kingdom .== "Bacteria"]
protists = embedding[:, organism_df.Kingdom .== "Protists"]
archaea = embedding[:, organism_df.Kingdom .== "Archaea"]

umapplot = Figure(size=(800,600));
ax = Axis(umapplot[1, 1], title = "UMAP Plot: MA distribution between Q$Q1 and Q$Q2", xlabel = "Dimension 1", ylabel = "Dimension 2", width = 600)
palet = palette(:seaborn_colorblind);
CairoMakie.scatter!(ax, animals[1, :], animals[2, :], color = palet[1], label = "Animals")
CairoMakie.scatter!(ax, plants[1, :], plants[2, :], color = palet[2], label = "Plants")
CairoMakie.scatter!(ax, fungi[1, :], fungi[2, :], color = palet[3], label = "Fungi")
CairoMakie.scatter!(ax, bacteria[1, :], bacteria[2, :], color = palet[4], label = "Bacteria")
CairoMakie.scatter!(ax, protists[1, :], protists[2, :], color = palet[5], label = "Protists")
CairoMakie.scatter!(ax, archaea[1, :], archaea[2, :], color = palet[9], label = "Archaea")
Legend(umapplot[1, 2], ax, "Kingdoms", tellwidth = false)
umapplot

CairoMakie.save("../figures/organism_comparison/umap_plot_Q$(Q1)_Q$(Q2).png", umapplot, px_per_unit=5)

# Make a t-SNE plot----------------------------------------------------------------------------------------------------------------------
tsne_embedding = tsne(d_minmax, 2, 0, 50000, 80; distance=true)
tsne_embedding = tsne_embedding'
animals = tsne_embedding[:, organism_df.Kingdom .== "Animals"]
plants = tsne_embedding[:, organism_df.Kingdom .== "Plants"]
fungi = tsne_embedding[:, organism_df.Kingdom .== "Fungi"]
bacteria = tsne_embedding[:, organism_df.Kingdom .== "Bacteria"]
protists = tsne_embedding[:, organism_df.Kingdom .== "Protists"]
archaea = tsne_embedding[:, organism_df.Kingdom .== "Archaea"]

tsneplot = Figure(size=(800,600));
ax = Axis(tsneplot[1, 1], title = "t-SNE Plot: MA distribution between Q$Q1 and Q$Q2", xlabel = "Dimension 1", ylabel = "Dimension 2", width = 600)
palet = palette(:seaborn_colorblind);
CairoMakie.scatter!(ax, animals[1, :], animals[2, :], color = palet[1], label = "Animals")
CairoMakie.scatter!(ax, plants[1, :], plants[2, :], color = palet[2], label = "Plants")
CairoMakie.scatter!(ax, fungi[1, :], fungi[2, :], color = palet[3], label = "Fungi")
CairoMakie.scatter!(ax, bacteria[1, :], bacteria[2, :], color = palet[4], label = "Bacteria")
CairoMakie.scatter!(ax, protists[1, :], protists[2, :], color = palet[5], label = "Protists")
CairoMakie.scatter!(ax, archaea[1, :], archaea[2, :], color = palet[9], label = "Archaea")
Legend(tsneplot[1, 2], ax, "Kingdoms", tellwidth = false)
tsneplot
CairoMakie.save("../figures/organism_comparison/tsne_plot_Q$(Q1)_Q$(Q2).png", tsneplot, px_per_unit=5)

# Minimum Spanning Tree----------------------------------------------------------------------------------------------------------------------
using Graphs
using GraphMakie, NetworkLayout

A = d_minmax
G = Graph(A)

M = SimpleGraphFromIterator(prim_mst(G, A))

kingdoms = ["Animals", "Plants", "Fungi", "Bacteria", "Protists", "Archaea"]
palet = palette(:seaborn_colorblind);
colors = Dict("Animals" => palet[1], "Plants" => palet[2], "Fungi" => palet[3],
              "Bacteria" => palet[4], "Protists" => palet[5], "Archaea" => palet[9])

node_color = [colors[k] for k in organism_df.Kingdom];

f = Figure(size=(800, 600));
ax = Axis(f[1, 1], title = "MST: MA distribution between Q$Q1 and Q$Q2", width = 600)
graphplot!(ax, M; node_color=node_color, node_size=10, layout=Stress(), nlabels = organisms, nlabels_fontsize=10)

legend_handles = [CairoMakie.scatter!(ax, [NaN], [NaN]; color=colors[k], label=k) for k in kingdoms]
Legend(f[1, 2], ax, "Kingdoms"; tellwidth=false)

f

CairoMakie.save("../figures/organism_comparison/mst_plot_Q$(Q1)_Q$(Q2).png", f, px_per_unit=5)

# Tree construction----------------------------------------------------------------------------------------------------------------------
tiplabels = organisms

using Clustering, NewickTreeTools

# Perform UPGMA (average linkage)
hc = hclust(d, linkage=:average)

open("../data/phylogeny/tree_Q$(Q1)_Q$(Q2).nw", "w") do io
    tree = NewickTreeTools.newick(hc, names=tiplabels)
    write(io, nwstr(tree) * "\n")
end

using BasicTreePlots

# Read tree and taxonomy
org_tree = readnw(String(read("../data/phylogeny/tree_Q$(Q1)_Q$(Q2).nw")))
taxonomy = organism_df
org2tax = map(eachrow(taxonomy)) do r
    (r.Organism => Dict(
        "Superkingdom" => r.Superkingdom,
        "Kingdom" => r.Kingdom,
        "Genus" => r.Genus,
        "Phylum" => r.Phylum,
        "Species" => r.Species,
    ))
end |> Dict

# Utility functions
specie2tax(s, tax = :Kingdom) = Dict(map(eachrow(taxonomy[!, [:Organism, :Superkingdom, :Kingdom, :Phylum, :Genus]])) do r
    r.Organism => r[tax]
end)[s]

function to_category(values::Vector{String})
    labels = sort(unique(values))
    mapper = Dict(zip(labels, 1:length(labels)))
    return map(v -> mapper[v], values), labels
end

function generate_colors(n)
    distinguishable_colors(n, RGB(1,1,1); dropseed=true, lchoices=range(0, stop=95, length=15))
end

# Plot setup
layoutstyle = :cladogram
f1 = Figure(size=(800,1200))
f1a1 = Axis(f1[1,1:3], xautolimitmargin = (0.1, 0.7), xtrimspine = true,
    title = "Hierarchical clustering tree: MA distributions between Q$Q1 and Q$Q2")

# Tree layout and tip names
leaves_xy, leaves_names = BasicTreePlots.nodepositions(org_tree, layoutstyle=layoutstyle) |> BasicTreePlots.tipannotations

# Assign Kingdom-based colors
kingdoms = String.(specie2tax.(leaves_names, :Kingdom))
kingdom_cat, kingdom_labels = to_category(kingdoms)
kingdom_colors = generate_colors(length(kingdom_labels))

# Assign Superkingdom-based colors (ensuring they do not overlap with Kingdom colors)
superkingdoms = String.(specie2tax.(leaves_names, :Superkingdom))
superkingdom_cat, superkingdom_labels = to_category(superkingdoms)
superkingdom_colors = generate_colors(length(superkingdom_labels) + length(kingdom_labels))[end-length(superkingdom_labels)+1:end]

# Add tree
treeplot!(
    f1a1,
    org_tree,
    layoutstyle=layoutstyle,
    tipannotationsvisible=false,
    openangle=deg2rad(10)
)

# Plot colored dots by Kingdom (main color)
CairoMakie.scatter!(
    f1a1,
    leaves_xy,
    color = kingdom_cat,
    colormap = kingdom_colors
)

# Add labels with Kingdom-based colors
max_depth = first(argmax(first, leaves_xy))
for (i, (xy, name)) in enumerate(zip(leaves_xy, leaves_names))
    text!(
        max_depth,
        last(xy);
        text=name*" * "*org2tax[name]["Species"],
        color=kingdom_colors[kingdom_cat[i]],
        fontsize=9,
        align = (:left, :center),
        offset = (50.0f0, 0.0f0)
    )
end

# Add Superkingdom and Kingdom bars with correct coloring
for (i, (tax, cats, labels, colormap)) in enumerate([
        (:Superkingdom, superkingdom_cat, superkingdom_labels, superkingdom_colors),
        (:Kingdom, kingdom_cat, kingdom_labels, kingdom_colors),
    ])
    coords = [(max_depth * (1 + i / 30), y) for (_, y) in leaves_xy]
    CairoMakie.scatter!(
        f1a1,
        coords,
        color = cats,
        colormap = colormap,
        marker = :rect,
        markersize = 14
    )

    Legend(f1[2,i],
        [
            MarkerElement(
                color = colormap[j],
                marker = :rect,
                markersize = 25
            ) for j in 1:length(labels)
        ],
        labels,
        "KEGG Tax\nLevel $i",
        nbanks = i,
        rowgap = 1,
        unique = true,
        merge = true
    )
end

# Final plot cleanup
hideydecorations!(f1a1)
hidespines!(f1a1, :r, :l, :t)
f1a1.xticks = range(0, trunc(max_depth;sigdigits=1), 3)
CairoMakie.ylims!(f1a1, (-1, last(argmax(last, leaves_xy)) + 1))
rowsize!(f1.layout, 2, Relative(1/5))
f1

CairoMakie.save("../figures/organism_comparison/tree_Q$(Q1)_Q$(Q2).png", f1; px_per_unit=5)
