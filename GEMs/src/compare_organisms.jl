using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
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
all_MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/complete_MAs.csv", DataFrame; header=true)
all_MAs = all_MAs[all_MAs.ma .!= "na", :]
all_MAs.ma = parse.(Int64, all_MAs.ma)

# Make histogram of all MA values----------------------------------------------------------------------------------------------------------------------
allhist = Figure();
allhist[1,1] = Axis(allhist, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Frequency")
hist!(allhist[1,1], all_MAs.ma, bins=75, color = :lightseagreen)
allhist;
CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_histogram_makie.png", allhist, px_per_unit=5)

# Make density plot of all MA values----------------------------------------------------------------------------------------------------------------------
alldens = Figure();
alldens[1,1] = Axis(alldens, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Density")
CairoMakie.density!(alldens[1,1], all_MAs.ma, color = (:lightseagreen, 0.5), strokewidth=2, strokecolor = :lightseagreen)
alldens;
CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_density_makie.png", alldens, px_per_unit=5)

# Study second peak in histogram----------------------------------------------------------------------------------------------------------------------
second_peak = all_MAs[(all_MAs.ma .>= 35 .&& all_MAs.ma .<= 50), :]
br_dict = Dict{String, Int}()
no_br = 0
for cpd in second_peak.cpd
    sleep(0.4)
    response = kegg_get(["cpd:"*cpd])[2][1]
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
organisms = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.lst")
organism_df = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.tsv", DataFrame; header=true)

# Link MA values to organism compounds and produce histograms----------------------------------------------------------------------------------------------------------------------
missing_ma = []
@showprogress for organism in organisms
    name = organism_df.Species[organism_df.Organism .== organism][1]
    lines = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
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
    histo[1,1] = Axis(histo, title="MA values for $name", xlabel="MA value", ylabel="Frequency")
    hist!(histo[1,1], [ma for ma in cpd_df.ma if !isnan(ma)], bins=40, color = :dodgerblue)
    histo
    CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/organisms_MA_hist/makie_MA_histogram_$organism.png", histo, px_per_unit=5)

    CSV.write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism.csv", cpd_df; header=true)
end
missings = unique(missing_ma)

# Preallocate the distance matrix----------------------------------------------------------------------------------------------------------------------
d = zeros(length(organisms), length(organisms))

# Set to desired quantiles
lower_quantile = 0.9
upper_quantile = 1

# Calculate the Wasserstein distance between the MA distributions between the quantiles of each pair of organisms
for i in 1:length(organisms)
    organism = organisms[i]
    for j in i+1:length(organisms)
        organism2 = organisms[j]
        if organism != organism2
            cpd_df_1 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism.csv", DataFrame; header=true)
            cpd_df_2 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism2.csv", DataFrame; header=true)

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
hm = Figure();
Axis(hm[1,1], title="Heatmap of MA distance matrix: MA distribution between Q$Q1 and Q$Q2", xlabel="Organisms", ylabel="Organisms")
CairoMakie.heatmap!(hm[1,1], d, colormap = :viridis)
Colorbar(hm[1, 2], label = "Wasserstein Distance", width = 15)
hm
CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/makie_MA_heatmap_Q$(Q1)_Q$(Q2).png",
hm, px_per_unit=5)

CSV.write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/MA_distance_matrix_Q$(Q1)_Q$(Q2).csv", Tables.table(d))

# Scale the distance matrix-------------------------------------------------------------------------------------------------------------
d_scaled = (d .- mean(d)) ./ std(d)

# Make an MDS plot----------------------------------------------------------------------------------------------------------------------
mds = fit(MDS, d_scaled; distances=true, maxoutdim=2)
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
CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/mds_plot_Q$(Q1)_Q$(Q2).png", mdsplot, px_per_unit=5)

# Make a UMAP plot---------------------------------------------------------------------------------------------------------------------- 
Random.seed!(666)
embedding = umap(d_scaled, 2; metric=:precomputed)
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

CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/umap_plot_Q$(Q1)_Q$(Q2).png", umapplot, px_per_unit=5)

# Make a t-SNE plot----------------------------------------------------------------------------------------------------------------------
tsne_embedding = tsne(d_scaled, 2, 0, 50000, 70; distance=true)
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
CairoMakie.save("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/tsne_plot_Q$(Q1)_Q$(Q2).png", tsneplot, px_per_unit=5)


# Tree construction----------------------------------------------------------------------------------------------------------------------
njclusts = regNJ(d)
# njclustsfast = fastNJ(d)
labels = organisms
# nwstring = newickstring(njclusts, labels)
nwstring = newickstring(njclusts, labels; labelinternalnodes=true)

write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk", nwstring)
# Phylo.upgma(d, labels)

matree = open(parsenewick, Phylo.path("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk"))
default(linecolor = :black, size = (400, 400))
Plots.plot(matree, size = (800, 1200), showtips = true)

# using Clustering

# function hclust_to_newick(tree::Clustering.Hclust, labels::Vector{String})
#     n = length(labels)

#     function recurse(node::Int)
#         if node < 0
#             return labels[-node]
#         else
#             l = tree.merge[node, 1]
#             r = tree.merge[node, 2]
#             return "(" * recurse(l) * "," * recurse(r) * ")"
#         end
#     end

#     # There are n-1 merges, so the final node is at index n-1
#     return recurse(n - 1) * ";"
# end

# ctree = hclust(d)

# newick_str = hclust_to_newick(ctree, labels)