using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI, StatsPlots, 
ExactOptimalTransport, Distributions, StatsBase, NeighborJoining, Phylo, Plots, MultivariateStats

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

all_MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/complete_MAs.csv", DataFrame; header=true)
all_MAs = all_MAs[all_MAs.ma .!= "na", :]
all_MAs.ma = parse.(Int64, all_MAs.ma)

histogram(all_MAs.ma, dpi = 600, bins=75, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Frequency", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_histogram.png")
density(all_MAs.ma, dpi = 600, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Density", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_density.png")

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

organisms = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.lst")
organism_df = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.tsv", DataFrame; header=true)

@showprogress for organism in organisms
    name = organism_df.Species[organism_df.Organism .== organism][1]
    lines = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
    cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

    cpd_df = DataFrame(cpd=cpds, ma=[NaN for _ in 1:length(cpds)])

    for i in 1:length(cpd_df.cpd)
        cpd = cpd_df.cpd[i]
        if cpd in all_MAs.cpd && !isnan(all_MAs[all_MAs.cpd .== cpd, :].ma[1])
            cpd_df[i, "ma"] = all_MAs[all_MAs.cpd .== cpd, :].ma[1]
        end
    end
    
    histogram(cpd_df.ma, dpi = 600, title="MA values for $name", xlabel="MA value", ylabel="Frequency", legend=false)
    savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/organisms_MA_hist/MA_histogram_$organism.png")

    CSV.write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism.csv", cpd_df; header=true)
end

d = zeros(length(organisms), length(organisms))

for i in 1:length(organisms)
    organism = organisms[i]
    for j in i+1:length(organisms)
        organism2 = organisms[j]
        if organism != organism2
            cpd_df_1 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism.csv", DataFrame; header=true)
            cpd_df_2 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism2.csv", DataFrame; header=true)

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

hm = heatmap(d, dpi=600, size=(600, 500), title="Wasserstein distance between organisms", xlabel="Organisms", ylabel="Organisms")
savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/MA_heatmap.png")
CSV.write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/MA_distance_matrix.csv", Tables.table(d))

d_scaled = (d .- mean(d)) ./ std(d)

mds = fit(MDS, d; distances=true, maxoutdim=2)
Y = predict(mds)

animals = Y[:, organism_df.Kingdom .== "Animals"]
plants = Y[:, organism_df.Kingdom .== "Plants"]
fungi = Y[:, organism_df.Kingdom .== "Fungi"]
bacteria = Y[:, organism_df.Kingdom .== "Bacteria"]
protists = Y[:, organism_df.Kingdom .== "Protists"]
archaea = Y[:, organism_df.Kingdom .== "Archaea"]

mdsplot = scatter(animals[1,:], animals[2,:], marker=:circle,linewidth=0, label = "Animals", dpi = 600);
scatter!(plants[1,:], plants[2,:], marker=:circle, linewidth=0, label = "Plants");
scatter!(fungi[1,:], fungi[2,:], marker=:circle, linewidth=0, label = "Fungi");
scatter!(bacteria[1,:], bacteria[2,:], marker=:circle, linewidth=0, label = "Bacteria");
scatter!(protists[1,:], protists[2,:], marker=:circle, linewidth=0, label = "Protists");
scatter!(archaea[1,:], archaea[2,:], marker=:circle, linewidth=0, label = "Archaea");
display(mdsplot)

njclusts = regNJ(d)
# njclustsfast = fastNJ(d)
labels = organisms
# nwstring = newickstring(njclusts, labels)
nwstring = newickstring(njclusts, labels; labelinternalnodes=true)

# write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk", nwstring)
Phylo.upgma(d, labels)

matree = open(parsenewick, Phylo.path("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk"))
default(linecolor = :black, size = (400, 400))
plot(matree, size = (800, 1200), showtips = true)


using Clustering

function hclust_to_newick(tree::Clustering.Hclust, labels::Vector{String})
    n = length(labels)

    function recurse(node::Int)
        if node < 0
            return labels[-node]
        else
            l = tree.merge[node, 1]
            r = tree.merge[node, 2]
            return "(" * recurse(l) * "," * recurse(r) * ")"
        end
    end

    # There are n-1 merges, so the final node is at index n-1
    return recurse(n - 1) * ";"
end

ctree = hclust(d)

newick_str = hclust_to_newick(ctree, labels)