using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI, StatsPlots, ExactOptimalTransport, Distributions, StatsBase, NeighborJoining, Phylo, Plots

function read_lines(file_path::String)
    open(file_path) do file
        return readlines(file)
    end
end

all_MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/complete_MAs.csv", DataFrame; header=true)
all_MAs = all_MAs[all_MAs.ma .!= "na", :]
all_MAs.ma = parse.(Int64, all_MAs.ma)

histogram(all_MAs.ma, dpi = 600, bins=75, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Frequency", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_histogram.png")
density(all_MAs.ma, dpi = 600, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Density", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_density.png")

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
            x_probs = collect(values(x_counts)) ./ length(x_clean)

            # Clean and aggregate y
            y_clean = filter(!isnan, cpd_df_2.ma)
            y_counts = countmap(y_clean)
            y_vals = collect(keys(y_counts))
            y_probs = collect(values(y_counts)) ./ length(y_clean)

            # Create distributions
            μ = DiscreteNonParametric(x_vals, x_probs)
            ν = DiscreteNonParametric(y_vals, y_probs)

            # Compute Wasserstein-1 distance
            w = wasserstein(μ, ν; p=Val(1))

            d[i, j] = w
            d[j, i] = w
        end
    end
end

hm = heatmap(d, dpi=600, size=(600, 500), title="Wasserstein distance between organisms", xlabel="Organisms", ylabel="Organisms")
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/MA_heatmap.png")

njclusts = regNJ(d)
labels = organisms
# nwstring = newickstring(njclusts, labels)
nwstring = newickstring(njclusts, labels; labelinternalnodes=true)

write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk", nwstring)

matree = open(parsenewick, Phylo.path("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/phylogeny/tree.nwk"))
default(linecolor = :black, size = (400, 400))
plot(matree, size = (800, 1200), showtips = true)