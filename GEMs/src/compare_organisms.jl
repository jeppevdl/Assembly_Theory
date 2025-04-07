using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI, StatsPlots, ExactOptimalTransport

function read_lines(file_path::String)
    open(file_path) do file
        return readlines(file)
    end
end

all_MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/all_MAs.csv", DataFrame; header=true)
all_MAs = all_MAs[all_MAs.ma .!= "na", :]
all_MAs.ma = parse.(Int64, all_MAs.ma)

# histogram(all_MAs.ma, dpi = 600, bins=75, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Frequency", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_histogram.png")
# density(all_MAs.ma, dpi = 600, title="MA values for all gathered compounds", xlabel="MA value", ylabel="Density", legend=false)
# savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/all_MA_hist/all_MA_density.png")

organisms = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.lst")

@showprogress for organism in organisms
    lines = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
    cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

    cpd_df = DataFrame(cpd=cpds, ma=[NaN for _ in 1:length(cpds)])

    for i in 1:length(cpd_df.cpd)
        cpd = cpd_df.cpd[i]
        if cpd in all_MAs.cpd && !isnan(all_MAs[all_MAs.cpd .== cpd, :].ma[1])
            cpd_df[i, "ma"] = all_MAs[all_MAs.cpd .== cpd, :].ma[1]
        end
    end
    
    histogram(cpd_df.ma, dpi = 600, title="MA values for $organism", xlabel="MA value", ylabel="Frequency", legend=false, xlims=(-1, 150), ylims=(0, 550))
    savefig("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/figures/organisms_MA_hist/MA_histogram_$organism.png")

    CSV.write("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/organism_MA_values/MA_$organism.csv", cpd_df; header=true)
end

# preallocate matrix for wasserstein distances
distances = zeros(length(organisms), length(organisms))
for i in 1:length(organisms)
    distances[i, i] = NaN
end
for i in 1:length(organisms)
    organism = organisms[i]
    for j in i+1:length(organisms)
        organism2 = organisms[j]
        if organism != organism2
            cpd_df_1 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/MA_$organism.csv", DataFrame; header=true)
            cpd_df_2 = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/MA_$organism2.csv", DataFrame; header=true)

            x = cpd_df_1.ma
            y = cpd_df_2.ma

            x_clean = filter(!isnan, x)
            y_clean = filter(!isnan, y)

            μ = fit(DiscreteNonParametric, x_clean, ones(length(x_clean)) ./ length(x_clean))
            ν = fit(DiscreteNonParametric, y_clean, ones(length(y_clean)) ./ length(y_clean))

            d = wasserstein(μ, ν; p=Val(1))
            println("Wasserstein-1 distance: ", d)
            distances[i, j] = d
            distances[j, i] = d
        end
    end
end

organism_df = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.tsv", DataFrame; header=true, delim="\t")
