using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")
using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI, StatsPlots, ExactOptimalTransport, Distributions

function read_lines(file_path::String)
    open(file_path) do file
        return readlines(file)
    end
end

MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/pathway_MA_values.tsv", DataFrame; header=true)

organism = "hsa"

lines = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
cpds = [split(line, ":")[2] for line in lines if startswith(line, "cpd")]

cpd_df = DataFrame(cpd=cpds, ma=[NaN for _ in 1:length(cpds)])

for i in 1:length(cpd_df.cpd)
    cpd = cpd_df.cpd[i]
    if cpd in MAs.cpd
        cpd_df[i, "ma"] = Float64(MAs[MAs.cpd .== cpd, :].ma[1])
    end
end

#number of NaN values in the ma column
count(isnan, cpd_df.ma)

missing_cpds = cpd_df.cpd[isnan.(cpd_df.ma)]
missing_cpds = "cpd:" .* missing_cpds

# collect missing mol files
@info "Retrieving mol files for missing compounds..."
allmol = Dict{String,Union{Missing,String}}()
@showprogress for batch in collect(Iterators.partition(missing_cpds, 10))
	try
		sleep(1)
		batch_mol = KEGGAPI.kegg_get(String.(batch), "mol", 5) |> r -> split(r[2][1], "\$\$\$\$\n", keepempty=false)
		map(zip(batch,batch_mol)) do (cpd,mol)
			allmol[cpd] = mol
		end
	catch e
		@warn e
		for cpd in batch
			allmol[cpd] = missing
		end
	end
end
@info "Found $(length(filter(x -> !ismissing(x), allmol))) out of $(length(missing_cpds)) missing compounds with mol files"

@info "Saving molfiles..."
@showprogress for (cpd, mol) in allmol
    if mol !== missing
		id = split(cpd, ":")[2]
        if !isfile("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/molfiles/$id.mol")
			open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
            open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/assembly_go/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
		end
    end
end

#generate id list for MA algorithm input
missing_MAs = collect(keys(filter(x -> ismissing(x), allmol)))
missing_MAs = [split(cpd, ":")[2] for cpd in missing_MAs]
open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/missing_MAs_$organism.json", "w") do f
    JSON3.write(f, missing_MAs)
end

# cpd_df_hsa
# cpd_df_eco

# x = cpd_df_hsa.ma
# y = cpd_df_gga.ma

# x_clean = filter(!isnan, x)
# y_clean = filter(!isnan, y)

# μ = fit(DiscreteNonParametric, x_clean, ones(length(x_clean)) ./ length(x_clean))
# ν = fit(DiscreteNonParametric, y_clean, ones(length(y_clean)) ./ length(y_clean))

# d = wasserstein(μ, ν; p=Val(1))
# println("Wasserstein-1 distance: ", d)