using Pkg
if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")
using CSV, DataFrames, ProgressMeter, JSON3, KEGGAPI

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

MAs = CSV.read("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/bash_MA_output/pathway_MA_values.tsv", DataFrame; header=true)

organisms = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/kegg-small.lst")
failed_ids = []
calc_ma = []

for organism in organisms
	sleep(10)
	@info "Processing organism: $organism"
	lines = read_lines("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/kegg-small/data/lookup/$organism.lookup_v20250320a.txt")
	cpds = [split(line, ":")[2] for line in lines]

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
	missing_cpds = [cpd for cpd in missing_cpds if !(cpd in calc_ma)]
	missing_cpds_prefixed = "cpd:" .* missing_cpds

	# collect missing mol files
	@info "Retrieving mol files for missing compounds..."
	allmol = Dict{String,Union{Missing,String}}()
	@showprogress for batch in collect(Iterators.partition(missing_cpds_prefixed, 10))
		try
			sleep(0.4)
			batch_mol = KEGGAPI.kegg_get(String.(batch), "mol", 5) |> r -> split(r[2][1], "\$\$\$\$\n", keepempty=false)
			if length(batch_mol) != length(batch)
				@warn "Mismatched molfile count: got $(length(batch_mol)) mol files for batch of $(length(batch)) compounds"
			else
				map(zip(batch,batch_mol)) do (cpd,mol)
					allmol[cpd] = mol
				end
			end
		catch e
			@warn e 
		end
	end
	@info "Found $(length(filter(x -> !ismissing(x), allmol))) out of $(length(missing_cpds)) compounds with mol files"

	@info "Saving molfiles..."
	@showprogress for (cpd, mol) in allmol
		if mol !== missing
			id = split(cpd, ":")[2]
			open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
			open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/bin/assembly_go/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
			
		end
	end

	#generate id list for MA algorithm input for cpds with successfully retrieved mol files
	id_list = collect(keys(filter(x -> !ismissing(x), allmol)))
	id_list = [split(cpd, ":")[2] for cpd in id_list]
	open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/MA_input_$organism", "w") do f
		for id in id_list
			write(f, "$id\n")
		end
	end
	for id in id_list
		if !(id in calc_ma)
			push!(calc_ma, id)
		end
	end

	failed_mol = [cpd for cpd in missing_cpds_prefixed if !haskey(allmol, cpd)]
	for cpd in failed_mol
		if !(cpd in failed_ids)
			push!(failed_ids, cpd)
		end
	end
end

#save calc_ma to file
sort!(calc_ma)
open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/calc_ma.txt", "w") do f
	for id in calc_ma
		write(f, "$id\n")
	end
end
#save failed_ids to file
open("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/data/pathway_complexities/failed_ids.txt", "w") do f
	for id in failed_ids
		write(f, "$id\n")
	end
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