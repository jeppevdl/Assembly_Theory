using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end
Pkg.activate(".")

using HTTP, JSON, KEGGAPI, MolecularGraph, ProgressMeter

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

println(PROGRAM_FILE * "$(@__FILE__) - v0.0")

# Get unique paths, unique compounds and relationshipo between both
@info "Retrieving pathways and it's compounds from KEGG..."
path2cpd, allpath, allcpd = begin
	pathways, compounds = KEGGAPI.link("cpd", "pathway").data
	rels = Dict{String,Vector{String}}(path => [] for path in pathways)
	for i in eachindex(pathways)
		push!(rels[pathways[i]], compounds[i])
	end
	rels, sort(unique(pathways)), sort(unique(compounds))
end
@info "Found $(length(allcpd)) compounds associated to $(length(allpath)) pathways"

# Get MOL file for compounds and convert to InCHI
# NOTE: DO NOT PARALLELIZE THIS STEP, IT WILL ERROR OUT AND LOCK YOU OUT OF USING THE KEGG API
@info "Computing InChI code for found compounds..."
allinchi = Dict{String,Union{Missing,String}}()
@showprogress for batch in collect(Iterators.partition(allcpd, 10))
	try
		sleep(0.4)
		batch_mol = KEGGAPI.kegg_get(String.(batch), "mol", 5) |> r -> split(r[2][1], "\$\$\$\$\n", keepempty=false)
		map(zip(batch,batch_mol)) do (cpd,mol)
			try
				allinchi[cpd] = inchi(string(mol))
			catch e
				# NOTE: Polymers (and some molecules) fail to convert due to them having an atom
				# of type "R" (repetition).
				@warn "Failed to convert '$cpd' to InChI."
				allinchi[cpd] = missing
			end
		end
	catch e
		@warn e
		map(batch) do cpd
			allinchi[cpd] = missing
		end
	end
end
@info "Collected InChI codes for $(length(filter(!ismissing, collect(values(allinchi))))) compounds."
@warn "Failed to collected InChI codes for $(length(filter(ismissing, collect(values(allinchi))))) compounds"

valid_inchi = [k => allinchi[k] for k in findall(!ismissing, allinchi)]
missing_inchi = [k => allinchi[k] for k in findall(ismissing, allinchi)]

open("data/inchi_valid.json", "w") do f
    JSON.print(f, valid_inchi)
end
@info "Wrote valid InChI codes: data/inchi_valid.json"

open("data/inchi_missing.json", "w") do f
    JSON.print(f, missing_inchi)
end
@info "Wrote missing InChI codes: data/inchi_missing.json"