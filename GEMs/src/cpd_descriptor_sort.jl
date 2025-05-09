using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using HTTP, JSON3, KEGGAPI, MolecularGraph, ProgressMeter, DataFrames, CSV

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

# use code below to get the pathways and compounds from KEGG in case of no prefetched molfiles 
#SKIP to next code block in other case--------------------------------------------------------------------------------------------------------

pathway = "pathway"

# Get unique paths, unique compounds and relationship between both
@info "Retrieving pathways and it's compounds from KEGG..."
path2cpd, allpath, allcpd = begin
	pathways, compounds = KEGGAPI.link("cpd", pathway).data
	rels = Dict{String,Vector{String}}(path => [] for path in pathways)
	for i in eachindex(pathways)
		push!(rels[pathways[i]], compounds[i])
	end
	rels, sort(unique(pathways)), sort(unique(compounds))
end
@info "Found $(length(allcpd)) compounds associated to $(length(allpath)) pathways"

# compounds = path2cpd["path:"*pathway]
compounds = allcpd

# get mol files for each compound
@info "Retrieving mol files for compounds..."
allmol = Dict{String,Union{Missing,String}}()
@showprogress for batch in collect(Iterators.partition(compounds, 10))
	try
		sleep(0.4)
		batch_mol = KEGGAPI.kegg_get(String.(batch), "mol", 5) |> r -> split(r[2][1], "\$\$\$\$\n", keepempty=false)
		map(zip(batch,batch_mol)) do (cpd,mol)
			allmol[cpd] = mol
		end
	catch e
		@warn e 
	end
end
@info "Found $(length(filter(x -> !ismissing(x), allmol))) out of $(length(allcpd)) compounds with mol files"

@showprogress for (cpd, mol) in allmol
    if mol !== missing
		id = split(cpd, ":")[2]
        if !isfile("../data/molfiles/$id.mol")
			open("../data/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
		end
		if !isfile("../bin/assembly_go/molfiles/$id.mol")
			open("../bin/assembly_go/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
		end
    end
end
#---------------------------------------------------------------------------------------------------------------------------------

# use this code block below if mol files are already downloaded and MA values are precalculated

all_ma = CSV.read("../data/bash_MA_output/complete_MAs.csv", DataFrame)
complexities = DataFrame(id=String[], n_atoms=Dict{Symbol,Int}[], hybr=Vector{Int}[], n_rings=Vector{Int}[], ma=Vector{Int64}[], mass = Vector{Float64}[])
@info "Calculating complexity of compounds..."
@showprogress for row in eachrow(all_ma)
	cpd = row.cpd
	ma = row.ma
	if ma != "na"
		try
			# println(id)
			graph = MolecularGraph.sdftomol("../data/molfiles/$cpd.mol")
		
			n_atoms = MolecularGraph.atom_counter(graph)
			hybr = MolecularGraph.hybridization(graph)
			n_rings = MolecularGraph.ring_count(graph)
			mass = MolecularGraph.exact_mass(graph)

			push!(complexities, (cpd, n_atoms, hybr, n_rings, ma, mass), promote=true)

		catch e
			@warn e
			println("Failed to calculate complexity for $cpd")
		end
	end
end
@info "Calculated complexity values for $(nrow(complexities)) compounds out of $(length(all_ma[all_ma.ma .!= "na", :]))"

#---------------------------------------------------------------------------------------------------------------------------------

function count_rings(complexities::DataFrame)
	unique_rings = unique.(filter.(x -> x > 0, complexities.n_rings))
    return length.(unique_rings)
end

complexities.counted_rings = count_rings(complexities)

function count_atoms(complexities::DataFrame)
	counted_atoms = []
	for compound in eachrow(complexities)
		sum = 0
		for (key, value) in compound.n_atoms
			if key != :H
				sum += value
			end
		end
		push!(counted_atoms, sum)
	end
	return counted_atoms
end

complexities.counted_atoms = count_atoms(complexities)

function count_hybr(complexities::DataFrame)
	sp = count.(x -> x == :sp, complexities.hybr)
	sp2 = count.(x -> x == :sp2, complexities.hybr)
	sp3 = count.(x -> x == :sp3, complexities.hybr)
	res = [sp, sp2, sp3]
	return res
end

complexities.sp = count_hybr(complexities)[1]
complexities.sp2 = count_hybr(complexities)[2]
complexities.sp3 = count_hybr(complexities)[3]

complexities.sp_perc = complexities.sp ./ (complexities.counted_atoms)
complexities.sp2_perc = complexities.sp2 ./ (complexities.counted_atoms)
complexities.sp3_perc = complexities.sp3 ./ (complexities.counted_atoms)
# Sort by complexity
sort!(complexities, [:ma, :counted_atoms, :sp3, :sp2, :sp, :counted_rings], rev=true)
@info "Done sorting compounds based on complexity"

CSV.write("../data/pathway_complexities/complete_descriptors.csv", complexities; header=true)