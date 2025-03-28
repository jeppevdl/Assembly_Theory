using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end
Pkg.activate(".")

using HTTP, JSON3, KEGGAPI, MolecularGraph, ProgressMeter, DataFrames

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

pathway = "map01230"

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

compounds = path2cpd["path:"*pathway]

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
        map(batch) do cpd
			allmol[cpd] = missing
        end
	end
end

# Calculate complexity of compounds
complexities = DataFrame(id=String[], n_atoms=Dict{Symbol,Int}[], hybr=Vector{Int}[], n_rings=Vector{Int}[])
@info "Calculating complexity of compounds..."
@showprogress for (cpd, mol) in allmol
    if mol !== missing
		id = split(cpd, ":")[2]
        if !isfile("molfiles/$id.mol")
			open("molfiles/$id.mol", "w") do f
				write(f, mol)
			end
		end
		try 
			println(id)
			graph = MolecularGraph.sdftomol("molfiles/$id.mol")
			
			n_atoms = MolecularGraph.atom_counter(graph)
			hybr = MolecularGraph.hybridization(graph)
			n_rings = MolecularGraph.ring_count(graph)

			push!(complexities, (id, n_atoms, hybr, n_rings), promote=true)

		catch e
			@warn e
		end
    end
end
@info "Calculated complexity values for $(nrow(complexities)) compounds out of $(length(allmol))"

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

# Sort by complexity
sort!(complexities, [:counted_atoms, :sp3, :sp2, :sp, :counted_rings], rev=true)
@info "Done sorting compounds based on complexity"

ma_values = JSON3.read("data/MA_values.json", Dict)

# Get MA values for each compound
@info "Retrieving MA values for compounds..."
complexities.ma = Vector{Union{Missing,Float64}}(missing, nrow(complexities))
@showprogress for i in 1:nrow(complexities)
	cpd = complexities.id[i]
	if haskey(ma_values, cpd)
		complexities.ma[i] = ma_values[cpd]["MA"]
	end
end

n = count(x -> !ismissing(x), complexities.ma)
@info "Found MA values for $n compounds"

@info "Calculating MA values using molecular assembly algorithm..."
@showprogress for i in nrow(complexities):-1:1
	if ismissing(complexities.ma[i])
		println(complexities.id[i])
		output = readchomp(`powershell -Command "complexity molfiles/$(complexities.id[i]).mol"`)
		m = match(r"complexity (\d+)", output)
		if m !== nothing
			println("MA: ", m[1])
			complexities.ma[i] = parse(Float64, m[1])
		end
	end
end
@info "Done calculating MA values"

open("data/complexities_$pathway.json", "w") do io
	JSON3.write(io, complexities)
end