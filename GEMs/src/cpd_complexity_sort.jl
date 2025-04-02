using Pkg
if basename(pwd()) != "GEMs"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs")
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
	batch_mol = nothing
	tries = 0
	while isnothing(batch_mol)
		try
			sleep(1)
			batch_mol = KEGGAPI.kegg_get(String.(batch), "mol", 5) |> r -> split(r[2][1], "\$\$\$\$\n", keepempty=false)
			map(zip(batch,batch_mol)) do (cpd,mol)
				allmol[cpd] = mol
			end
		catch e
			@warn e
			@info "Sleeping for a minute"
            sleep(60)
			if tries > 2
				@error "Failed 2 times to retrieve information, aborting"
				map(zip(batch,batch_mol)) do (cpd,mol)
					allmol[cpd] = missing
				end
				exit(1)
			end
			tries += 1
		end
	end
end
@info "Found $(length(filter(x -> !ismissing(x), allmol))) out of $(length(allcpd)) compounds with mol files"

# Calculate complexity of compounds
complexities = DataFrame(id=String[], n_atoms=Dict{Symbol,Int}[], hybr=Vector{Int}[], n_rings=Vector{Int}[])
@info "Calculating complexity of compounds..."
@showprogress for (cpd, mol) in allmol
    if mol !== missing
		id = split(cpd, ":")[2]
        if !isfile("data/molfiles/$id.mol")
			open("data/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
			open("assembly_go/molfiles/$id.mol", "w") do f
				write(f, mol)
			end
		end
		try 
			println(id)
			graph = MolecularGraph.sdftomol("data/molfiles/$id.mol")
			
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
# code below to read in precalculated complexities
# json_data = JSON3.read(open("data/pathway_complexities/complexities_$pathway.json", "r"))
# columns = json_data[:columns] 
# complexities = DataFrame(columns, json_data[:colindex][:names])

# skip this if we calculate all compounds with the same algorithm
# # Get MA values for each compound
# @info "Retrieving MA values for compounds..."
# complexities.ma = Vector{Union{Missing,Float64}}(missing, nrow(complexities))
# @showprogress for i in 1:nrow(complexities)
# 	cpd = complexities.id[i]
# 	if haskey(ma_values, cpd)
# 		complexities.ma[i] = ma_values[cpd]["MA"]
# 	end
# end

# n = count(x -> !ismissing(x), complexities.ma);
# @info "Found MA values for $n compounds"

# cd("assembly_go")

# function run_with_timeout(command::Cmd, timeout::Float64)

#     stdout_buffer = IOBuffer()
#     stderr_buffer = IOBuffer()

#     p = run(pipeline(command; stdout = stdout_buffer, stderr = stderr_buffer), wait = false)

#     t_start = time()
#     # check if output is returned or time limit is exceeded
#     while isempty(readchomp(seekstart(stdout_buffer))) && (time() - t_start < timeout)
#         sleep(0.1)
#     end

#     if isopen(p)
#         kill(p, Base.SIGINT)
#     end

#     # HACK: It seems one has to access the data before converting it to a String, unsure why
#     println(stdout_buffer.data)
#     println(stderr_buffer.data)

#     stdout_str = deepcopy(String(stdout_buffer.data))
#     stderr_str = deepcopy(String(stderr_buffer.data))

#     return stdout_str, stderr_str
# end

# @info "Calculating MA values using molecular assembly algorithm..."
# @showprogress for i in nrow(complexities):-1:1
# 	println(complexities.id[i])
# 	cmd = `./assembly.exe molfiles/$(complexities.id[i]).mol`
# 	stdout, stderr = run_with_timeout(cmd, 10.0)
# 	#...
# end
# @info "Done calculating MA values"

open("data/complexities_$pathway.json", "w") do io
	JSON3.write(io, complexities)
end