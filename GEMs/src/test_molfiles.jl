using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")
using Glob, ProgressMeter

# CHECK IF MOLFILES IN A DIRECTORY ARE NAMED CORRECTLY AND ASSIGNED TO THE CORRECT ENTRY ID

# Function to extract the entry ID from the <ENTRY> tag
function extract_entry_id(filepath::String)
    open(filepath, "r") do file
        for line in eachline(file)
            if occursin("cpd:", line)
                m = match(r"cpd:(C\d+)", String(line))
                return m.captures[1]
            end
        end
    end
end

# Function to check the filenames and their corresponding entry IDs
mismatches = []
function check_molfile_entries(directory::String)
    mol_files = Glob.glob("C*.mol", directory)
    @showprogress for filepath in mol_files
        filename = split(basename(filepath), ".")[1]
        entry_id = extract_entry_id(filepath)
        if entry_id === nothing
            println("\nNo <ENTRY> tag found in $filepath")
        elseif filename != entry_id
            push!(mismatches, (filename, entry_id))
            println("Mismatch in $filepath: filename = $filename, <ENTRY> = $entry_id")
        end
    end
end

# Run the check
molfile_dir = "..\\bin\\assembly_go\\molfiles"
check_molfile_entries(molfile_dir)

# code below to get the mismatched files from KEGG
using KEGGAPI
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

input = ["cpd:"*m[1] for m in mismatches]
allmol = Dict{String,Union{Missing,String}}()
@showprogress for batch in collect(Iterators.partition(input, 1))
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

@showprogress for (cpd, mol) in allmol
    if mol !== missing
        id = split(cpd, ":")[2]
        open("../data/molfiles/$id.mol", "w") do f
            write(f, mol)
        end
        open("../bin/assembly_go/molfiles/$id.mol", "w") do f
            write(f, mol)
        end
    end
end

deletes = [cpd for cpd in input if !haskey(allmol, cpd)]
for cpd in deletes
    id = split(cpd, ":")[2]
    #remove file
    rm("../data/molfiles/$id.mol")
    rm("../bin/assembly_go/molfiles/$id.mol")
end