#!/usr/bin/env julia
#title      : get_kegg_reaction_data
#description: Generate JSON with KEGG REACTION data
#author     : Carlos Vigil-Vásquez
#date       : 2025/03/20
#version    : 20250320a
#notes      : Requires instantiation of project before using (`julia --project=. -e "using Pkg; Pkg.instantiate()`)
#copyright  : Copyright (C) 2025 Carlos Vigil-Vásquez  (carlos.vigil.v@gmail.com)
#license    : Permission to copy and modify is granted under the MIT license

using Pkg
if basename(pwd()) != "kegg-small"
    error("Not in correct directory")
end
Pkg.activate(".")

println(basename(@__FILE__()) * " - v20250320a")
t0 = time()

using KEGGAPI
using JSON
using ProgressMeter

# Constants
CHUNKSIZE = 500
MAXTRIES = 5

# Get reactions lists
rn_list = KEGGAPI.list("rn").data |> first |> x -> map(e -> "rn:$e", x)
@info "Found $(length(rn_list)) reaction entries for data retrieval."

# Get reaction information in chunks
rn_data = Dict{String, Dict{String, Vector{String}}}()
pbar = Progress(length(rn_list))
for chunk in [rn_list[i:min(i + CHUNKSIZE - 1, end)] for i in 1:CHUNKSIZE:length(rn_list)]
    # Try getting information from KEGG's API, if it fails wait a minute before trying again
    response = nothing
    tries = 0
    while isnothing(response)
        try
            response = KEGGAPI.kegg_get(chunk)
        catch e
            @warn e
            @info "Sleeping for a minute"
            sleep(60)
            if tries > MAXTRIES
                @error "Failed $MAXTRIES times to retrieve information, aborting"
                exit(1)
            end
        end
    end

    # Process information
    for entry in response[2]
        try
            lines = split(entry, "\n")

            # Parse lines
            data = Dict{String, Vector{String}}()
            lastkey = ""
            entry = nothing
            for (k, v) in split.(lines, r"\s+", limit = 2)
                # Ensure we have somewhere to put the data
                if k == ""
                    k = lastkey
                end
                if !haskey(data, k)
                    data[k] = String[]
                end

                # Parse specific entries
                if k == "ENTRY"
                    entry = first(split(v, " ", keepempty = false))
                elseif k in ["ENZYME", "REACTION"]
                    v = split(v, " ", keepempty = false)
                    append!(data[k], v)
                else
                    push!(data[k], v)
                end
                lastkey = k
            end

            # Only add information if entry is valid
            if !isnothing(entry)
                rn_data["rn:$entry"] = data
            end
            next!(pbar)
        catch e
            @warn "Error for reaction '$entry':" e
        end
    end
end
@info "Total reactions retrieved from KEGG: $(length(rn_data))"

# Save information as JSON file
outpath = "data/metadata/kegg_rn.v20250320a.json"
open(outpath, "w+") do file
    write(file, JSON.json(rn_data))
end
@info "Saved data at '$outpath'"

# Generate link to latest version of database
linkpath = "data/metadata/kegg_rn.latest.json"
islink(linkpath) && rm(linkpath); symlink(outpath, linkpath)
