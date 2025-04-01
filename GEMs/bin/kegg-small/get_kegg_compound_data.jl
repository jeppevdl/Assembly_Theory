#!/usr/bin/env julia
#title      : get_kegg_compound_data
#description: Generate JSON with KEGG COMPOUND data
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
using Dates

# Constants
CHUNKSIZE = 500
MAXTRIES = 5

# Get compounds lists
cpd_list = KEGGAPI.list("cpd").data |> first |> x -> map(e -> "cpd:$e", x)
@info "Found $(length(cpd_list)) compound entries for data retrieval."

# Get compound information in chunks
cpd_data = Dict{String, Dict{String, Vector{String}}}()
pbar = Progress(length(cpd_list))
for chunk in [cpd_list[i:min(i + CHUNKSIZE - 1, end)] for i in 1:CHUNKSIZE:length(cpd_list)]
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
            for kv in split.(lines, r"\s+", limit = 2)
                # Handle entries with subkeys
                if length(kv) == 1
                    k = kv[1]
                    v = ""
                else
                    k, v = kv
                end

                # Ensure we have somewhere to put the data
                if k == ""
                    k = lastkey
                end
                if !haskey(data, k)
                    data[k] = String[]
                end

                # Parse specific entries
                if k == "ENTRY"
                    entry = String(first(split(v, " ", keepempty = false)))
                    data[k] = [entry]
                elseif k in ["ENZYME", "REACTION"]
                    v = split(v, " ", keepempty = false)
                    append!(data[k], v)
                elseif v != ""
                    push!(data[k], v)
                end

                # Set current key as last key
                lastkey = k
            end

            # Only add information if entry is valid
            if !isnothing(entry)
                cpd_data["cpd:$entry"] = data
            end

            # Advance progress bar
            next!(pbar)
        catch e
            @warn "Error for compound '$entry':" e
        end
    end
end
@info "Total compounds retrieved from KEGG: $(length(cpd_data))"

# Save information as JSON file
outpath = "./data/metadata/kegg_cpd.v20250320a.json"
open(outpath, "w+") do file
    write(file, JSON.json(cpd_data))
end
@info "Saved data at '$outpath'"

# Generate link to latest version of database
linkpath = "./data/metadata/kegg_cpd.latest.json"
islink(linkpath) && rm(linkpath); symlink(abspath(outpath), abspath(linkpath))
