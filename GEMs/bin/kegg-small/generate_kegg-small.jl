#!/usr/bin/env julia
#title      : generate_kegg-small
#description: Generate the "kegg-small" subset
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
using DataFrames
using CSV

"""
    split_phylogeny(r)

Split a phylogeny string into levels.

This function takes a phylogeny string `r` separated by semicolons and splits it into
individual levels. If the resulting array has fewer than 4 elements, it is padded with
"-" to ensure a length of 4.

# Arguments
- `r::String`: A phylogeny string with levels separated by semicolons.

# Returns
- `Vector{String}`: An array of phylogeny levels, always with a length of 4.
"""
function split_phylogeny(r::String)::Vector{String}
    lvls = String.(split(r, ";"))
    if length(lvls) < 4
        append!(lvls, repeat(["-"], 4 - length(lvls)))
    end
    return lvls
end

@info "Getting data from KEGG API ($(trunc(time() - t0; digits = 3)) seconds)"
result = KEGGAPI.list("organism")

@info "Converting to dataframe ($(trunc(time() - t0; digits = 3)) seconds)"
df = DataFrame(result.data, result.colnames)

@info "Splitting taxonomy into columns ($(trunc(time() - t0; digits = 3)) seconds)"
transform!(df, :Phylogeny => ByRow(r -> split_phylogeny(r)) => [:Superkingdom, :Kingdom, :Phylum, :Genus])

@info "Filtering out organisms with incomplete taxonomy ($(trunc(time() - t0; digits = 3)) seconds)"
filter!(r -> r.Genus != "-", df)

@info "Grouping by phylum ($(trunc(time() - t0; digits = 3)) seconds)"
gdf = groupby(df, :Phylum)

@info "Keeping one orgnanism per group ($(trunc(time() - t0; digits = 3)) seconds)"
sdf = combine(gdf, x -> first(x, 1))

@info "Cleaning-up table ($(trunc(time() - t0; digits = 3)) seconds)"
sdf = sdf[!, ["Organism", "Species", "Superkingdom", "Kingdom", "Phylum", "Genus", "T. number"]]

@info "Writing table ($(trunc(time() - t0; digits = 3)) seconds)"
CSV.write("data/kegg-small.tsv", sdf; delim="\t")
