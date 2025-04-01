#!/usr/bin/env julia
#title      : AutoKEGGRec.jl
#description: Reconstruct metabolic networks for organisms found in KEGG
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

VERSION = "20250320a"

println(basename(@__FILE__()) * " - v$VERSION")
t0 = time()

using ArgParse
using Graphs
using JSON
using KEGGAPI
using DelimitedFiles
using CSV, DataFrames
using ProgressMeter

function generate_mapper(from, to)
    # Retrieve data
    data = KEGGAPI.link(to, from).data

    # Map "to" to "from"
    mapper = Dict{String, Vector{String}}(from => [] for from in unique(data[1]))
    for (from, to) in zip(data...)
        push!(mapper[from], to)
    end
    return mapper
end

function create_latest_link(path, version = VERSION)
    linkpath = replace(path, "v$version" => "latest")
    islink(linkpath) && rm(linkpath)
    symlink(abspath(path), abspath(linkpath))
    return nothing
end

function main(args)
    # Argument parsing
    settings = ArgParseSettings(
        prog = basename(@__FILE__()),
        description = "Reconstruct metabolic networks for organisms found in KEGG",
        version = VERSION,
        add_version = true,
    )

    add_arg_group!(settings, "I/O option:")
    @add_arg_table! settings begin
        "-i", "--codes"
        help = "List of KEGG ORGANISM codes"
        arg_type = String
        required = true
        action = :store_arg
        "-d", "--delimiter"
		help = "List delimiter (default: `\n`)"
        arg_type = String
        default = "\n"
        action = :store_arg
        "-o", "--output"
		help = "Output directory (default: `./data`)"
        arg_type = String
		default="./data"
        action = :store_arg
    end

    args = parse_args(args, settings)

    # Main body
    @info "Loading organisms to reconstruct metabolic network"
    allorganisms = open(args["codes"], "r") do f
        return String.(split(String(read(f)), args["delimiter"], keepempty = false))
    end
    @info "Found $(length(allorganisms)) organisms as input"
    @info "First 5 organisms: $(join(first(allorganisms, 5), "; "))"

    @info "Creating output directories"
    !isdir(args["output"]) ? mkdir(args["output"]) : @info "Found output directory, skipping creation"
    for p in ["adjmat", "lookup"]
        !isdir(args["output"] * "/$p") ? mkdir(args["output"] * "/$p") : error("Found output directory")
    end

    # Create mappings
    @info "Creating mappings"
    @info "Mapping enzymes to reactions"
    ec2rn = generate_mapper("ec", "rn")
    @info "Mapping reactions to compounds"
    rn2cpd = generate_mapper("rn", "cpd")

    # Get organism information
    stats = []
    @showprogress desc = "Reconstructing metabolic networks..." showspeed = true for organism in allorganisms
        org_genes2ec = generate_mapper(organism, "ec")
        org_ec = values(org_genes2ec) |> collect |> x -> vcat(x...) |> unique |> sort
        org_ec2rn = [get(ec2rn, ec, missing) for ec in org_ec]
        org_rn = unique(vcat(filter(!ismissing, org_ec2rn)...))
        org_cpd = unique(vcat(filter(!ismissing, [get(rn2cpd, rn, missing) for rn in org_rn])...))

        push!(
            stats,
            Dict(
                "org" => organism,
                "ec_coding_genes" => length(org_genes2ec),
                "ec" => length(org_ec),
                "ec_with_rn" => length(filter(!ismissing, org_ec2rn)),
                "rn" => length(org_rn),
                "cpd" => length(org_cpd),
            )
        )

        # Construct metabolic network
        V = sort([org_rn; org_cpd])
        V_idx = enumerate(V) |> e -> reverse.(e) |> Dict
        G = Graph(length(V))
        for rn in org_rn
            rn_pos = get(V_idx, rn, missing)
            cpds = get(rn2cpd, rn, missing)
            if !ismissing(rn_pos) && !ismissing(cpds)
                cpds_pos = map(cpd -> get(V_idx, cpd, missing), cpds)
                rn_e = [(rn_pos, cpd_pos) for cpd_pos in cpds_pos if !ismissing(cpd_pos)]
                map(e -> add_edge!(G, e...), rn_e)
            end
        end

        # Save information
        adjmat_path = args["output"] * "/adjmat/$organism.adjmat_v$VERSION.txt"
        open(adjmat_path, "w+") do f
            writedlm(f, adjacency_matrix(G))
        end
        create_latest_link(adjmat_path)
        lookup_path = args["output"] * "/lookup/$organism.lookup_v$VERSION.txt"
        open(lookup_path, "w+") do f
            writedlm(f, V)
        end
        create_latest_link(lookup_path)
    end

    open(args["output"] * "/report_AutoKEGGRec_v$VERSION.json", "w+") do f
        write(f, JSON.json(stats))
    end
    create_latest_link(args["output"] * "/report_AutoKEGGRec_v$VERSION.json")

    return 0
end

main(ARGS)
