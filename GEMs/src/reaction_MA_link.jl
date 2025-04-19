using Pkg

# cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/AT_exploration/")
# Pkg.activate(".")
# using Catalyst, ModelingToolkit, DifferentialEquations

if basename(pwd()) != "src"
    cd("C:/Users/jeppe/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs/src")
end
Pkg.activate(".")

using JSON3, KEGGAPI, ProgressMeter, DataFrames, CSV, Plots, StatsPlots, Statistics

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

function parse_side_vec(eqs::String)
    parts = split(eqs, r"\+")
    compounds = String[]
    coefficients = Int[]
    for part in parts
        clean = strip(part)
        m = match(r"(?:(\d+)\s+)?(C\d{5})", clean)
        if m !== nothing
            coeff = m.captures[1] === nothing ? 1 : parse(Int, m.captures[1])
            cid = m.captures[2]
            push!(compounds, cid)
            push!(coefficients, coeff)
        end
    end
    return compounds, coefficients
end

MAs = CSV.read("../data/bash_MA_output/complete_MAs.csv", DataFrame;)
rc2cpd_pairs = JSON3.read("../data/reactions/rc_cpd_pairs.json", Dict{String, Vector{String}})

missingma = []

pathway_metadata = CSV.read("../data/reactions/pathway_metadata.csv", DataFrame)
pathway_list = pathway_metadata.id

stats_dict = Dict{String, Dict{String, Any}}()
unique_cpd_dict = Dict{String, Dict{String, Int}}()

all_ma_values = Float64[]
all_ma_groups = String[]

all_diff_values = Float64[]
all_diff_groups = String[]

for pathway in pathway_list
    @info "Processing pathway $pathway"

    if isfile("../data/reactions/rn_data/kegg_rn_$pathway.json")
        rn_data = open("../data/reactions/rn_data/kegg_rn_$pathway.json", "r") do io
            JSON3.read(io, Dict{String, Dict{String, Vector{String}}})
        end
    else
        pathway_reactions = KEGGAPI.link("rn", pathway).data[2]

        rn_data = Dict{String, Dict{String, Vector{String}}}()

        @showprogress for chunk in pathway_reactions
            # Try getting information from KEGG's API, if it fails wait a minute before trying again
            sleep(0.4)
            response = nothing
            tries = 0
            while isnothing(response)
                try
                    response = KEGGAPI.kegg_get([chunk])
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
                catch e
                    @warn "Error for reaction '$entry':" e
                end
            end
        end
        @info "Total reactions retrieved from KEGG: $(length(rn_data)) out of $(length(pathway_reactions))"

        # Save information as JSON file
        open("../data/reactions/rn_data/kegg_rn_$pathway.json", "w") do io
            JSON3.write(io, rn_data)
        end
    end

    rn_ma = Dict{String, Dict{String, Any}}()
    
    num_reactions = 0
    num_with_rclass = 0

    no_ma = []
    all_in = []
    all_out = []

    unique_cpds = Dict{String, Int}()
    stats = Dict{String, Any}()

    # rxs = []
    # t = default_t()

    for entry in keys(rn_data)
        # println("Processing $entry")
        data = rn_data[entry]

        num_reactions += 1
        if haskey(data, "RCLASS") && !isempty(data["RCLASS"])
            num_with_rclass += 1
        end

        if haskey(data, "ENZYME")
            ec = data["ENZYME"][1][1] |> x -> parse(Int, string(x))
        else
            ec = NaN
        end

        main_reaction_pairs = []
        reaction_classes = []
        if haskey(data, "RCLASS")
            classes = data["RCLASS"]
            for class in classes
                elements = split(class, "  ")
                cl = elements[1]
                push!(reaction_classes, cl)
                pairs = elements[2:end]
                main_pairs = rc2cpd_pairs[cl]
                for pair in pairs
                    if pair in main_pairs
                        push!(main_reaction_pairs, cl*":"*pair)
                    end
                end
            end
        else
            println("No RCLASS for $entry")
            continue
        end

        equation = data["EQUATION"][1]
        lhs, rhs = split(equation, "<=>")
        ins, _ = parse_side_vec(String(lhs))
        outs, _ = parse_side_vec(String(rhs))

        # species_in = [sp * "(t)" for sp in ins] .|> Meta.parse
        # species_in = [only(@eval @species $(element)) for element in species_in]

        # species_out = [sp * "(t)" for sp in outs] .|> Meta.parse
        # species_out = [only(@eval @species $(element)) for element in species_out]

        # push!(rxs, Reaction(5.0, species_in, species_out, input_coeffs, output_coeffs))

        ma_in = []
        ma_out = []

        main_pairs = [string(split(elem, ":")[2]) for elem in main_reaction_pairs]
        main_cpds = unique(vcat([split(elem, "_") for elem in main_pairs]...))
        
        for cpd_in in ins 
            if cpd_in in main_cpds
                if cpd_in in MAs.cpd
                    if MAs.ma[MAs.cpd .== cpd_in][1] != "na"
                        push!(ma_in, MAs.ma[MAs.cpd .== cpd_in][1] |> x -> parse(Int, string(x)))
                        if !haskey(unique_cpds, cpd_in)
                            unique_cpds[cpd_in] = MAs.ma[MAs.cpd .== cpd_in][1] |> x -> parse(Int, string(x))
                        end
                    else
                        push!(ma_in, NaN)
                        println("Missing MA for $cpd_in in $entry")
                        push!(no_ma, cpd_in)
                    end
                else
                    push!(ma_in, NaN)
                    println("Missing MA for $cpd_in in $entry")
                    push!(no_ma, cpd_in)
                end
            end
        end

        for cpd_out in outs
            if cpd_out in main_cpds
                if cpd_out in MAs.cpd
                    if MAs.ma[MAs.cpd .== cpd_out][1] != "na"
                        push!(ma_out, MAs.ma[MAs.cpd .== cpd_out][1] |> x -> parse(Int, string(x)))
                        if !haskey(unique_cpds, cpd_out)
                            unique_cpds[cpd_out] = MAs.ma[MAs.cpd .== cpd_out][1] |> x -> parse(Int, string(x))
                        end
                    else
                        push!(ma_out, NaN)
                        println("Missing MA for $cpd_out in $entry")
                        push!(no_ma, cpd_out)
                    end
                else
                    push!(ma_out, NaN)
                    println("Missing MA for $cpd_out in $entry")
                    push!(no_ma, cpd_out)
                end
            end
        end
        
        if any(isnan, ma_in) || any(isnan, ma_out)
            diff = NaN
            mean_ma = NaN
        else
            diff = sum(ma_out) - sum(ma_in)
            mean_ma = mean([mean(ma_in), mean(ma_out)])
        end
        
        all_in = vcat(all_in, ins)
        all_out = vcat(all_out, outs)
        
        rn_ma[entry] = Dict("ENZYME" => ec, "EQUATION" => equation, "IN" => ins, "OUT" => outs, 
        "MA_IN" => ma_in, "MA_OUT" => ma_out, "DIFF" => diff, "MEAN_MA" => mean_ma, "MAIN_PAIRS" => main_reaction_pairs, "RCLASS" => reaction_classes)
    end
    
    push!(missingma, no_ma)

    unique_cpd_dict[pathway] = unique_cpds

    stats["num_reactions"] = num_reactions
    stats["num_with_rclass"] = num_with_rclass
    stats["rclass_coverage"] = num_with_rclass / num_reactions

    stats["origin"] = pathway_metadata.origin[pathway_metadata.id .== pathway][1]
    stats["super_pathway"] = pathway_metadata.super_pathway[pathway_metadata.id .== pathway][1]

    if length(rn_ma) != 0 && !(pathway in ["map00510", "map00532", "map00563"])
        mean_diff = mean([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["mean_diff"] = mean_diff

        mean_ma = mean(collect(values(unique_cpds)))
        stats["mean_ma"] = mean_ma

        median_diff = median([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["median_diff"] = median_diff

        median_ma = median(collect(values(unique_cpds)))
        stats["median_ma"] = median_ma

        min_diff = minimum([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["min_diff"] = min_diff

        min_ma = minimum(collect(values(unique_cpds)))
        stats["min_ma"] = min_ma

        max_diff = maximum([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["max_diff"] = max_diff

        max_ma = maximum(collect(values(unique_cpds)))
        stats["max_ma"] = max_ma

        q1_diff = quantile([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])], 0.25)
        stats["q1_diff"] = q1_diff

        q1_ma = quantile(collect(values(unique_cpds)), 0.25)
        stats["q1_ma"] = q1_ma

        q3_diff = quantile([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])], 0.75)
        stats["q3_diff"] = q3_diff

        q3_ma = quantile(collect(values(unique_cpds)), 0.75)
        stats["q3_ma"] = q3_ma

        h1 = histogram([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])], dpi = 600, size = (800,600), legend=false, xlabel = "MA difference between input and output", ylabel = "Frequency", title = "Distribution of MA difference between input and output for $pathway", alpha = 0.75, color = :blue);
        savefig(h1, "../figures/pathway_diff_dist/$(pathway)_diff.png")
        h2 = histogram(collect(values(unique_cpds)), dpi = 600, size = (800,600), legend=false, xlabel = "MA value", ylabel = "Frequency", title = "Distribution of MA values for $pathway", alpha = 0.75, color = :red)
        savefig(h2, "../figures/pathway_ma_dist/$(pathway)_ma.png")
        
        ma_bp_data = collect(values(unique_cpds))
        diff_bp_data = [abs(value["DIFF"]) for value in values(rn_ma) if !isnan(value["DIFF"])]

        group_label = pathway_metadata.origin[pathway_metadata.id .== pathway][1]

        append!(all_ma_values, ma_bp_data)
        append!(all_ma_groups, fill(group_label, length(ma_bp_data)))

        append!(all_diff_values, diff_bp_data)
        append!(all_diff_groups, fill(group_label, length(diff_bp_data)))
    end
    stats_dict[pathway] = stats
end

ma_bp = boxplot(all_ma_groups, all_ma_values; dpi = 600, size = (800,800),xlabel = "Pathway lineage of origin", ylabel = "MA value", 
    title = "Distribution of MA values for all pathways", alpha = 0.75, legend = false, xticks = :auto, xrotation = 45)
savefig(ma_bp, "../figures/pathway_ma_dist/boxplot_ma_grouped.png")

diff_bp = boxplot(all_diff_groups, all_diff_values; dpi = 600, size = (800,800), xlabel = "Pathway lineage of origin", ylabel = "MA difference",
    title = "Distribution of MA differences for all pathways", alpha = 0.75, legend = false, xticks = :auto, xrotation = 45)
savefig(diff_bp, "../figures/pathway_diff_dist/boxplot_diff_grouped.png")

CSV.write("../data/reactions/pathway_stats.csv", stats_dict)

allmissing = []
for element in missingma
    if length(element) > 0
        for cpd in element
            if cpd in allmissing
                continue
            else
                push!(allmissing, cpd)
            end
        end
    end
end
println(allmissing)

# group rn_ma by "ENZYME"
# rn_ma_grouped = Dict{String, Vector{Dict{String, Any}}}()
# for entry in keys(rn_ma)
#     data = rn_ma[entry]
#     if haskey(rn_ma_grouped, string(data["ENZYME"]))
#         push!(rn_ma_grouped[string(data["ENZYME"])], data)
#     else
#         rn_ma_grouped[string(data["ENZYME"])] = [data]
#     end
# end

# ecdata = rn_ma_grouped["4"]
# ecvector = []
# for entry in ecdata
#     push!(ecvector, entry["DIFF"])
# end
# ecvector = [element for element in ecvector if !isnan(element)]
# density(ecvector)


# all_in = Set(all_in)
# all_out = Set(all_out)

# starting_cpds = setdiff(all_in, all_out)
# ending_cpds = setdiff(all_out, all_in)

# @named rn = ReactionSystem(rxs, t);
# rn = complete(rn)

# u0_str = [replace(string(sp), "(t)" => "") for sp in unknowns(rn)]
# u0 = [Symbol(sp) => (sp in starting_cpds ? 500. : 0.) for sp in u0_str]

# oprob = ODEProblem(rn, u0, (0, 50.), [])
# sol = solve(oprob, Tsit5())

# p = plot(title = "Concentration of starting compounds over time", xlabel = "Time", ylabel = "Concentration", dpi = 600);
# colors = cgrad(:jet, length(starting_cpds), categorical = true);

# for (i, start_cpd) in enumerate(["C00135", "C00047"])
#     idx = findfirst(x -> x == start_cpd, u0_str)
#     solution = []
#     for j in 1:length(sol.t)
#         push!(solution, sol.u[j][idx])
#     end
#     plot!(p, sol.t, solution, label = start_cpd, color=colors[i])
# end

# p2 = plot(title = "Concentration of output compounds over time", xlabel = "Time", ylabel = "Concentration", dpi = 600);
# colors2 = cgrad(:jet, length(collect(ending_cpds)), categorical = true);

# for (i, end_cpd) in enumerate(["C00135", "C00047", "C00078"])
#     idx = findfirst(x -> x == end_cpd, u0_str)
#     solution = []
#     for j in 1:length(sol.t)
#         push!(solution, sol.u[j][idx])
#     end
#     plot!(p2, sol.t, solution, label = end_cpd, color=colors2[i])
# end