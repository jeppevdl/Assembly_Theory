using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using JSON3, KEGGAPI, ProgressMeter, DataFrames, CSV, StatsPlots, Statistics, CairoMakie

# function for gathering KEGG information
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

# function for parsing one side of a reaction equation
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

# function for calculating IQR of a vector
function iqr(v)
    q25 = quantile(v, 0.25)
    q75 = quantile(v, 0.75)
    return q75 - q25
end

# read in all calculated MA values
MAs = CSV.read("../data/bash_MA_output/complete_MAs.csv", DataFrame;)

# read in all reaction classes with their respective main compound pairs
rc2cpd_pairs = JSON3.read("../data/reactions/rc_cpd_pairs.json", Dict{String, Vector{String}})

# read in pathway metadata
pathway_metadata = CSV.read("../data/reactions/pathway_metadata.csv", DataFrame)
pathway_list = pathway_metadata.id
skipped_pathways = ["map00510", "map00532", "map00563", "map00530", "map00190","map00271",
"map00272","map00031","map00150","map00252","map01030","map01031","map01032","map00533","map00602","map00601",
"map00534", "map00603", "map00531", "map00512", "map00604"]

# keep track of missing MA values
missingma = []

# dictionaries for storing stats, unique compounds and reaction class information
stats_dict = Dict{String, Dict{String, Any}}()
unique_cpd_dict = Dict{String, Dict{String, Int}}()
rclass_dict = Dict{String, Dict{String, Dict{String, Any}}}()

origin_df = DataFrame(origin = unique(pathway_metadata.origin), vals = [[] for _ in 1:length(unique(pathway_metadata.origin))],
    diff = [[] for _ in 1:length(unique(pathway_metadata.origin))])

super_df = DataFrame(super_pathway = unique(pathway_metadata.super_pathway), vals = [[] for _ in 1:length(unique(pathway_metadata.super_pathway))],
    diff = [[] for _ in 1:length(unique(pathway_metadata.super_pathway))])

colors = palette(:seaborn_colorblind);

for pathway in pathway_list
    @info "Processing pathway $pathway"
    if pathway in skipped_pathways
        @info "Skipping pathway $pathway"
        continue
    end

    if isfile("../data/reactions/rn_data/kegg_rn_$pathway.json")
        rn_data = open("../data/reactions/rn_data/kegg_rn_$pathway.json", "r") do io
            JSON3.read(io, Dict{String, Dict{String, Vector{String}}})
        end
    else
        pathway_reactions = KEGGAPI.link("rn", String(pathway)).data[2]

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

    # preallocate dictionary for reaction data linked with MA values
    rn_ma = Dict{String, Dict{String, Any}}()
    
    num_reactions = 0
    num_with_rclass = 0

    no_ma = []
    all_in = []
    all_out = []

    # keep track of unique compounds and stats for each pathway
    unique_cpds = Dict{String, Int}()
    stats = Dict{String, Any}()

    # rxs = []
    # t = default_t()

    # go over reactions in pathway
    for entry in keys(rn_data)

        data = rn_data[entry]

        # check if reaction has a class
        num_reactions += 1
        if haskey(data, "RCLASS") && !isempty(data["RCLASS"])
            num_with_rclass += 1
        end

        # check if reaction has an enzyme
        if haskey(data, "ENZYME")
            ec = data["ENZYME"][1][1] |> x -> parse(Int, string(x))
        else
            ec = NaN
        end

        # parse the reaction equation
        equation = data["EQUATION"][1]
        lhs, rhs = split(equation, "<=>")
        ins, _ = parse_side_vec(String(lhs))
        outs, _ = parse_side_vec(String(rhs))

        # species_in = [sp * "(t)" for sp in ins] .|> Meta.parse
        # species_in = [only(@eval @species $(element)) for element in species_in]

        # species_out = [sp * "(t)" for sp in outs] .|> Meta.parse
        # species_out = [only(@eval @species $(element)) for element in species_out]

        # push!(rxs, Reaction(5.0, species_in, species_out, input_coeffs, output_coeffs))

        main_reaction_pairs = []
        reaction_classes = []

        # go over reaction classes present in the pathway and calculate statistics
        if haskey(data, "RCLASS")
            classes = data["RCLASS"]
            for class in classes
                elements = split(class, "  ")
                cl = elements[1]
                push!(reaction_classes, cl)
                pairs = elements[2:end]
                main_pairs = rc2cpd_pairs[cl]
                for pair in pairs
                    if pair in main_pairs # only consider main reaction pairs
                        push!(main_reaction_pairs, cl*":"*pair)
                        if haskey(rclass_dict, cl)
                            if !haskey(rclass_dict[cl], pair)
                                p1 = split(pair, "_")[1]
                                p2 = split(pair, "_")[2]
                                if (p1 in ins && p2 in outs) || (p1 in outs && p2 in ins)
                                    if MAs.ma[MAs.cpd .== p1][1] != "na" && MAs.ma[MAs.cpd .== p2][1] != "na"
                                        ma1 = MAs.ma[MAs.cpd .== p1][1] |> x -> parse(Int, string(x))
                                        ma2 = MAs.ma[MAs.cpd .== p2][1] |> x -> parse(Int, string(x))
                                        mean_ma_pair = mean([ma1, ma2])
                                        diff_ma_pair = abs(ma1 - ma2)
                                        rclass_dict[cl][pair] = Dict("MA" => [ma1, ma2], "MEAN_MA" => mean_ma_pair, "DIFF" => diff_ma_pair, "ENZYME" => ec)
                                    end
                                end
                            end
                        else
                            p1 = split(pair, "_")[1]
                            p2 = split(pair, "_")[2]
                            if (p1 in ins && p2 in outs) || (p1 in outs && p2 in ins)
                                if MAs.ma[MAs.cpd .== p1][1] != "na" && MAs.ma[MAs.cpd .== p2][1] != "na"
                                    ma1 = MAs.ma[MAs.cpd .== p1][1] |> x -> parse(Int, string(x))
                                    ma2 = MAs.ma[MAs.cpd .== p2][1] |> x -> parse(Int, string(x))
                                    mean_ma_pair = mean([ma1, ma2])
                                    diff_ma_pair = abs(ma1 - ma2)
                                    rclass_dict[cl] = Dict(pair => Dict("MA" => [ma1, ma2], "MEAN_MA" => mean_ma_pair, "DIFF" => diff_ma_pair, "ENZYME" => ec))
                                end
                            end
                        end
                    end
                end
            end
        else
            println("No RCLASS for $entry")
            continue
        end

        main_pairs = [string(split(elem, ":")[2]) for elem in main_reaction_pairs]
        main_cpds = unique(vcat([split(elem, "_") for elem in main_pairs]...))
        
        ma_in = []
        ma_out = []

        # gather MA values for input and output compounds that belong to main pairs and add to unique compounds dictionary
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
        
        # calculate MA difference and mean MA
        if any(isnan, ma_in) || any(isnan, ma_out)
            diff = NaN
            mean_ma = NaN
        else
            diff = sum(ma_out) - sum(ma_in)
            mean_ma = mean([mean(ma_in), mean(ma_out)])
        end
        
        all_in = vcat(all_in, ins)
        all_out = vcat(all_out, outs)
        
        # store reaction data in rn_ma dictionary
        rn_ma[entry] = Dict("ENZYME" => ec, "EQUATION" => equation, "IN" => ins, "OUT" => outs, 
        "MA_IN" => ma_in, "MA_OUT" => ma_out, "DIFF" => diff, "MEAN_MA" => mean_ma, 
        "MAIN_PAIRS" => main_reaction_pairs, "RCLASS" => reaction_classes)
    end
    
    push!(missingma, no_ma)
    CSV.write("../data/reactions/rn_ma_stats/rn_stats_$pathway.csv", rn_ma)
    unique_cpd_dict[pathway] = unique_cpds

    stats["num_reactions"] = num_reactions
    stats["num_with_rclass"] = num_with_rclass
    stats["rclass_coverage"] = num_with_rclass / num_reactions

    stats["origin"] = pathway_metadata.origin[pathway_metadata.id .== pathway][1]
    stats["super_pathway"] = pathway_metadata.super_pathway[pathway_metadata.id .== pathway][1]

    pw_name = pathway_metadata.name[pathway_metadata.id .== pathway][1]

    # calculate statistics for the pathway, make histograms and boxplot data based on grouping
    if length(rn_ma) != 0 && !(pathway in skipped_pathways) # skip certain pathways
        mean_diff = mean([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["mean_diff"] = mean_diff

        # NO duplicates are used for MA stats
        mean_ma = mean(collect(values(unique_cpds)))
        stats["mean_ma"] = mean_ma

        sdev_diff = std([value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])])
        stats["sdev_diff"] = sdev_diff

        sdev_ma = std(collect(values(unique_cpds)))
        stats["sdev_ma"] = sdev_ma

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

        h1 = Figure(size = (800, 600));
        ax1 = Axis(h1[1, 1], xlabel = "Absolute MA difference between input and output", ylabel = "Frequency", 
        title = "Distribution of MA difference between input and output for reactions in $pw_name")
        CairoMakie.hist!(ax1, [value["DIFF"] for value in values(rn_ma) if !isnan(value["DIFF"])],
        color = colors[1])
        CairoMakie.save("../figures/pathway_diff_dist/$(pathway)_diff.png", h1, px_per_unit=5)
        
        h2 = Figure(size = (800, 600));
        ax2 = Axis(h2[1, 1], xlabel = "MA", ylabel = "Frequency", title = "Distribution of MA values for $pw_name")
        CairoMakie.hist!(ax2, collect(values(unique_cpds)), color = colors[2])
        CairoMakie.save("../figures/pathway_ma_dist/$(pathway)_ma.png", h2, px_per_unit=5)
        
        ma_bp_data = collect(values(unique_cpds))
        diff_bp_data = [abs(value["DIFF"]) for value in values(rn_ma) if !isnan(value["DIFF"])]

        group_label_origin = pathway_metadata.origin[pathway_metadata.id .== pathway][1]
        group_label_super = pathway_metadata.super_pathway[pathway_metadata.id .== pathway][1]

        # if group_label_super == "Metabolism of cofactors & vitamins"
        #     #print compound entries with ma larger than 59 from unique_cpds
        #     for (cpd, ma) in unique_cpds
        #         if ma > 59
        #             println("Compound $cpd has MA $ma")
        #         end
        #     end
        # end

        append!(origin_df[!, :vals][origin_df.origin .== group_label_origin][1], ma_bp_data)
        append!(origin_df[!, :diff][origin_df.origin .== group_label_origin][1], diff_bp_data)

        append!(super_df[!, :vals][super_df.super_pathway .== group_label_super][1], ma_bp_data)
        append!(super_df[!, :diff][super_df.super_pathway .== group_label_super][1], diff_bp_data)

    end
    stats_dict[pathway] = stats
end

# save stats to CSV
CSV.write("../data/reactions/pathway_stats.csv", stats_dict)
CSV.write("../data/reactions/pathway_cpd_stats.csv", unique_cpd_dict)
CSV.write("../data/reactions/rclass_stats.csv", rclass_dict)

# boxplots and histograms grouped by origin
ma_bp_origin = Figure(size = (800,600));
ax1 = Axis(ma_bp_origin[1,1], title="MA values for all pathways grouped by pathway origin", xlabel="Pathway origin", ylabel="MA", width=600, xticklabelsvisible=false)
i=1
medians = []
for row in eachrow(origin_df)
    ma_values = row.vals
    CairoMakie.boxplot!(ax1, i*ones(length(ma_values)), ma_values; color = colors[i], label=row.origin, width=0.5)
    n = length(ma_values)
    CairoMakie.text!(ax1, i, -4, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(ma_values))

    histo = Figure(size=(800,600));
    ax = Axis(histo[1, 1], xlabel = "MA", ylabel = "Frequency", title = "Distribution of MA values for all pathways in $(row.origin) origin")
    CairoMakie.hist!(ax, ma_values, color = colors[i])
    CairoMakie.save("../figures/pathway_ma_dist/origin/$(row.origin)_ma.png", histo, px_per_unit=5)
    i += 1
end

CairoMakie.Legend(ma_bp_origin[1,2], ax1, "Pathway origin")
ma_bp_origin

CairoMakie.save("../figures/pathway_ma_dist/origin/MA_boxplot_origin.png", ma_bp_origin, px_per_unit=5)
CSV.write(("../figures/pathway_ma_dist/origin/MA_medians.csv"), DataFrame(origin = origin_df.origin, median = medians))

diff_bp_origin = Figure(size = (800,600));
ax2 = Axis(diff_bp_origin[1,1], title="Absolute differences in MA between compound pairs in pathway reactions grouped by pathway origin", xlabel="Pathway origin", ylabel="Absolute MA difference", width=600, xticklabelsvisible=false)
i=1
medians = []
for row in eachrow(origin_df)
    diff_values = row.diff
    CairoMakie.boxplot!(ax2, i*ones(length(diff_values)), diff_values; color = colors[i], label=row.origin, width=0.5)
    n = length(diff_values)
    CairoMakie.text!(ax2, i, -2, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(diff_values))

    histo = Figure(size=(800,600));
    ax = Axis(histo[1, 1], xlabel = "Absolute MA difference", ylabel = "Frequency", title = "Distribution of absolute MA differences for all pathway reactions in $(row.origin) origin")
    CairoMakie.hist!(ax, diff_values, color = colors[i])
    CairoMakie.save("../figures/pathway_diff_dist/origin/$(row.origin)_diff.png", histo, px_per_unit=5)
    i += 1
end

CairoMakie.Legend(diff_bp_origin[1,2], ax2, "Pathway origin")
diff_bp_origin

CairoMakie.save("../figures/pathway_diff_dist/origin/diff_boxplot_origin.png", diff_bp_origin, px_per_unit=5)
CSV.write(("../figures/pathway_diff_dist/origin/diff_medians.csv"), DataFrame(origin = origin_df.origin, median = medians))

# boxplots and histograms grouped by super pathway (leave out glycans)
super_df = super_df[super_df.super_pathway .!= "Glycan biosynthesis & metabolism", :]
ma_bp_super = Figure(size = (1000,800));
ax3 = Axis(ma_bp_super[1,1], title="MA values for all pathways grouped by super pathway", xlabel="Super pathway", ylabel="MA", width=600, xticklabelsvisible=false)
i=1
medians = []
for row in eachrow(super_df)
    ma_values = row.vals
    CairoMakie.boxplot!(ax3, i*ones(length(ma_values)), ma_values; color = colors[i], label=row.super_pathway, width=0.5)
    n = length(ma_values)
    CairoMakie.text!(ax3, i, -2, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(ma_values))

    histo = Figure(size=(800,600));
    ax = Axis(histo[1, 1], xlabel = "MA", ylabel = "Frequency", title = "Distribution of MA values for all pathways in $(row.super_pathway)")
    CairoMakie.hist!(ax, ma_values, color = colors[i])
    CairoMakie.save("../figures/pathway_ma_dist/super_pathway/$(row.super_pathway)_ma.png", histo, px_per_unit=5)
    i += 1
end

CairoMakie.Legend(ma_bp_super[1,2], ax3, "Super pathway")
ma_bp_super

CairoMakie.save("../figures/pathway_ma_dist/super_pathway/MA_boxplot_super.png", ma_bp_super, px_per_unit=5)
CSV.write(("../figures/pathway_ma_dist/super_pathway/MA_medians.csv"), DataFrame(super_pathway = super_df.super_pathway, median = medians))

diff_bp_super = Figure(size = (1000,800));
ax4 = Axis(diff_bp_super[1,1], title="Absolute differences in MA between compound pairs in pathway reactions grouped by super pathway", xlabel="Super pathway", ylabel="Absolute MA difference", width=600, xticklabelsvisible=false)
i=1
medians = []
for row in eachrow(super_df)
    diff_values = row.diff
    CairoMakie.boxplot!(ax4, i*ones(length(diff_values)), diff_values; color = colors[i], label=row.super_pathway, width=0.5)
    n = length(diff_values)
    CairoMakie.text!(ax4, i, -2, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(diff_values))

    histo = Figure(size=(800,600));
    ax = Axis(histo[1, 1], xlabel = "Absolute MA difference", ylabel = "Frequency", title = "Distribution of absolute MA differences for all pathway reactions in $(row.super_pathway)")
    CairoMakie.hist!(ax, diff_values, color = colors[i])
    CairoMakie.save("../figures/pathway_diff_dist/super_pathway/$(row.super_pathway)_diff.png", histo, px_per_unit=5)
    i += 1
end

CairoMakie.Legend(diff_bp_super[1,2], ax4, "Super pathway")
diff_bp_super

CairoMakie.save("../figures/pathway_diff_dist/super_pathway/diff_boxplot_super.png", diff_bp_super, px_per_unit=5)
CSV.write(("../figures/pathway_diff_dist/super_pathway/diff_medians.csv"), DataFrame(super_pathway = super_df.super_pathway, median = medians))

# Create bubble chart for reaction classes
function count_unique_ec_classes(pairs)
    ec_classes = Int[]
    for pairinfo in values(pairs)
        if haskey(pairinfo, "ENZYME")
            enz = pairinfo["ENZYME"]
            if !(ismissing(enz) || isnan(enz))
                push!(ec_classes, Int(enz))
            end
        end
    end
    return length(unique(ec_classes))
end

# Collect data
x_vals = Float64[]
y_vals = Float64[]
sizes = Float64[]
colors = Int[]
labels = String[]

for (rclass, pair_dict) in rclass_dict
    diffs = [abs(pairinfo["DIFF"]) for pairinfo in values(pair_dict)]
    mas = vcat([pairinfo["MA"] for pairinfo in values(pair_dict)]...)
    push!(x_vals, median(diffs))
    push!(y_vals, median(mas))
    push!(sizes, iqr(diffs))
    push!(colors, length(diffs))
    push!(labels, rclass)
end

scaled_sizes = sizes
log_colors = log10.(colors .+ 1)

f = Figure(size = (900, 650));
ax = Axis(f[1, 1],
    xlabel = "Median abs(MA difference)",
    ylabel = "Median MA",
    title = "Reaction Class Bubble Chart"
)

CairoMakie.scatter!(
    ax, x_vals, y_vals;
    markersize = scaled_sizes,
    color = log_colors,
    colormap = :viridis,
    colorrange = (minimum(log_colors), maximum(log_colors)),
    transparency = true
)

Colorbar(f[1, 2],
    limits = (minimum(log_colors), maximum(log_colors)),
    colormap = :viridis,
    label = "log₁₀(Compound pairs + 1)"
)

for (x, y, size, label, color) in zip(x_vals, y_vals, sizes, labels, log_colors)
    if color >= log10(6) && (x > 15 || y > 35)
        CairoMakie.text!(ax, string(label), position = (x, y), align = (:left, :center), fontsize = 10)
    end
end

f

CairoMakie.save("../figures/reaction_analysis/reaction_class_bubble_chart.png", f, px_per_unit=5)

#study enzyme classes
class_names = ["Oxidoreductases", "Transferases", "Hydrolases", "Lyases", "Isomerases", "Ligases", "Translocases"]
ec_info = DataFrame(id = 1:7, name = class_names, all_ma = [[] for _ in 1:7], mean_ma = [[] for _ in 1:7], diff = [[] for _ in 1:7], num_pairs = zeros(7))
for (class, pair_dict) in rclass_dict
    for (pair, pair_info) in pair_dict
        if haskey(pair_info, "ENZYME")
            ec = pair_info["ENZYME"]
            if !isnan(ec)
                append!(ec_info[ec_info.id .== ec, "all_ma"][1], pair_info["MA"])
                append!(ec_info[ec_info.id .== ec, "mean_ma"][1], pair_info["MEAN_MA"])
                append!(ec_info[ec_info.id .== ec, "diff"][1], pair_info["DIFF"])
            end
        end
    end
end
ec_info.num_pairs = ec_info.num_pairs .+ length.(ec_info.all_ma)./2

f = Figure(size=(900,600));
ax = Axis(f[1,1], title="Distribution of mean MA values of compound pairs per enzyme class", xlabel = "Enzyme class", ylabel = "Mean MA",
xticks = 1:7, width = 600)
i=1
colors = palette(:seaborn_colorblind)[1:7];
medians = []
for row in eachrow(ec_info)  
    mean_mas = row.mean_ma
    CairoMakie.boxplot!(ax, i*ones(length(mean_mas)), mean_mas; color = colors[i], label=row.name, width=0.5)
    n = length(mean_mas)
    CairoMakie.text!(ax, i, maximum(mean_mas) + 3, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(mean_mas))
    i += 1
end
CairoMakie.Legend(f[1,2], ax, "Enzyme class")
f

CairoMakie.save("../figures/reaction_analysis/ec_mean_ma.png", f, px_per_unit=5)
CSV.write(("../figures/reaction_analysis/ec_mean_ma_medians.csv"), DataFrame(ec_class = ec_info.name, median = medians))

f2 = Figure(size=(900,600));
ax2 = Axis(f2[1,1], title="Distribution of absolute MA differences of compound pairs per enzyme class", xlabel = "Enzyme class", ylabel = "Absolute MA difference",
xticks = 1:7, width = 600)
i=1
colors = palette(:seaborn_colorblind)[1:7];
medians = []
for row in eachrow(ec_info)  
    diffs = row.diff
    if length(diffs) < 5
        # Show small groups as points
        CairoMakie.scatter!(ax2, i*ones(length(diffs)), diffs; color = colors[i], marker=:circle, markersize=10, label=row.name)
        n = length(diffs)
        if n > 0
            CairoMakie.text!(ax2, i, maximum(diffs) + 3, text = "n=$n", align = (:center, :bottom), color = :black)
            push!(medians, median(diffs))
        else
            push!(medians, NaN)
        end
        i += 1
        continue
    end
    CairoMakie.boxplot!(ax2, i*ones(length(diffs)), diffs; color = colors[i], label=row.name, width=0.5)
    n = length(diffs)
    CairoMakie.text!(ax2, i, maximum(diffs) + 3, text = "n=$n", align = (:center, :bottom), color = :black)
    push!(medians, median(diffs))
    i += 1
end
CairoMakie.Legend(f2[1,2], ax2, "Enzyme class")
f2

CairoMakie.save("../figures/reaction_analysis/ec_diff_ma.png", f2, px_per_unit=5)
CSV.write(("../figures/reaction_analysis/ec_diff_ma_medians.csv"), DataFrame(ec_class = ec_info.name, median = medians))

# allmissing = []
# for element in missingma
#     if length(element) > 0
#         for cpd in element
#             if cpd in allmissing
#                 continue
#             else
#                 push!(allmissing, cpd)
#             end
#         end
#     end
# end
# println(allmissing)

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