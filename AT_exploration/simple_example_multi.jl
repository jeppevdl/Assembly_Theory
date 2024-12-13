using Catalyst, OrdinaryDiffEq, Plots, GraphRecipes, JumpProcesses, DataFrames, Statistics, Base.Threads
#functions--------------------------------------------------------------------------------

#function to create a reaction network
function create_reactions(building_blocks, upper_bound, rate; selection_target="", selection_rate=2, upreg_reactions = [])
    species = building_blocks
    # create all possible species under upper bound by combining building blocks
    for i in 2:upper_bound 
        combinations = vec(map(collect, Iterators.product(ntuple(_ -> building_blocks, i)...)))
        spec = [join(c) for c in combinations]
        species = vcat(species, spec)
    end

    # create expressions
    species = [sp * "(t)" for sp in species] .|> Meta.parse
    species = [only(@eval @species $(element)) for element in species]

    # create all possible reactions
    reactions = []
    for str1 in species
        for str2 in species
            len_comb = length(string(str1)*string(str2)) - 6
            if len_comb <= upper_bound
                comb = replace(string(str1), "(t)" => "") * string(str2) .|> Meta.parse
                comb = only(@eval @species $(comb))
                if selection_target != "" && [string(str1), string(str2)] in upreg_reactions
                    println(string(str1), string(str2))
                    if string(str1) == string(str2) 
                        push!(reactions, Reaction(selection_rate, [str1], [comb], [2], [1]))
                    else
                        push!(reactions, Reaction(selection_rate, [str1, str2], [comb], [1, 1], [1]))
                    end
                elseif string(str1) == string(str2) 
                    push!(reactions, Reaction(rate, [str1], [comb], [2], [1]))
                else
                    push!(reactions, Reaction(rate, [str1, str2], [comb], [1, 1], [1]))
                end
            end
        end
    end
    return reactions, species
end

# Get assembly index of target `t` (from assembly_indexV3.jl)
function assembly_index(t; return_visited = false)
    A, visited = _assembly_index(t, String[])
    if return_visited
        return A, visited
    else
        return A
    end
end

# Function calculating the assembly index of a compound `s`
function _assembly_index(s, visited) # aka "AssemblyIndex"
    if s in visited
        # no construction cost for previously visited compounds
        return 0, visited 
    elseif length(s) <= 3
        # strings of length <= 3 can only be constructed in a set amount of combinations
        return length(s) - 1, [visited; s]
    else
        # not yet visited "complex" compound => split into substructures
        substructure_pairs = getsubstructures(s)

        As = Vector{Int64}(undef, length(substructure_pairs))
        visiteds = Vector{Vector{String}}(undef, length(substructure_pairs))

        for (sp_idx, sp) in enumerate(substructure_pairs)
            A1, visited1 = _assembly_index(sp[1], visited)
            A2, visited2 = _assembly_index(sp[2], [visited1; sp[1]])
            As[sp_idx] = A1 + A2 + 1
            visiteds[sp_idx] = [visited2; sp[2]]
        end

        min_idx = argmin(As)
        return As[min_idx], visiteds[min_idx]
    end
end

function getsubstructures(s)
    return [[s[1:i], s[i+1:end]] for i in 1:(length(s)-1)]
end

function assembly_vector(sol, df)
    assembly_vector = []
    for i in 1:length(sol.t)
        NT = sum(sol.u[i])
        A = 0
        for (j, ai) in enumerate(df.ai)
            if sol.u[i][j] != 0
                a = exp(ai) * (sol.u[i][j]-1) / NT
                A += a
            end
        end
        append!(assembly_vector, A)
    end
    return assembly_vector
end

function target_vector(sol, species_strings, target)
    target_idx = findfirst(item -> item == target, species_strings)
    target_vector = []
    for i in 1:length(sol.t)
        append!(target_vector, sol.u[i][target_idx])
    end
    return target_vector
end
    
#-----------------------------------------------------------------------------------------
# setting up model

building_blocks = ["A", "B", "N"]
upper_bound = 6
rate = 1
selection_rate = 10

#set to "" for no selection
selection_target = ""

upregulated = [["B(t)", "A(t)"], ["N(t)", "A(t)"], ["NA(t)", "NA(t)"], ["BA(t)", "NANA(t)"]]
t = default_t()

#creating the reactions can take a long time
reactions, species = create_reactions(building_blocks, upper_bound, rate; selection_target, selection_rate, upreg_reactions = upregulated)
@named jumpmodel = ReactionSystem(reactions, t)

u0 = vcat([Symbol(replace(string(sp), "(t)" => "")) => 5000 for sp in species[1:length(building_blocks)]], [Symbol(replace(string(sp), "(t)" => "")) => 0 for sp in species[length(building_blocks)+1:end]])
tspan = (0.0, 0.02)
ps = []

prob = DiscreteProblem(complete(jumpmodel), u0, tspan, ps)

jump_prob = JumpProblem(complete(jumpmodel), prob, Direct())

#-----------------------------------------------------------------------------------------
#create dataframe with names, lengths and assembly indices for plotting
species_strings = [replace(string(sp), "(t)" => "") for sp in jumpmodel.unknowns]
sol_lengths = map(length, species_strings)
sol_ai = map(assembly_index, species_strings)

df = DataFrame(name = species_strings, length = sol_lengths, ai = sol_ai)

num_sims = 50

# lengths = []
# ais = []
# assemblies = []
# targets = []
# times = []
target_values = []

p_a = plot(title = "Assembly Through Time", xlabel = "Time", ylabel = "Assembly", dpi = 600, legend=false, ylims = (0, 80));
p_t = plot(title = "Target Species Through Time", xlabel = "Time", ylabel = "Target Species", dpi = 600, legend=false, ylims = (0, 400));

@time @threads for _ in 1:num_sims
    
    sol = solve(jump_prob, SSAStepper())
    # push!(times, sol.t)

    # len_data = []
    # for (i, len) in enumerate(unique(sol_lengths))

    #     group_indices = findall(df.length .== len)

    #     group_sum = [sum([sol.u[t][idx] for idx in group_indices]) for t in 1:length(sol.t)]

    #     push!(len_data, group_sum)
    # end
    # push!(lengths, len_data)
    
    # assembly_index_data = []
    # for (i, ai) in enumerate(unique(sol_ai))

    #     group_indices = findall(df.ai .== ai)

    #     group_sum = [sum([sol.u[t][idx] for idx in group_indices]) for t in 1:length(sol.t)]

    #     push!(assembly_index_data, group_sum)
    # end
    # push!(ais, assembly_index_data)
    
    av = assembly_vector(sol, df)
    plot!(p_a, sol.t, av, alpha = 0.3, color = :blue);
    # push!(assemblies, av)

    target_of_interest = "BANANA"
    tv = target_vector(sol, species_strings, target_of_interest)
    append!(target_values, tv[end])
    plot!(p_t, sol.t, tv, alpha = 0.3, color = :red);
    # push!(targets, tv)
end

display(p_a)
# savefig(p_a, "AT_exploration/fig/multiple_runs/simple_example_assembly_50sims.png")
average = mean(target_values)
annotate!(p_t, 0.012, 45, text("Average number of target species formed: $average", 6, :red, :left));
display(p_t)
# savefig(p_t, "AT_exploration/fig/multiple_runs/simple_example_target_50sims.png")