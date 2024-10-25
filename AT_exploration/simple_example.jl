using Catalyst, OrdinaryDiffEq, Plots, Graphs, GraphRecipes, JumpProcesses, DataFrames
#functions--------------------------------------------------------------------------------

#function to create a reaction network
function create_reactions(building_blocks, upper_bound, rate)
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
                if string(str1) == string(str2) 
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
#-----------------------------------------------------------------------------------------
# setting up model

building_blocks = ["A", "B", "N"]
upper_bound = 8
rate = 1
t = default_t()

#creating the reactions can take a long time
reactions, species = create_reactions(building_blocks, upper_bound, rate);
@named jumpmodel = ReactionSystem(reactions, t)

u0 = vcat([Symbol(replace(string(sp), "(t)" => "")) => 1000 for sp in species[1:length(building_blocks)]], [Symbol(replace(string(sp), "(t)" => "")) => 0 for sp in species[length(building_blocks)+1:end]])
tspan = (0.0, 1.0)
ps = []

prob = DiscreteProblem(complete(jumpmodel), u0, tspan, ps)

jump_prob = JumpProblem(complete(jumpmodel), prob, Direct())

sol = solve(jump_prob, SSAStepper())

#-----------------------------------------------------------------------------------------
#create dataframe with names, lengths and assembly indices for plotting
species_symbols = [Symbol(replace(string(sp), "(t)" => "")) for sp in jumpmodel.unknowns]
sol_lengths = map(length, map(string, species_symbols))
sol_ai = map(assembly_index, map(string, species_symbols))

df = DataFrame(name = species_symbols, length = sol_lengths, ai = sol_ai)

#Plot by length
p = plot(title = "Sum of Variables by Length", xlabel = "Timestep", ylabel = "Sum", dpi = 600);

for len in unique(sol_lengths)

    group_indices = findall(df.length .== len)

    group_sum = [sum([sol.u[t][idx] for idx in group_indices]) for t in 1:length(sol.t)]

    plot!(p, 1:length(sol.u), group_sum, label = "Length $len")
end

display(p)
savefig(p, "AT_exploration/fig/simple_example_8.png")

#Plot by assembly index
p_ai = plot(title = "Sum of Variables by Assembly Index", xlabel = "Timestep", ylabel = "Sum", dpi = 600);

for ai in unique(sol_ai)

    group_indices = findall(df.ai .== ai)

    group_sum = [sum([sol.u[t][idx] for idx in group_indices]) for t in 1:length(sol.t)]

    plot!(p_ai, 1:length(sol.u), group_sum, label = "Assembly Index $ai")
end

display(p_ai)
savefig(p_ai, "AT_exploration/fig/simple_example__ai.png")

#Plot assembly through time
assembly_vector = []
for i in 1:length(sol.t)
    NT = sum(sol.u[i])
    A = 0
    for (j, ai) in enumerate(df.ai)
        if sol.u[i][j] != 0
            a = ℯ^(ai) * (sol.u[i][j]-1) / NT
            A += a
        end
    end
    append!(assembly_vector, A)
end

p_a = plot(title = "Assembly through time", xlabel = "Timestep", ylabel = "Assembly", dpi = 600, legend=false);
plot!(p_a, 1:length(assembly_vector), assembly_vector);
display(p_a)
savefig(p_a, "AT_exploration/fig/simple_example__assembly_8.png")