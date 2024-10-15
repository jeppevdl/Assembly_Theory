using Catalyst, OrdinaryDiffEq, Plots, Graphs, GraphRecipes, JumpProcesses, DataFrames

t = default_t() # onafhankelijke variabele

len_str = 7 # definieer maximum lengte van gecombineerde string
species = ["A(t)", "B(t)"] .|> Meta.parse

for i in 2:len_str 
    combinations = vec(map(collect, Iterators.product(ntuple(_ -> ["A", "B"], i)...)))
    spec = [join(c) * "(t)" for c in combinations] .|> Meta.parse
    species = vcat(species, spec)
end

species = [only(@eval @species $(element)) for element in species]

reactions = []
for (i, str1) in enumerate(species)
    for (j, str2) in enumerate(species)

        comb = replace(string(str1), "(t)" => "") * string(str2) .|> Meta.parse
        comb = only(@eval @species $(comb))

        if string(str1) == string(str2) && length(string(comb)) <= len_str + 3 # limiet op mogelijke combinaties omdat species op voorhand gedefinieerd moeten zijn, 
                                                                               # hier ook +3 omdat comb nog "(t)" bevat (= 3 karakters)
            push!(reactions, Reaction(0.5, [str1], [comb], [2], [1]))
        elseif length(string(comb)) <= len_str + 3
            push!(reactions, Reaction(0.5, [str1, str2], [comb], [1, 1], [1]))
        end

    end
end

species_symbols = [Symbol(replace(string(s), "(t)" => "")) for s in species]
@named jumpmodel = ReactionSystem(reactions, t)
u0 = vcat([:A => 500, :B => 500], [species_symbols[i] => 0 for i in 3:length(species_symbols)])
tspan = (0.0, 1.0)
ps = []

prob = DiscreteProblem(complete(jumpmodel), u0, tspan, ps)

jump_prob = JumpProblem(complete(jumpmodel), prob, Direct())

sol = solve(jump_prob, SSAStepper())

# plot(sol)

sol_lengths = map(length, map(string, species_symbols))
df = DataFrame(name = species_symbols, length = sol_lengths)

p = plot(title = "Sum of Variables by Length", xlabel = "Time", ylabel = "Sum");

for len in unique(sol_lengths)

    group_indices = findall(df.length .== len)

    group_sum = [sum([sol.u[t][idx] for idx in group_indices]) for t in 1:length(sol.u)]

    plot!(p, 1:length(sol.u), group_sum, label = "Length $len")
end

display(p)
# savefig(p, "AT_abundance7.png")


#--------------------------------------------------------------------------------
# Create the directed graph
g = DiGraph(length(species_symbols))

str_species = [string(s) for s in species]
# Add edges based on the reactions
for r in reactions
    for reactant in r.substrates
        for product in r.products
            # Get the index of the reactant and product in the species_symbols array
            reactant_idx = findfirst(x -> x == string(reactant), str_species)
            product_idx = findfirst(x -> x == string(product), str_species)
            
            # Ensure both reactant and product are found
            if reactant_idx !== nothing && product_idx !== nothing
                add_edge!(g, reactant_idx, product_idx)
            else
                println("Warning: Could not find species for reaction involving $reactant and $product")
            end
        end
    end
end

# Plot the graph
graphplot(g, names=species_symbols, method=:tree, 
    nodeshape=:rect, nodesize=0.2, nodecolor=:lightblue)