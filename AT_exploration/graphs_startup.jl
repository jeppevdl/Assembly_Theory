using Pkg; Pkg.activate(".")
using Graphs
using GraphRecipes, Plots

num_C = 5
g = SimpleDiGraph(num_C) # Gerichte graaf met `num_C` nodes
for curr_Cs in 1:num_C
    for leq_Cs in 1:(curr_Cs-1)
        add_edge!(g, leq_Cs, curr_Cs) # voegt een edge toe bij graaf `g` tussen 2 gespecifieerde nodes 
    end
end

graphplot(g, names = ["C"^i for i in 1:num_C])