using JumpProcesses, Catalyst, OrdinaryDiffEq, Plots, StochasticDiffEq

t = default_t() # onafhankelijke variabele

num_C = 7
C_vec = ["C"^i * "(t)" for i in 1:num_C] .|> Meta.parse # vector met species namen maken en dan omzetten naar Expr (nodig voor macro shenanigans)

species = [only(@eval @species $(C_species)) for C_species in C_vec] # macro magic

reactions = []
for i in 1:num_C
    for j in i:num_C-i
        if i == j
            push!(reactions, Reaction(0.5, [species[i]], [species[2*i]], [2], [1]))
        else
            push!(reactions, Reaction(1, [species[i], species[j]], [species[i+j]], [1, 1], [1]))
        end
    end
end

# klassieke Catalyst dingen
@named jumpmodel = ReactionSystem(reactions, t)

u0 = vcat([:C => 500], [Symbol("C"^i) => 0 for i in 2:num_C])
tspan = (0.0, 1.0)
ps = []

prob = DiscreteProblem(complete(jumpmodel), u0, tspan, ps)

jump_prob = JumpProblem(complete(jumpmodel), prob, Direct())

sol = solve(jump_prob, SSAStepper())

plot(sol)