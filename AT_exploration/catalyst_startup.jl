# using Pkg; Pkg.activate(".")
using Catalyst, OrdinaryDiffEq, Plots

t = default_t() # onafhankelijke variabele

num_C = 15
C_vec = ["C"^i * "(t)" for i in 1:num_C] .|> Meta.parse # vector met species namen maken en dan omzetten naar Expr (nodig voor macro shenanigans)

species = [only(@eval @species $(C_species)) for C_species in C_vec] # macro magic

reactions = []
for i in 1:num_C
    for j in i:num_C-i
        if i == j
            push!(reactions, Reaction(1, [species[i]], [species[2*i]], [2], [1]))
        else
            push!(reactions, Reaction(1, [species[i], species[j]], [species[i+j]], [1, 1], [1]))
        end
    end
end

reactions

# klassieke Catalyst dingen
@named model = ReactionSystem(reactions, t)
u0 = vcat([:C => 100.0], [Symbol("C"^i) => 0.0 for i in 2:num_C])
# u0 = [:C => 100.0, :CC => 50.0, :CCC => 25.0, :CCCC => 12.5, :CCCCC => 0.0]
tspan = (0., 1.)
ps = []
ode = ODEProblem(complete(model), u0, tspan, ps)

# RUN IT!
sol = solve(ode)
plot(sol)