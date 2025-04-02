using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end

Pkg.activate(".")
using COBREXA, Tulip, Escher, CairoMakie, Colors

# open the SBML file and load the contents
model = load_model("data/e_coli_core.xml")

# visualize the metabolism 
escherplot(
    "data/e_coli_core_map.json"; 
    reaction_show_text = true,
    reaction_edge_color = :grey,
    metabolite_show_text = true,
    metabolite_node_colors = Dict("glc__D_e" => :red),
    metabolite_node_color = :lightskyblue,
)
hidexdecorations!(current_axis())
hideydecorations!(current_axis())
current_figure()

# perform flux balance analysis 
fluxes = flux_balance_analysis_dict(model, Tulip.Optimizer)

tol = 1e-3 # threshold if reaction is active

# ids from e coli model are different in e coli map so the red fluxes don't show
escherplot(
    "data/e_coli_core_map.json";
    reaction_edge_colors = Dict(id => :red for (id, flux) in fluxes if abs(flux) > tol),
)
hidexdecorations!(current_axis())
hideydecorations!(current_axis())
current_figure()

# convert to a model type that is efficient to modify
m = convert(StandardModel, model)

# find the model objective value if oxygen or carbon dioxide transports are disabled
COBREXA.screen(m, # the base model
    variants=[ # this specifies how to generate the desired model variants
        [], # one with no modifications, i.e. the base case
        [with_changed_bound("R_O2t", lower=0.0, upper=0.0)], # disable oxygen
        [with_changed_bound("R_CO2t", lower=0.0, upper=0.0)], # disable CO2
        [with_changed_bound("R_O2t", lower=0.0, upper=0.0),
	        with_changed_bound("R_CO2t", lower=0.0, upper=0.0)], # disable both
    ],
    # this specifies what to do with the model variants (received as the argument `x`)
    analysis = x ->
        flux_balance_analysis_dict(x, Tulip.Optimizer)["R_BIOMASS_Ecoli_core_w_GAM"],
)

# load the task distribution package, add several worker nodes, and load
# COBREXA and the solver on the nodes
using Distributed
addprocs(4)
@everywhere using COBREXA, Tulip

# get a list of the workers
worker_list = workers()

# run the processing in parallel for many model variants
res = screen(m,
    variants=[
	# create one variant for each reaction in the model, with that reaction knocked out
        [with_changed_bound(reaction_id, lower=0.0, upper=0.0)]
	for reaction_id in reactions(m)
    ],
    analysis = model -> begin
	# we need to check if the optimizer even found a feasible solution,
	# which may not be the case if we knock out important reactions
    	sol = flux_balance_analysis_dict(model, Tulip.Optimizer)
	isnothing(sol) ? nothing : sol["R_BIOMASS_Ecoli_core_w_GAM"]
    end,
    # run the screening in parallel on all workers in the list
    workers = worker_list,
)

Dict(reactions(m) .=> res)