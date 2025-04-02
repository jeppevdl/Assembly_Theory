using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end

Pkg.activate(".")

using HTTP, JSON, KEGGAPI, Graphs, GraphPlot, SimpleWeightedGraphs, GraphMakie, CairoMakie, NetworkLayout

reference_map = "map01310"
reactions = KEGGAPI.link("rn", reference_map).data[2]

# url = "https://rest.kegg.jp/link/rn/$reference_map"
# response = HTTP.get(url)
# reactions = [split(line, "\t")[2][4:end] for line in split(String(response.body), "\n") if !isempty(line)]

reaction_dict = Dict{String, Any}()
for reaction in reactions
    url2 = "https://rest.kegg.jp/get/$reaction"
    response2 = HTTP.get(url2)
    equation = [split(line, "EQUATION")[2][5:end] for line in split(String(response2.body), "\n") if startswith(line, "EQUATION")]
    reaction_dict[reaction] = equation
end

G = SimpleWeightedGraph(length(reaction_dict))
reaction_ids = collect(keys(reaction_dict))

for i in 1:length(reaction_ids)
    for j in i+1:length(reaction_ids)
        rxn_id1 = reaction_ids[i]
        rxn_id2 = reaction_ids[j]
        equation1 = reaction_dict[rxn_id1][1]
        equation2 = reaction_dict[rxn_id2][1]

        in1 = [compound.match for compound in eachmatch(r"C\d{5}", split(equation1, "<=>")[1]) if compound.match != "C00001"]
        out1 = [compound.match for compound in eachmatch(r"C\d{5}", split(equation1, "<=>")[2]) if compound.match != "C00001"]
        in2 = [compound.match for compound in eachmatch(r"C\d{5}", split(equation2, "<=>")[1]) if compound.match != "C00001"]
        out2 = [compound.match for compound in eachmatch(r"C\d{5}", split(equation2, "<=>")[2]) if compound.match != "C00001"]

        if length(intersect(in1, out2)) > 0 || length(intersect(out1, in2)) > 0 || length(intersect(in1, in2)) > 0 || length(intersect(out1, out2)) > 0
            add_edge!(G, i, j, 1)
        end
    end
end

fig, ax, plt = graphplot(G, layout=Stress(), node_size=10, edge_width=0.5, nlabels = reaction_ids);
resize!(fig, 1000, 1000)
fig