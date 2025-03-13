using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end

Pkg.activate(".")

using HTTP, JSON

reference_map = "map00620"
url = "https://rest.kegg.jp/link/rn/$reference_map"
response = HTTP.get(url)

reactions = [split(line, "\t")[2][4:end] for line in split(String(response.body), "\n") if !isempty(line)]

reaction_dict = Dict{String, Any}()
for reaction in reactions
    url2 = "https://rest.kegg.jp/get/$reaction"
    response2 = HTTP.get(url2)
    equation = [split(line, "EQUATION")[2][5:end] for line in split(String(response2.body), "\n") if startswith(line, "EQUATION")]
    reaction_dict[reaction] = equation
end

MA_dict = JSON.parsefile("data/MA_dict.json")

for reaction in reactions
    equation = reaction_dict[reaction][1]
    ids = [m.match for m in eachmatch(r"C\d{5}", equation)]
    println("\n")
    println(reaction)
    println(equation)
    for id in ids
        if haskey(MA_dict, id)
            println(id, " ", MA_dict[id]["MA"])
        end
    end
end
