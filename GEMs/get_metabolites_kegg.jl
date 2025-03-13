using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end

Pkg.activate(".")

using HTTP, JSON

inchi_dict = Dict{String, String}()
not_found = Dict{String, Int}()
function get_compound_inchi(pathway_id, dict::Dict{String, String}, not_found::Dict{String, Int})
    url = "https://rest.kegg.jp/link/cpd/$pathway_id"
    response = HTTP.get(url)
    compounds = [split(line, "\t")[2][5:end] for line in split(String(response.body), "\n") if !isempty(line) && !haskey(dict, split(line, "\t")[2][5:end])]

    batch_size = 10
    for batch in Iterators.partition(compounds, batch_size)
        compound_list = join(batch, "+")
        url2 = "https://rest.kegg.jp/conv/chebi/compound:$compound_list"
        response2 = HTTP.get(url2)

        for line in split(String(response2.body), "\n")
            if !isempty(line)
                parts = split(line, "\t")
                compound = parts[1][5:end]
                chebi_id = parts[2][7:end]

                url_inchi = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId=CHEBI:$chebi_id"
                response_inchi = HTTP.get(url_inchi)
                inchi_match = match(r"<inchi>(.*?)</inchi>", String(response_inchi.body))

                if !isnothing(inchi_match)
                    inchi = inchi_match.captures[1]
                    dict[compound] = inchi
                else
                    println("No InChI found for $compound")
                    not_found[compound] = get(not_found, compound, 0) + 1
                end
            end
        end
    end
    return dict, not_found
end

reference_pathways = HTTP.get("https://rest.kegg.jp/list/pathway")
reference_pathway_ids = [split(line, "\t")[1] for line in split(String(reference_pathways.body), "\n") if !isempty(line)]

for pathway_id in reference_pathway_ids[1:100]
    inchi_dict, not_found = get_compound_inchi(pathway_id, inchi_dict, not_found)
end

open("data/inchi_dict_subset.json", "w") do f
    JSON.print(f, inchi_dict)
end

open("data/not_found_subset.json", "w") do f
    JSON.print(f, not_found)
end


n = 30
for i in 1:n

end