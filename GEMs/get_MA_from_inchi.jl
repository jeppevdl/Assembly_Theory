using Pkg
if basename(pwd()) != "GEMs"
    cd("GEMs")
end

Pkg.activate(".")

using HTTP, JSON, JSON3, DataFrames

inchi_dict = Dict(JSON.parsefile("data/inchi_dict_subset.json"))
inchi_df = DataFrame(cpd = collect(keys(inchi_dict)), inchi = collect(values(inchi_dict)))

inchi_list = inchi_df.inchi

full_dictionary = Dict{String, Any}()

function batch_lookup_inchis(inchi_list::Vector{Any}, batch_size::Int)
    results = Dict{String, Any}()
    base_url = "http://molecular-assembly.com/batch_lookup"

    for batch in Iterators.partition(inchi_list, batch_size)
        dict_inchis = Dict{String, Any}("n" => length(batch))
        
        for (i, inchi) in enumerate(batch)
            j = i - 1
            dict_inchis["i$j"] = inchi
        end

        println("Query parameters: ", dict_inchis)

        response = HTTP.get(base_url, query = dict_inchis)
        response_data = JSON3.read(String(response.body))
        for entry in response_data
            inchi_id = entry["inchi"]
            cpd_id = inchi_df.cpd[inchi_df.inchi .== inchi_id][1]
            ma = entry["MA"]
            method = entry["method"]
            results[cpd_id] = Dict("MA" => ma, "method" => method, "inchi" => inchi_id)
        end
    end

    return results
end

batch_size = 5

for i in 1:batch_size:length(inchi_list)
    results = batch_lookup_inchis(inchi_list[i:min(i + batch_size - 1, length(inchi_list))], batch_size)
    full_dictionary = merge(full_dictionary, results)
end

open("data/MA_dict.json", "w") do f
    JSON.print(f, full_dictionary)
end
