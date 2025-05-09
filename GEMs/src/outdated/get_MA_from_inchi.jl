using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end

Pkg.activate(".")

using HTTP, JSON, JSON3, DataFrames, ProgressMeter

inchi_json = JSON.parsefile("../data/inchi_stuff/valid_inchi.json")
inchi_dict = Dict(k => v for d in inchi_json for (k, v) in d)
inchi_df = DataFrame(cpd = collect(keys(inchi_dict)), inchi = collect(values(inchi_dict)))

inchi_list = inchi_df.inchi

full_dictionary = Dict{String, Any}()

function batch_lookup_inchis(inchi_list::Vector{String}, batch_size::Int)
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
        sleep(1)

        response_data = JSON3.read(String(response.body))
        println("number of values returned: ", length(response_data), " out of ", batch_size)
        
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

batch_size = 10

@showprogress for batch in collect(Iterators.partition(inchi_list, batch_size))
    results = batch_lookup_inchis(Vector{String}(batch), batch_size)
    full_dictionary = merge(full_dictionary, results)
end

@info "Collected MA values for $(length(collect(values(full_dictionary)))) compounds."

open("../data/inchi_stuff/MA_values.json", "w") do f
    JSON.print(f, full_dictionary)
end
