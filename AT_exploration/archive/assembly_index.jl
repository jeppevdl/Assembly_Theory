#from https://www.mdpi.com/1099-4300/24/7/884


#basic algorithm
function main(B, t, upper_bound)
    global A = length(B) + upper_bound
    assembly_index(B, t)
    return A - length(B)
end

function assembly_index(S, t)
    global A

    #use different stopping criterion?
    if A <= length(S) + 1
        return
    end
    
    for s1 in S
        for s2 in S
            combined = s1 * s2
            if combined == t && A > length(S) + 1
                A = length(S) + 1
                if A <= length(S) + 1
                    return
                end
            elseif combined in species
                assembly_index([S; combined], t)
            end
        end
    end
end


B = ["B", "A", "N"] #basic building blocks
# t = "ANNNBNN" #target string
t = rand(B, 8) |> join
upper_bound = length(t) #upper bound of complexity

#create list of species
species = ["B", "A", "N"] 
for i in 2:upper_bound 
    combinations = vec(map(collect, Iterators.product(ntuple(_ -> ["B", "A", "N"], i)...)))
    spec = [join(c) for c in combinations]
    species = vcat(species, spec)
end

#calculate assembly index
main(B, t, upper_bound)