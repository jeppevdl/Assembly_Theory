# Get assembly index of target `t`
function assembly_index(t; return_visited = false) # aka "Main"
    A, visited = _assembly_index(t, String[])
    if return_visited
        return A, visited
    else
        return A
    end
end

# Function calculating the assembly index of a compound `s`
function _assembly_index(s, visited) # aka "AssemblyIndex"
    if s in visited
        # no construction cost for previously visited compounds
        return 0, visited 
    elseif length(s) <= 3
        # strings of length <= 3 can only be constructed in a set amount of combinations
        return length(s) - 1, [visited; s]
    else
        # not yet visited "complex" compound => split into substructures
        substructure_pairs = getsubstructures(s)

        As = Vector{Int64}(undef, length(substructure_pairs))
        visiteds = Vector{Vector{String}}(undef, length(substructure_pairs))

        for (sp_idx, sp) in enumerate(substructure_pairs)
            A1, visited1 = _assembly_index(sp[1], visited)
            A2, visited2 = _assembly_index(sp[2], [visited1; sp[1]])
            As[sp_idx] = A1 + A2 + 1
            visiteds[sp_idx] = [visited2; sp[2]]
        end

        min_idx = argmin(As)
        return As[min_idx], visiteds[min_idx]
    end
end

function getsubstructures(s)
    return [[s[1:i], s[i+1:end]] for i in 1:(length(s)-1)]
end

assembly_index("BANANA", return_visited = true)
assembly_index("ABRACADABRA", return_visited = true)

# # tests
assembly_index("AAAA") == 2
assembly_index("BANANA") == 4
# assembly_index("redrumredrumredrumredrumredrumredrumredrumredrumredrumredrum") == 9 # still takes too long

# # timing

@time assembly_index("AAA")
@time assembly_index("BANANA")
@time assembly_index("ABNNNBA")
@time assembly_index("ABRACADABRA")