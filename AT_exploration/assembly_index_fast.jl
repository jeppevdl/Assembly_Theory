struct Assembly
    val::Int64
    visited::Vector{String}
end

# Get assembly index of target `t`
function assembly_index(t; return_all = false) # aka "Main"
    assemblies = Dict{String, Assembly}()
    A, visited = _assembly_index(t, String[], assemblies)
    if return_all
        return A, visited, assemblies
    else
        return A
    end
end

# Function calculating the assembly index of a compound `s`
function _assembly_index(s, visited, assemblies) # aka "AssemblyIndex"
    if s in visited
        return 0, visited
    elseif haskey(assemblies, s)
        return assemblies[s].val, [visited; s]
    elseif length(s) <= 3
        # strings of length l <= 3 can only be constructed in l-1 combinations
        return length(s) - 1, [visited; s]
    else
        # not yet visited "complex" compound => split into substructures
        substructure_pairs = getsubstructures(s)

        As = Vector{Int64}(undef, length(substructure_pairs))
        visiteds = Vector{Vector{String}}(undef, length(substructure_pairs))
        
        for (sp_idx, sp) in enumerate(substructure_pairs)

            if occursin(sp[1], sp[2]) && !(sp[1] in get(assemblies, sp[2], Assembly(0, String[])).visited) # abra in cadabra but A(cadabra) was calculated without abra
                delete!(assemblies, sp[2]) # we need to recalculate cadabra
                    #! should be only if you can do better than just doing abra + cadabra with current values
                A1, visited1 = _assembly_index(sp[1], visited, assemblies) # calc abra as normal
                A2, visited2 = _assembly_index(sp[2], [visited1; sp[1]], assemblies) # calc cadabra including abra as existing building block

                As[sp_idx] = A1 + A2 + 1
                visiteds[sp_idx] = visited2
            elseif occursin(sp[2], sp[1]) && !(sp[2] in get(assemblies, sp[1], Assembly(0, String[])).visited) # reverse situation of previous
                delete!(assemblies, sp[1])
                A2, visited2 = _assembly_index(sp[2], visited, assemblies)
                A1, visited1 = _assembly_index(sp[1], [visited2; sp[2]], assemblies)

                As[sp_idx] = A1 + A2 + 1
                visiteds[sp_idx] = visited1
            else
                A1, visited1 = _assembly_index(sp[1], visited, assemblies)
                A2, visited2 = _assembly_index(sp[2], visited, assemblies)

                As[sp_idx] = A1 + A2 + 1
                visiteds[sp_idx] = unique([visited1; visited2])
            end
        end

        min_idx = argmin(As)
        assemblies[s] = Assembly(As[min_idx], visiteds[min_idx])
        return As[min_idx], visiteds[min_idx]
    end
end

function getsubstructures(s)
    return [[s[1:i], s[i+1:end]] for i in 1:(length(s)-1)]
end

assembly_index("BANANA", return_all = true)
assembly_index("ABRACADABRA", return_all = true)
assembly_index("CADABRA", return_all = true)

# # tests

assembly_index("AAAA") == 2
assembly_index("BANANA") == 4
assembly_index("ABCDEFXYZABCDEFMN") == 11
assembly_index("redrum") == 5
assembly_index("redrumredrum") == 6


# # timing

@time assembly_index("AAA")
@time assembly_index("BANANA")
@time assembly_index("ABNNNBA")
@time assembly_index("BANANANANANANA")
# @time assembly_index("redrumredrumredrumredrumredrumredrumredrumredrumredrumredrum") # still takes too long