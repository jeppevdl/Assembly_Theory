# Get assembly index of target `t`
function assembly_index(t)
    A = [Inf]
    _assembly_index([t], A)
    return only(A)
end

function _assembly_index(P, A)
    remnant = P[end]
    for (subLeft, subRight) in getsubstructures(remnant)
        if subRight == subLeft
            newP = [P; subLeft]
            newremnant = subRight
        end # else???

        curr_index = calc_index(newP)
        if curr_index < A[1]
            A .= curr_index
        end
        newP = [newP; newremnant]
        _assembly_index(newP, A)
    end
end

function getsubstructures(s) # ???
    return [[s[1:i], s[i+1:end]] for i in 1:(length(s)-1)]
end

function calc_index(s) # ???
    if length(s) <= 3
        return length(s) - 1
    else
        return Inf
    end
end


assembly_index("BANANA")
assembly_index("ABRACADABRA", return_all = true)
assembly_index("CADABRA", return_all = true)
# # tests
assembly_index("AAAA") == 2
assembly_index("BANANA") == 4

assembly_index("ABCABC")
assembly_index("redrum")
a, b, c = assembly_index("redrumredrum", return_all = true)
# assembly_index("redrumredrumredrumredrumredrumredrumredrumredrumredrumredrum") # still takes too long

# # timing

@time assembly_index("AAA")
@time assembly_index("BANANA", return_visited = true)
@time assembly_index("ABNNNBA")

#still not fully correct
@time assembly_index("ABCDEFXYZABCDEFMN", return_visited = true)