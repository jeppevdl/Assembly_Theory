# Calculates assembly index `A` of target `t` starting from elementary building blocks `B`
function assembly_index(B, t)
    upper_bound = length(t)
    A = [length(B) + upper_bound]
    _assembly_index!(B, t, A)
    return A[1] - length(B)
end

# Recursive function which sets assembly index to lowest encountered
function _assembly_index!(S, t, A)

    if length(S) + 1 >= A[1]
        return
    end

    for s1 in S
        for s2 in S
 
            combined = s1*s2
            if !(combined in S) # combination not already explored
                if combines_to_target(s1, s2, t) && (length(S) + 1) < A[1]
                   A .= length(S) + 1
                elseif combines_to_something(s1, s2, t)
                    _assembly_index!([S; combined], t, A)
                end
            end

        end
    end

end

combines_to_target = (s1, s2, t) -> s1*s2 ==  t
combines_to_something = (s1, s2, t) -> (length(s1)+length(s2)) < length(t)

# basic building blocks
B = ["B", "A", "N"]

# # tests
assembly_index(B, "AAAA") == 2
assembly_index(B, "BANANA") == 4

# # timing

@time assembly_index(B, "AAA")
@time assembly_index(B, "BANANA")
@time assembly_index(B, "ABNNNBA")
# @time assembly_index(B, "ABRACADABRA") # only V3 can handle