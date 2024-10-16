function main(B, t)
    upper_bound = length(t)
    global A = length(B) + upper_bound
    assembly_index(B, t)

    return A - length(B)
end

function assembly_index(S, t)
    global A

    if A <= length(S) + 1
        return
    end
    
    for s1 in S
        for s2 in S
            combined = s1 * s2
            if combined == t && A > length(S) + 1
                A = length(S) + 1
                return

            elseif length(combined) < length(t)
                assembly_index([S; combined], t)
            end
        end
    end
end

B = ["B", "A", "N"] #basic building blocks

# # tests
main(B, "AAAA") == 2
main(B, "BANANA") == 4

# # timing

@time main(B, "AAA")
@time main(B, "BANANA")
@time main(B, "ABNNNBA")
# @time main(B, "ABRACADABRA") # only V3 can handle