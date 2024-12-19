using DynamicPolynomials
using Combinatorics
using LinearAlgebra
import Base.string

@polyvar A0 A1 A2 A3 A4 A12 A13 A14 A23 A24 A34 A123 A124 A134 A234 A1234

names =["","1","2","3","4","12","13","14","23","24","34","123","124","134","234","1234"] 
vars = [A0,A1,A2,A3,A4,A12,A13,A14,A23,A24,A34,A123,A124,A134,A234,A1234]
subsets = vcat([[]], [parse.(Int64, split(name, "")) for name in names[2:end]])
mapping = Dict(zip(names, vars))

hyperdeterminant = A3^2*A12^2 - 2*A2*A3*A12*A13 + A2^2*A13^2 - 2*A1*A3*A12*A23 - 2*A1*A2*A13*A23 + 4*A1*A2*A3*A123 + A1^2*A23^2 + 4*A0*A12*A13*A23 - 2*A0*A3*A12*A123 - 2*A0*A2*A13*A123 - 2*A0*A1*A23*A123 + A0^2*A123^2


# poly1 =
# A1234^2*A1*A0 - A1234*A123*A14*A0 - A1234*A123*A1*A4 - A1234*A124*A13*A0 - A1234*A124*A1*A3 -
# A1234*A134*A12*A0 - A1234*A134*A1*A2 - A1234*A12*A34*A1 - A1234*A13*A24*A1 - A1234*A14*A23*A1 -
# A123*A124*A13*A4 - A123*A124*A14*A3 - A123*A134*A12*A4 - A123*A134*A14*A2 - A123*A234*A14*A1 -
# A123*A12*A14*A34 - A123*A13*A14*A24 - A124*A134*A12*A3 - A124*A134*A13*A2 - A124*A234*A13*A1 - 
# A124*A12*A13*A34 - A124*A13*A14*A23 - A134*A234*A12*A1 - A134*A12*A13*A24 - A134*A12*A14*A23 + 
# 2*A1234*A12*A13*A4 + 2*A1234*A12*A14*A3 + 2*A1234*A13*A14*A2 + 2*A123*A124*A134*A0 + 2*A123*A124*A34*A1 +
# 2*A123*A134*A24*A1 + 2*A124*A134*A23*A1 + 2*A234*A12*A13*A14 + A1234*A234*A1^2 + A123^2*A14*A4 +
# A123*A14^2*A23 + A124^2*A13*A3 + A124*A13^2*A24 + A134^2*A12*A2 + A134*A12^2*A34
# 
# poly2 = (A0^2*A1234^2 + A4^2*A123^2 + A3^2*A124^2 + A34^2*A12^2 + A2^2*A134^2 + A24^2*A13^2 + A23^2*A14^2 + A234^2*A1^2) + 2*(- A0*A34*A12*A1234 - A0*A24*A13*A1234 - A0*A23*A14*A1234 - A0*A234*A1*A1234 + A0*A234*A14*A123 + A0*A234*A13*A124 + A0*A234*A134*A12 - A4*A3*A124*A123 - A4*A2*A134*A123 + A4*A23*A1*A1234 - A4*A23*A14*A123 + A4*A23*A13*A124 + A4*A23*A134*A12 - A4*A234*A1*A123 - A3*A2*A134*A124 + A3*A24*A1*A1234 + A3*A24*A14*A123 - A3*A24*A13*A124 + A3*A24*A134*A12 - A3*A234*A1*A124 + A34*A2*A1*A1234 + A34*A2*A14*A123 + A34*A2*A13*A124 - A34*A2*A134*A12 - A34*A24*A13*A12 - A34*A23*A14*A12 - A2*A234*A1*A134 - A24*A23*A14*A13)
# 
# 
# # Test
# # display(subs(hyperdeterminant, get.([mapping], names, "na") => minors))
# # display(subs(poly1, get.([mapping], names, "na") => minors))
# # display(subs(poly2, get.([mapping], names, "na") => minors))
# # exit()
# #
bools = [false, true]
polys = []
for p in permutations([1,2,3,4])
    for z in Iterators.product(bools,bools,bools,bools)
        permutation = []
        for (name, subset) = zip(names, subsets)
            permuted_subset = []
            for (i, zi) = enumerate(z)
                if (p[i] in subset) ⊻ zi
                    push!(permuted_subset, i)
                end
            end
            push!(permutation, join(sort(permuted_subset), ""))
        end
        push!(polys, subs(hyperdeterminant, vars => get.([mapping], permutation, "na")))
        # push!(polys, subs(poly1, vars => get.([mapping], permutation, "na")))
        # push!(polys, subs(poly2, vars => get.([mapping], permutation, "na")))
    end
end

polys = unique(polys)

function string(m::DynamicPolynomials.Monomial)
    parts = []
    for (v, i) = zip(m.vars, m.z)
        push!(parts, string(v) + "^" + string(i))
    end
    return "*".join(parts)
end

file = open("hyperdeterminants.txt","a")

@polyvar c11 c21 c31 c41 c12 c22 c32 c42
for poly in polys
    write(file, string(subs(poly, [A14,A24,A34,A23] => [-A12 - A13 - c12 + c22 + c32 + c42, A13 + c12 - c22 + c32 - c42, A12 + c12 + c22 - c32 - c42, -A12 - A13 + 2*c42])))
    write(file, ",\n")
end

# using JLD2
# using FileIO
# 
# # Create a file
# file = File(format"JLD2", "relations.jld2")
# 
# # Save data into the file
# save(file, "relations", polys)
# 
# # Load the file
# data = load(file)
# 
# # Display user-visible data
# display(data["relations"][1])
