### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ d84119d6-b81c-4b69-bf12-c2ce24d1d321
# ╠═╡ disabled = true
#=╠═╡
begin
	using Pkg
	Pkg.activate("/Users/michielstock/SVN_GITHUB/Teaching/ModSim")
end
  ╠═╡ =#

# ╔═╡ 5b89018a-b969-11ef-2972-d5cc70066b6b
using Plots, Catalyst, SparseArrays

# ╔═╡ 7711e87c-af4b-4a11-857e-194895a1cb48
using DifferentialEquations, ModelingToolkit

# ╔═╡ 9c3cffed-5cda-4137-bf6c-40f5dbdc1a8c
building_blocks = ["A", "B"]

# ╔═╡ 872860c4-0a5a-4998-88c0-5c1027f20440
# can these products merge?
rule(S1, S2) = count(!=('A'), S1) + count(!=('A'), S2) ≤ 2

# ╔═╡ 3a506f1e-d327-466e-aedb-8dbbf6cc20e2
n_assembly_steps = 4

# ╔═╡ d9a3c06d-e74e-4b4b-a7ee-f81f19b310b1
begin
	assembly_indices = Dict(bb=>0 for bb in building_blocks)
	r = (S1="S2", S2="S2", P="P")
	reactions = Dict{Tuple{String,String},String}()
	for t in 1:n_assembly_steps
		assemblies = collect(keys(assembly_indices))
		for (s1, s2) in Iterators.product(assemblies, assemblies)
			if rule(s1, s2)
				p = s1 * s2
				!haskey(assembly_indices, p) && (assembly_indices[p] = t)
				reactions[(s1, s2)] = p
			end
		end
	end

	assemblies = keys(assembly_indices) |> collect
	sort!(assemblies, by=length)
	species_index = Dict(s=>i for (i, s) in enumerate(assemblies))

	
	n_reactions = length(reactions)
	n_species = length(assembly_indices)
end;

# ╔═╡ 1db71ccd-1441-4a8a-be6d-b1ae682728d3
assembly_indices

# ╔═╡ 15be741e-213f-4432-a6c6-68e79d8f45b3
length(assemblies)

# ╔═╡ 95fc190d-7440-410c-a4a7-3f7669254ec5
let
	scatter([assembly_indices[a] for a in assemblies], length.(assemblies), xlab="assembly index", ylab="length")
end

# ╔═╡ 6a245117-97e2-4575-ac10-185d747150cc
reactions

# ╔═╡ 17c37a1c-dcaf-4613-8d54-ee3048670bd7
length(reactions)

# ╔═╡ 4b69faef-b885-4011-bed9-3c94e8154511
k = 1  # reaction rate

# ╔═╡ c24df670-2012-4b25-ae1c-9b20479cb10a
begin
	S = spzeros(Int, n_species, n_reactions)
	for (i, ((s1, s2), p)) in enumerate(reactions)
		si, sj = species_index[s1], species_index[s2]
		S[si, i] += -1
		S[sj, i] += -1
		pi = species_index[p]
		S[pi, i] = 1
	end
	S
end

# ╔═╡ 37580803-39ae-45c7-9c27-faad1edfd481
histogram(collect(values(assembly_indices)), yscale=:log10, xlab="assembly index")

# ╔═╡ 1e731bb0-b7d0-4e6d-a0cf-5e549068ce53
assembly_equation(copies) = sum(exp(assembly_indices[spi]) * (ni-1) for (spi, ni) in copies) / sum(values(copies))

# ╔═╡ 59afc414-b771-4b56-8f09-9f734073a058
assembly_indices

# ╔═╡ a5ee13e9-11ba-43bb-893d-3e2ae9d0cd87
copies = Dict{String,Int}()

# ╔═╡ 0e5ed6c6-8293-464b-bc5b-c7b76ce7e421
function simulate!(copies, n_steps=1; bb0=2)
	for t in 1:n_steps
		for bb in building_blocks
			copies[bb] = bb0
		end
		# need to sample according to copy number
		spi, ni = rand(copies)
		spj, nj = rand(copies)
		spi != spj || ni ≥ 2 || continue
		if haskey(reactions, (spi, spj))
			p = reactions[(spi, spj)]
			copies[p] = get(copies, p, 0) + 1
			copies[spi] -= 1
			copies[spj] -= 1
			copies[spi] == 0 && delete!(copies, spi)
			copies[spj] == 0 && delete!(copies, spj)
		end
	end
	return copies
end

# ╔═╡ daa90cf6-7684-4b74-b6c6-6049cffede30
#simulate!(copies, 10000)

# ╔═╡ 762e14c1-3be1-4ba6-8236-cc417f1edcd1
begin
	# formation reactions
	rxs = [eval(Meta.parse("@reaction k, $si + $sj --> $p")) for ((si, sj), p) in reactions]
	# constant addition of building blocks
	for bb in building_blocks
		push!(rxs, eval(Meta.parse("@reaction a, 0 --> $bb")))
	end
	# decay
	for s in assemblies
		push!(rxs, eval(Meta.parse("@reaction d, $s --> 0")))
	end
end

# ╔═╡ 38d8f63d-c7cd-41cb-8ff3-21070a12b546
t = default_t()

# ╔═╡ 9b90ebcc-f29f-424c-bc0a-821d1687a16e
begin
	@named rn = ReactionSystem(rxs, t);
	rn = complete(rn)
end;

# ╔═╡ 9f91ebb7-a492-44b2-a549-1f446c46ae34
# ╠═╡ disabled = true
#=╠═╡
g = Graph(rn)  # dragons, only execute when small
  ╠═╡ =#

# ╔═╡ 9aa3dc11-afc8-4b28-8146-0963ee0e5c5c
oprob = ODEProblem(rn, zeros(n_species), (0, 1000.), [:k=>1, :a=>1, :d=>0.1])

# ╔═╡ 2ecbaf82-d534-430b-a4b2-a7974ad7ce3b
sol = solve(oprob, Tsit5())

# ╔═╡ c3b20185-7ef0-4602-861d-5418bfb8a98f
histogram(sol[end], ylims = (0, 10))

# ╔═╡ 01628c9a-3c07-43bd-aa29-7db146633034
end_concentrations = sort!([(c, string(s)) for (c, s) in zip(sol[end], species(rn))], rev=true)

# ╔═╡ 62f7e767-3806-47f9-baa8-43c033528d07
assembly_indices

# ╔═╡ Cell order:
# ╠═5b89018a-b969-11ef-2972-d5cc70066b6b
# ╠═d84119d6-b81c-4b69-bf12-c2ce24d1d321
# ╠═7711e87c-af4b-4a11-857e-194895a1cb48
# ╠═9c3cffed-5cda-4137-bf6c-40f5dbdc1a8c
# ╠═872860c4-0a5a-4998-88c0-5c1027f20440
# ╠═3a506f1e-d327-466e-aedb-8dbbf6cc20e2
# ╠═d9a3c06d-e74e-4b4b-a7ee-f81f19b310b1
# ╠═1db71ccd-1441-4a8a-be6d-b1ae682728d3
# ╠═15be741e-213f-4432-a6c6-68e79d8f45b3
# ╠═95fc190d-7440-410c-a4a7-3f7669254ec5
# ╠═6a245117-97e2-4575-ac10-185d747150cc
# ╠═17c37a1c-dcaf-4613-8d54-ee3048670bd7
# ╠═4b69faef-b885-4011-bed9-3c94e8154511
# ╠═c24df670-2012-4b25-ae1c-9b20479cb10a
# ╠═37580803-39ae-45c7-9c27-faad1edfd481
# ╠═1e731bb0-b7d0-4e6d-a0cf-5e549068ce53
# ╠═59afc414-b771-4b56-8f09-9f734073a058
# ╠═a5ee13e9-11ba-43bb-893d-3e2ae9d0cd87
# ╠═0e5ed6c6-8293-464b-bc5b-c7b76ce7e421
# ╠═daa90cf6-7684-4b74-b6c6-6049cffede30
# ╠═762e14c1-3be1-4ba6-8236-cc417f1edcd1
# ╠═38d8f63d-c7cd-41cb-8ff3-21070a12b546
# ╠═9b90ebcc-f29f-424c-bc0a-821d1687a16e
# ╠═9f91ebb7-a492-44b2-a549-1f446c46ae34
# ╠═9aa3dc11-afc8-4b28-8146-0963ee0e5c5c
# ╠═2ecbaf82-d534-430b-a4b2-a7974ad7ce3b
# ╠═c3b20185-7ef0-4602-861d-5418bfb8a98f
# ╠═01628c9a-3c07-43bd-aa29-7db146633034
# ╠═62f7e767-3806-47f9-baa8-43c033528d07
