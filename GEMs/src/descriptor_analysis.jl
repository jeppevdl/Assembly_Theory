using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using DataFrames, CSV, CairoMakie, Plots, Statistics

# read in descriptor data
df = CSV.read("../data/pathway_complexities/complete_descriptors.csv", DataFrame)
palet = palette(:seaborn_colorblind);

# plot exact molecular mass vs MA
mass_ma = Figure();
ax1 = Axis(mass_ma[1, 1], title = "Mass vs MA", xlabel = "Exact mass (g/mol)", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.scatter!(ax1, df.mass, df.ma, markersize = 3, color = palet[1])
CairoMakie.text!(ax1, 0.5, 140, text = "Pearson correlation: $(round(cor(df.mass, df.ma), digits = 3))", color = :red, fontsize = 20)
mass_ma
CairoMakie.save("../figures/descriptor_analysis/mass_ma.png", mass_ma, px_per_unit=5)

# plot number of atoms vs MA
atoms_ma = Figure();
ax2 = Axis(atoms_ma[1, 1], title = "Number of atoms vs MA", xlabel = "Number of atoms", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.scatter!(ax2, df.counted_atoms, df.ma, markersize = 3, color = palet[2])
CairoMakie.text!(ax2, 0.5, 140, text = "Pearson correlation: $(round(cor(df.counted_atoms, df.ma), digits = 3))", color = :red, fontsize = 20)
atoms_ma
CairoMakie.save("../figures/descriptor_analysis/atoms_ma.png", atoms_ma, px_per_unit=5)

# plot number of rings vs MA
rings_ma = Figure();
ax3 = Axis(rings_ma[1, 1], title = "Number of rings vs MA", xlabel = "Number of rings", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.boxplot!(ax3, df.counted_rings, df.ma, markersize = 3, color = palet[3])
for i in unique(df.counted_rings)
    n = sum(df.counted_rings .== i)
    CairoMakie.text!(ax3, i, maximum(df.ma) + 1, text = "n=$n", align = (:center, :bottom), color = :black, fontsize = 15.5)
end
rings_ma
CairoMakie.save("../figures/descriptor_analysis/rings_ma.png", rings_ma, px_per_unit=5)

# plot number of sp-hybridized atoms vs MA
sp_ma = Figure();
ax4 = Axis(sp_ma[1, 1], title = "sp-hybridization vs MA", xlabel = "Number of sp-hybridized atoms", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.boxplot!(ax4, df.sp, df.ma, markersize = 3, color = palet[4])
for i in unique(df.sp)
    n = sum(df.sp .== i)
    CairoMakie.text!(ax4, i, maximum(df.ma) + 1, text = "n=$n", align = (:center, :bottom), color = :black, fontsize = 15.5)
end
sp_ma
CairoMakie.save("../figures/descriptor_analysis/sp_ma.png", sp_ma, px_per_unit=5)

# plot number of sp2-hybridized atoms vs MA
sp2_ma = Figure();
ax5 = Axis(sp2_ma[1, 1], title = "sp2-hybridization vs MA", xlabel = "Number of sp2-hybridized atoms", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.scatter!(ax5, df.sp2, df.ma, markersize = 3, color = palet[5])
CairoMakie.text!(ax5, 0.5, 140, text = "Pearson correlation: $(round(cor(df.sp2, df.ma), digits = 3))", color = :red, fontsize = 20)
sp2_ma
CairoMakie.save("../figures/descriptor_analysis/sp2_ma.png", sp2_ma, px_per_unit=5)

# plot number of sp3-hybridized atoms vs MA
sp3_ma = Figure();
ax6 = Axis(sp3_ma[1, 1], title = "sp3-hybridization vs MA", xlabel = "Number of sp3-hybridized atoms", ylabel = "MA",
titlesize = 25, xlabelsize = 22, ylabelsize = 22, xticklabelsize = 20, yticklabelsize = 20)
CairoMakie.scatter!(ax6, df.sp3, df.ma, markersize = 3, color = palet[6])
CairoMakie.text!(ax6, 0.5, 140, text = "Pearson correlation: $(round(cor(df.sp3, df.ma), digits = 3))", color = :red, fontsize = 20)
sp3_ma
CairoMakie.save("../figures/descriptor_analysis/sp3_ma.png", sp3_ma, px_per_unit=5)