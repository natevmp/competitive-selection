using CairoMakie
using JLD2

l0 = 102
l1 = 251
l10 = 149/2
xTot_t = 10 * [0, 0, 11, 0, 13, 239, 336, 456, 426] ./ l10
n_t = [0, 0, 1, 0, 1, 18, 12, 15, 13]
_t = [0, 29, 38, 48, 63, 75, 76, 77, 81]
jldsave("./Data/variantNumbersSizesMitchell.jld2"; xTot_t=xTot_t./100, n_t, _t)
# load("./Data/variantNumbersSizesMitchell.jld2", "xTot_t", "n_t", "_t")

fig = scatterlines(_t, xTot_t, axis=(xlabel="age", ylabel="clonal fraction"))
xlims!(0,100)
# display(fig)

fig2 = scatterlines(_t, n_t, axis=(xlabel="age", ylabel="number of clones"))
xlims!(0,100)
# display(fig2)

colors = Makie.wong_colors()
fig3 = Figure()
Axis(fig3[1,1], xlabel="age",ylabel="clonal fraction", yticklabelcolor=colors[1], ylabelcolor=colors[1])
scatterlines!(_t, xTot_t, color=colors[1])
xlims!(0,100)
Axis(fig3[1,1], ylabel="number of clones", yaxisposition=:right, yticklabelcolor=colors[2], ylabelcolor=colors[2])
scatterlines!(_t, n_t, color=colors[2])
xlims!(0,100)
display(fig3)
