include("../src/mullerPlot.jl")
using .MullerPlot
using JLD2, FileIO

n_t_vid = load("variantsizesOverTime.jld2", "n_t_vid")
# this is a 101x5 matrix. The 101 rows are indexed by the measurement timepoints and the 5 columns indexed by the different variants. Each variant has an id corresponding to its row-number, whereby the first column represents the wild type. So for example, `n_t_vid[:,3]` is a vector containing all time measurements of the variant with id=3. It's probably best to sort them by order of appearance, though this is not strictly necessary.

_child_Vid = [
    [2,3],
    [4],
    [5],
    [],
    [],
]
# this vector contains the direct children of each variant. So element `i` is a list of variant id's corresponding the the children of variant `i`. So in this case, the wildtype (vid=1) has two children (variants that occurred within `i`), variant2 and variant3. Going down the list, variant2 has only one child (variant4), and variant3 also has one child (variant5)

xL_t_vid, xU_t_vid = mullerBounds(n_t_vid, _child_Vid)
# this function creates the bounds for each variant, to be used with a band plot.


## ======== Figure ========
using CairoMakie

_t = 0:size(n_t_vid,1)-1
fig = Figure(resolution=(550,250))
axis = Axis(
    fig[1,1],
    backgroundcolor=:grey65,
    palette=(patchcolor=cgrad(:Egypt),),
    xgridvisible=false,
    ygridvisible=false,
    xlabel="time",
)
hidedecorations!(axis)
for i in 2:length(_child_Vid)
    band!(
        _t, xL_t_vid[:,i], xU_t_vid[:,i],
        show_grid=false,
        color=Cycled(i)
    )
end
ylims!(0,1)
xlims!(0,100)
display(fig)
