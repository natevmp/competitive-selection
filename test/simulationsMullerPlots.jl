include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using CSV, DataFrames
using Statistics
using CairoMakie

##

## ------------ Run sims -----------

paramSet = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 2000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
    :μ => 0,
    :s => 0.1,
    :σ => 0.01,
)
@time solEns, _ = CompetitiveSelection.evolvePopSim(
    paramSet;
    _trackerVariant=[
        (0.01,0.15),
        (5.,0.15),
        (25.,0.19),
        (50.,0.29),
        ],
    noDiffusion=true
)
@time solEnsSingleClone, _ = CompetitiveSelection.evolvePopSim(
    paramSet;
    _trackerVariant=[
        (0.01,0.15),
        ],
    noDiffusion=true
)



## --------------------------

sol = solEns[1]
size_vid_T = solEns[1].u
solSC = solEnsSingleClone[1]
sizeSC_t = [x_t[1] for x_t in solEnsSingleClone[1].u]

## ============== Muller plot functions =================

_t = range(0,100,length=length(sol))
n_t_vid0 = Array{Float64,2}(undef, length(sol), length(sol[1])+1)
nSC_t_vid0 = Array{Float64,2}(undef, length(solSC), length(solSC[1])+1)
for i in 1:length(sol[1])
    n_t_vid0[:,1+i] = [solT[i] for solT in sol.u]
end
n_t_vid0[:,1] = [1-sum(n_t_vid0[t,2:end]) for t in 1:size(n_t_vid0,1)]

##

include("../src/mullerPlot.jl")
using .MullerPlot

# using CSV
# jldsave("variantsizesOverTime.jld2", n_t_vid=n_t_vid)

nSC_t_vid0[:,2] = sizeSC_t
nSC_t_vid0[:,1] = 1 .- sizeSC_t

_child_Vid0 = [
    [2,3],
    [4],
    [5],
    Int[],
    Int[],
]

_child_Vid0SC = [[2,], Int[]]

xL_t_vid0, xU_t_vid0 = mullerBounds(n_t_vid0, _child_Vid0)
xL_t_vid0, xU_t_vid0 = mullerBounds(n_t_vid0, _child_Vid0)
xSCL_t_vid0, xSCU_t_vid0 = mullerBounds(nSC_t_vid0, _child_Vid0SC)

tIndBirth_vid = MullerPlot.variantBirthTimes(n_t_vid0)[2:end]
tBirth_vid = _t[tIndBirth_vid]
## ------------------ Figures ----------------------

using ColorSchemes
# colormap = ColorSchemes.tab10
colormap = ColorSchemes.tol_light

scale = 0.4
lw = 1.5
fig1 = Figure(resolution=(scale*1100,scale*800))
sizeWT_t = [1-sum(size_vid) for size_vid in size_vid_T]
axbottom = Axis(
    fig1[2,1], 
    palette=(color=cgrad(:Egypt),),
    xlabel="time",
    ylabel="size of clone",
    # xticksvisible=false,
    # xticklabelsvisible=false,
    ygridvisible=true,
    xgridvisible=true,
)
ylims!(0,1)
xlims!(0,100)
# hidedecorations!(axbottom, grid=false)
lines!(
    _t, sizeWT_t,
    color=:grey65,
    linewidth=lw,
)
for i in 1:length(size_vid_T[1])
    sizeV_t = [size_vid[i] for size_vid in size_vid_T]
    lines!(
        sizeV_t,
        # color=Cycled(i),
        color=colormap[i],
        linewidth=lw,
    )
end

# fig2 = Figure(resolution=(1100,500))
axtop = Axis(
    fig1[1,1],
    backgroundcolor=:grey65,
    # palette=(patchcolor=cgrad(:Egypt),),
    xgridvisible=false,
    ygridvisible=false,
    # xlabel="time",
    # ylabel="size of clone"
)
hidedecorations!(axtop)
for i in 1:(size(_child_Vid0,1)-1)
    band!(
        _t, xL_t_vid0[:,1+i], xU_t_vid0[:,1+i],
        show_grid=false,
        # color=Cycled(i)
        color=colormap[i]
    )
    scatter!(
        [tBirth_vid[i],],
        [xU_t_vid0[tIndBirth_vid[i], 1+i]+0.015,],
        color=colormap[i],
        # color=Cycled(i),
        # rotations=[-0.2,],
        # marker_offset=(0,0.1),
        # marker='↓',
        # markersize=25,
        # marker_offset=(0,10),
        # marker='▿',
        marker=:dtriangle,
        markersize=24,
    )
end
ylims!(0,1)
xlims!(0,100)
display(fig1)

save("Figures/Poster/clonalInterference_4Clones.svg", fig1, pt_per_unit = 1)
save("Figures/Conceptual/clonalInterference_4Clones.png", fig1, px_per_unit = 4)


## -------------------------- Single Clone Plot ---------------------

fig2 = Figure(resolution=(scale*1100,scale*1000))
sizeWT_t = 1 .- sizeSC_t
axbottom = Axis(
    fig2[2,1], 
    xlabel="time",
    ylabel="size of clone",
    # xticksvisible=false,
    # xticklabelsvisible=false,
    ygridvisible=true,
    xgridvisible=true,
)
ylims!(0,1)
xlims!(0,100)
# hidedecorations!(axbottom, grid=false)
lines!(
    _t, sizeWT_t,
    color=:grey65,
    linewidth=lw,
)
lines!(
    sizeSC_t,
    color=colormap[1],
    linewidth=lw,
)

# fig2 = Figure(resolution=(1100,500))
axtop = Axis(
    fig2[1,1],
    backgroundcolor=:grey65,
    xgridvisible=false,
    ygridvisible=false,
    # xlabel="time",
    # ylabel="size of clone"
)
hidedecorations!(axtop)
for i in 1:(size(_child_Vid0SC,1)-1)
    band!(
        _t, xSCL_t_vid0[:,1+i], xSCU_t_vid0[:,1+i],
        show_grid=false,
        color=colormap[i]
    )
end
ylims!(0,1)
xlims!(0,100)
display(fig2)

# save("Figures/Conceptual/clonalInterference_SingleClone.svg", fig2, pt_per_unit = 2)

## ========================================
#region = Combined Plots 


scale = 0.6
lw = 1.5
colormap = ColorSchemes.tol_light
wtColor=:grey75
mksize = 20
mkoffset = 0.05
fontsize=20

figFinal = Figure(resolution=(scale*1250,scale*800), fontsize=fontsize)

axbottomComp = Axis(
    figFinal[2,2], 
    ygridvisible=true,
    xgridvisible=true,
)
ylims!(0,1)
xlims!(0,100)
lines!(
    _t, sizeWT_t,
    color=wtColor,
    linewidth=lw,
)
for i in 1:length(size_vid_T[1])
    sizeV_t = [size_vid[i] for size_vid in size_vid_T]
    lines!(
        sizeV_t,
        color=colormap[i],
        linewidth=lw,
    )
end

axtopComp = Axis(
    figFinal[1,2],
    backgroundcolor=wtColor,
    xgridvisible=false,
    ygridvisible=false,
)
hidedecorations!(axtopComp)
for i in 1:(size(_child_Vid0,1)-1)
    band!(
        _t, xL_t_vid0[:,1+i], xU_t_vid0[:,1+i],
        show_grid=false,
        color=colormap[i]
    )
    scatter!(
        [tBirth_vid[i],],
        [xU_t_vid0[tIndBirth_vid[i], 1+i]+mkoffset,],
        color=colormap[i],
        marker=:dtriangle,
        markersize=mksize,
    )
end
ylims!(0,1)
xlims!(0,100)

sizeWTSC_t = 1 .- sizeSC_t
axbottomSingle = Axis(
    figFinal[2,1], 
    # ylabel="size of clone",
    ygridvisible=true,
    xgridvisible=true,
)
ylims!(0,1)
xlims!(0,100)
# hidedecorations!(axbottom, grid=false)
lines!(
    _t, sizeWTSC_t,
    color=wtColor,
    linewidth=lw,
)
lines!(
    sizeSC_t,
    color=colormap[1],
    linewidth=lw,
)

axtopSingle = Axis(
    figFinal[1,1],
    backgroundcolor=wtColor,
    xgridvisible=false,
    ygridvisible=false,
)
hidedecorations!(axtopSingle)
for i in 1:(size(_child_Vid0SC,1)-1)
    band!(
        _t, xSCL_t_vid0[:,1+i], xSCU_t_vid0[:,1+i],
        show_grid=false,
        color=colormap[i]
    )
    scatter!(
        [tBirth_vid[i],],
        [xU_t_vid0[tIndBirth_vid[i], 1+i]+mkoffset,],
        color=colormap[i],
        marker=:dtriangle,
        markersize=mksize,
    )
end
ylims!(0,1)
xlims!(0,100)

xAxis = Label(figFinal[3,:], "time", fontsize=1.2fontsize, padding=(0,0,0,-10))
yAxis = Label(figFinal[:,0], "clone size", fontsize=1.2fontsize,rotation=pi/2, padding=(0,-10,0,0))

display(figFinal)
save("Figures/Conceptual/clonalInterference_Combined.svg", figFinal, pt_per_unit = 2)

#endregion