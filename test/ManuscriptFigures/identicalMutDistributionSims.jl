include("../../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../../src/dataStructuring.jl")
using .DataStructuring

## ----------------------------------------
#region - test particle with fixed s and fixed t0

params = Dict(
    :N => Int(2E5),
    :α => 1, # wild type growth rate
    :s => 0.11, # mean of mutant fitness distribution
    :σ => 0.02,    # std of mutant fitness distribution
    :T => 90,   # evolve time
    :μ => 3.5,  # mutant arrival rates
    :sType => "gamma",  # shape of fitness distribution
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 10000
# runs = 1000
sLowTracker = 0.08
sHighTracker = 0.15
t0Tracker = 10

@time solEnsLow, simArgsLow = CompetitiveSelection.evolvePopSim(params; runs, _trackerVariant=[(t0Tracker, sLowTracker)]);
@time solEnsHigh, simArgsHigh = CompetitiveSelection.evolvePopSim(params; runs, _trackerVariant=[(t0Tracker, sHighTracker)]);
3
#endregion

## ----------------------------------------
#region - calculate distribution of fitnesses

include("../../src/dataStructuring.jl")
using .DataStructuring

function getFitnessDistTime(solEns, simArgs, params, tWindow_T, tPoints, vafThreshold)
    sFit_vid_T = Vector{Float64}[]
    for tWindow in tWindow_T
        _tMeasure = range(tWindow[1],tWindow[2],length=tPoints)
        sFit_vid = Float64[]
        for (sid, sol) in enumerate(solEns)
            vid = simArgs[sid, :_trackerID][1]
            x_tMeasure = [sol(t)[vid] for t in _tMeasure]
            sum(x_tMeasure.>=vafThreshold) < 3 && continue
            sFit = DataStructuring.fitModelsLightweight([x_tMeasure], _tMeasure, [1/vafThreshold for i in 1:tPoints], params[:N], ["positive logistic", "negative logistic"]; errorQuantile=0.99)[1]
            ismissing(sFit) && continue
            push!(sFit_vid, sFit)
        end
        push!(sFit_vid_T, sFit_vid)
    end
    return sFit_vid_T
end

##
tPoints = 4
tWindow_T = [(50,60),(60,70),(70,80),(80,90)]
vafThreshold = 0.001

@time sLowFit_vid_T = getFitnessDistTime(solEnsLow, simArgsLow, params, tWindow_T, tPoints, vafThreshold)
@time sHighFit_vid_T = getFitnessDistTime(solEnsHigh, simArgsHigh, params, tWindow_T, tPoints, vafThreshold)
#endregion

## ----------------------------------------
#region - save/load fitness dists

using JLD2, FileIO

# jldsave("./Data/observedFitnessIdenticalMuts.jld2",
#     sLowFit_vid_T=sLowFit_vid_T,
#     sHighFit_vid_T=sHighFit_vid_T,
#     tWindow_T=tWindow_T,
#     vafThreshold=vafThreshold,
#     params=params,
#     sLowTracker=sLowTracker,
#     sHighTracker=sHighTracker,
# )

(
    sLowFit_vid_T,
    sHighFit_vid_T,
    vafThreshold,
    tWindow_T,
    sLowTracker,
    sHighTracker,
    params,
) = load("./Data/observedFitnessIdenticalMuts.jld2",
    "sLowFit_vid_T",
    "sHighFit_vid_T",
    "vafThreshold",
    "tWindow_T",
    "sLowTracker",
    "sHighTracker",
    "params",
)
#endregion

## ========================================
#region = Plotting

using CairoMakie
using Statistics

sLowAv_t = mean.(sLowFit_vid_T)
sHighAv_t = mean.(sHighFit_vid_T)

#endregion

## ----------------------------------------
#region - high fitness variant

offset = -30
_bins = [30, 30, 45, 65,]

_yticks = [
    "50-60 years",
    "60-70 years",
    "70-80 years",
    "80-90 years",
]
fig = Figure(fontsize=26, resolution=(1000,600))
Axis(
    fig[1,1],
    xlabel="observed logistic fitness",
    # ylabel="density of variants",
    # yticklabelsvisible=false,
    # yticksvisible=false,
    yticks=([0+10, offset+10, 2offset+10, 3offset+10],_yticks),
    ygridvisible=false,
    # xgridvisible=false,
    title="high innate fitness",
)
for i in eachindex(sHighFit_vid_T)
    hist!(sHighFit_vid_T[i],
        label=string(tWindow_T[i][1])*"-"*string(string(tWindow_T[i][2]))*" years",
        normalization=:pdf,
        offset=offset*(i-1),
        bins=_bins[i],
        color=(:black,0.35),
        # color=(Makie.wong_colors()[1],0.5)
    )
    scatter!(sHighAv_t[i], offset*(i-1)+10,
        markersize=20,
        marker=:diamond,
        color=:black,
        # color=(Makie.wong_colors()[1],1),
    )
end
vlines!(sHighTracker;
    # label="mean",
    color=:black,
    linewidth=3,
    # color=(Makie.wong_colors()[1],1),
)
xlims!(-0.1,0.25)
# axislegend(position=:lt)
# display(fig)
#endregion

## ----------------------------------------
#region - separate low/high fitness

offset = -60
# fig = Figure(fontsize=26)
Axis(
    fig[1,2],
    xlabel="observed logistic fitness",
    # ylabel="density of variants",
    yticklabelsvisible=false,
    yticksvisible=false,
    ygridvisible=false,
    # xgridvisible=false,
    title="low innate fitness",
)
for i in eachindex(sLowFit_vid_T)
    hist!(sLowFit_vid_T[i],
        label=string(tWindow_T[i][1])*"-"*string(string(tWindow_T[i][2]))*" years",
        normalization=:pdf,
        offset=offset*(i-1),
        bins=30,
        color=(:black,0.35),
    )
    scatter!(sLowAv_t[i], offset*(i-1)+10,
        markersize=20,
        color=:black,
        marker=:diamond,
    )
end
vlines!(sLowTracker;
    # label="mean",
    color=:black,
    linewidth=3,
)
xlims!(-0.15,0.2)
ylims!(3*offset-12, 74)
display(fig)

figname="identicalMutFitnessDistributionV1.svg"
# save("./Figures/ManuscriptDrafts/"*figname, fig)
#endregion

## ----------------------------------------
#region - separate low/high fitness Vertical

offset = 1
barscale = 0.5
fig = Figure(fontsize=28, resolution=(800,600))
_bins = [25, 30, 44, 50,]
_yearticks = [
    "50-60 y.o.",
    "60-70 y.o.",
    "70-80 y.o.",
    "80-90 y.o.",
]
ga = fig[1:2,1] = GridLayout()
Axis(
    ga[1,1],
    xgridvisible=false,
    xticklabelsvisible=false,
    xticksvisible=false,
)
for i in eachindex(sHighFit_vid_T)
    hist!(sHighFit_vid_T[i],
        normalization=:pdf,
        offset=offset*(i-1),
        bins=_bins[i],
        color=:grey65,
        direction=:x,
        scale_to=-barscale,
    )
    scatter!(offset*(i-1), sHighAv_t[i],
        markersize=20,
        color=:grey25,
        marker=:diamond,
        # label="s = 0.15"
    )
end
hlines!(sHighTracker;
    color=:grey25,
    linewidth=3,
)
xlims!(-0.7,3.6)
ylims!(-0.08,0.3)

_bins = [20, 20, 20, 20,]
Axis(
    ga[2,1],
    xticks=([0, offset, 2offset, 3offset],_yearticks),
    xgridvisible=false,
)
for i in eachindex(sLowFit_vid_T)
    hist!(sLowFit_vid_T[i],
        normalization=:pdf,
        offset=offset*(i-1),
        bins=_bins[i],
        direction=:x,
        scale_to=-barscale,
        color=:grey65,
    )
    scatter!( offset*(i-1)-0.2, sLowAv_t[i],
        markersize=20,
        color=:grey25,
        marker=:diamond,
    )
end
hlines!(sLowTracker;
    # label="mean",
    # color=:black,
    color=:grey25,
    linewidth=3,
    # label="s = 0.08",
)
xlims!(-0.7,3.6)
ylims!(-0.15,0.18)

rowgap!(ga, 10)

yLabel = Label(fig[1:2,0], "observed logistic fitness", rotation=pi/2, padding=(0,-20,0,0))
Legend(fig[0,1], [legendInnate, legendObs, LegendMean], ["innate fitness", "observed fitness", "mean fitness"], orientation=:horizontal, padding=(0,0,-15,0), framevisible=false)

display(fig)
figname="identicalMutFitnessDistributionV4.svg"
save("./Figures/ManuscriptDrafts/"*figname, fig)

#endregion

## ----------------------------------------
#region - Simultaneous low & high

offset = 1
barscale = 0.5
fig = Figure(fontsize=28, resolution=(800,600))
_bins = [25, 30, 44, 66,]
_yearticks = [
    "50-60 y.o.",
    "60-70 y.o.",
    "70-80 y.o.",
    "80-90 y.o.",
]
Axis(
    fig[1,1],
    ylabel="observed logistic fitness",
    # xticks=([0+10, offset+10, 2offset+10, 3offset+10],_yearticks),
    xticks=([0, offset, 2offset, 3offset],_yearticks),
    xgridvisible=false,
)
for i in eachindex(sHighFit_vid_T)
    hist!(sHighFit_vid_T[i],
        normalization=:pdf,
        offset=offset*(i-1),
        bins=_bins[i],
        color=(Makie.wong_colors()[1],0.5),
        direction=:x,
        scale_to=barscale,
    )
    # scatter!(sHighAv_t[i], offset*(i-1)+10,
    scatter!(offset*(i-1)+0.2, sHighAv_t[i],
        markersize=20,
        color=Makie.wong_colors()[1],
        marker=:diamond,
        # label="s = 0.15"
    )
end
hlines!(sHighTracker;
    # label="mean",
    # color=:black,
    color=Makie.wong_colors()[1],
    linewidth=3,
    # label="s = 0.15",
)

for i in eachindex(sLowFit_vid_T)
    hist!(sLowFit_vid_T[i],
        normalization=:pdf,
        offset=offset*(i-1),
        bins=30,
        direction=:x,
        scale_to=-barscale,
        color=(Makie.wong_colors()[2],0.5),
    )
    # scatter!(sLowAv_t[i], offset*(i-1)+10,
    scatter!( offset*(i-1)-0.2, sLowAv_t[i],
        markersize=20,
        color=(Makie.wong_colors()[2],1),
        marker=:diamond,
    )
end
hlines!(sLowTracker;
    # label="mean",
    # color=:black,
    color=(Makie.wong_colors()[2],1),
    linewidth=3,
    # label="s = 0.08",
)
# legend plots
legendInnate = hlines!(
    0.3,
    color=:grey,
    linewidth=3,
    label="innate",
)
legendObs = barplot!(
    50,50,
    color=(:grey,0.5),
    label="observed",
)
LegendMean = scatter!(
    50,50,
    marker=:diamond,
    color=:grey,
    markersize=16,
    label="mean",
)
xlims!(-0.7,3.6)
ylims!(-0.1,0.2)
axislegend(position=:lb)
# Legend(fig[1,2], [legendInnate, legendObs, LegendMean], ["innate", "observed", "mean"])
display(fig)

figname="identicalMutFitnessDistributionV3.svg"
# save("./Figures/ManuscriptDrafts/"*figname, fig)
#endregion

