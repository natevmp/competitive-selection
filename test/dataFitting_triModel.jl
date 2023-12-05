using CSV, DataFrames, Statistics, StatsBase
include("../src/dataStructuring.jl")
using .DataStructuring
using ProgressMeter
using CairoMakie
include("../src/analysisTools.jl")
using .AnalysisTools

using JLD2
SAVEFIG = false

## ========================================
#region = Data Analysis
using FileIO

df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
Nf = 2E5
DataStructuring.structureData!(df)
# dfVid = DataStructuring.sortDataPerVariant(df)
# DataStructuring.fitAllModels!(dfVid, Nf, ["positive logistic", "negative logistic", "constant"]; errorQuantile=0.99)
# DataStructuring.fitAllModels!(dfVid, Nf, ["positive logistic", "negative logistic"]; errorQuantile=0.99)
# DataStructuring.testFit!(dfVid, 0.99)
# jldsave("Data/dfVidFittedTriModel.jld2"; df=df, dfVid=dfVid)
dfVid = load("Data/dfVidFittedTriModel.jld2", "dfVid")

#endregion

## ========================================
#region = Data masks

_posMask = dfVid[!,:fitType] .== "positive logistic"
_negMask = dfVid[!,:fitType] .== "negative logistic"
_satMask = dfVid[!,:fitType] .== "constant"
_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]
# _λMask = [dfVCur[:γ]<grThreshold for dfVCur in eachrow(dfVid)]
_fitMask = dfVid[!,:goodFit]
sum(_fitMask)

size(dfVid[_satMask.&&_fitMask,:],1) / size(dfVid,1)
size(dfVid[_posMask.&&_fitMask,:],1) / size(dfVid,1)
size(dfVid[_negMask.&&_fitMask,:],1) / size(dfVid,1)
size(dfVid[.!_fitMask,:],1) / size(dfVid,1)

#endregion

## ========================================
#region = Plot colors

colorsTheme = Makie.wong_colors()
colorsTheme = cgrad(:roma10, rev=true)
colorCycle(i) = Makie.wong_colors()[i]
colorType(typeS) = 
    if typeS=="positive logistic"
        colorCycle(1)
    elseif typeS=="negative logistic"
        colorCycle(2)
    elseif typeS=="constant"
        colorCycle(3)
    else
        nothing
    end
#endregion

## ========================================
#region = Size distributions

fbins = 15
_f = zeros(fbins)
xMin = 0.0002
xMax = 0.25
# dfVid = load("./Data/dfDataFittedVid.jld2", "dfVid")
x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
xMask_sid = x_sid.>0
_tEdges = [50,65,80,95]
n_f_T = Vector{Vector{Float64}}(undef, length(_tEdges)-1)
for tid in 1:(length(_tEdges)-1)
    _fEdges, nData_f = AnalysisTools.sizeDistDens(
        x_sid[xMask_sid], t_sid[xMask_sid], (_tEdges[tid],_tEdges[tid+1]); bins=fbins, xMin=xMin, xMax=xMax
    )
    _f .= _fEdges[2:end] .- step(_fEdges)/2
    n_f_T[tid] = nData_f
end

figSizeDist = Figure(fontsize=22)
Axis(figSizeDist[1,1],
    xscale=log10,
    yscale=log10,
    xlabel="variant size",
    ylabel="(normalized) number of variants"
)

ylims!(1E-3, 1)
for tid in eachindex(n_f_T)
    lines!(_f, n_f_T[tid] ./ (sum(nData_f)),
        label="ages "*string(_tEdges[tid])*"-"*string(_tEdges[tid+1])
    )
end
axislegend()
display(figSizeDist)

#endregion

## --------------------------------------------------
#region - Stacked bars trajectory types

tbl = (x = [1, 1, 2,2, 3, 3],
    height = [
        size(dfVid[_posMask.&&_fitMask,:],1) / size(dfVid,1),
        size(dfVid[_posMask.&&(.!_fitMask),:],1) / size(dfVid,1),
        size(dfVid[_negMask.&&_fitMask,:],1) / size(dfVid,1),
        size(dfVid[_negMask.&&(.!_fitMask),:],1) / size(dfVid,1),
        size(dfVid[_satMask.&&_fitMask,:],1) / size(dfVid,1),
        size(dfVid[_satMask.&&(.!_fitMask),:],1) / size(dfVid,1),
        ],
    grp = [1, 2, 1, 2, 1, 2],
    labels = ["successful fit", "unsuccessful fit"],
    elements = [PolyElement(polycolor = (:grey35,1)), PolyElement(polycolor = (:grey35,0.4))],
)

propFit = size(dfVid[_fitMask,:],1) / size(dfVid,1)
propPos = size(dfVid[_posMask.&&_fitMask,:],1) / size(dfVid[_fitMask,:],1)
propNeg = size(dfVid[_negMask.&&_fitMask,:],1) / size(dfVid[_fitMask,:],1)
propSat = size(dfVid[_satMask.&&_fitMask,:],1) / size(dfVid[_fitMask,:],1)
println("positive trajectories: ", propPos)
println("negative trajectories: ", propNeg)
println("saturating trajectories: ", propSat)


rescale = 0.9
figTypes = Figure(
    fontsize=26,
    resolution=(rescale*750,rescale*600)
)
ax = Axis(
    xgridvisible=false,
    ygridvisible=false,
    figTypes[1,1], 
    xticks=(1:3, ["increasing", "decreasing", "saturating"]),
    ylabel="proportion of trajetories"
)
barplot!(tbl.x, tbl.height,
    stack = tbl.grp,
    color = [colorCycle(1), (colorCycle(1),0.5), colorCycle(2), (colorCycle(2),0.5), colorCycle(3), (colorCycle(3),0.5)],
)
axislegend(ax, tbl.elements, tbl.labels)
display(figTypes)

# figname="trajectoryFitTypes.png"
figname="trajectoryFitTypes.svg"
figloc="./Figures/Poster/"
SAVEFIG && save(figloc*figname, figTypes)
#endregion

## --------------------------------------------------
#region - trajectory type per measurement age: Barplot

measureDuration_vid = [_t[end]-_t[1] for _t in dfVid[!, :_t]]
measureAge_vid = [mean(_t) for _t in dfVid[!,:_t]]


_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfVid[!,:vid],
    bins=13
)

_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nPos_t = [sum((vid->in(vid,vid_t)).(dfVid[_posMask,:vid])) for vid_t in vid_id_T]
nSat_t = [sum((vid->in(vid,vid_t)).(dfVid[_satMask,:vid])) for vid_t in vid_id_T]
nNeg_t = [sum((vid->in(vid,vid_t)).(dfVid[_negMask,:vid])) for vid_t in vid_id_T]
nVars_t = length.(vid_id_T)

weight_t = nVars_t ./ sum(nVars_t)

figTypeAge = Figure(fontsize=26)
ax = Axis(figTypeAge[1,1],
    xlabel="mean age at measurement",
    ylabel="proportion of trajectories"
)
# hist!(measureAge_vid[_posMask], normalization=:probability)
barplot!(
    vcat(_t,_t,_t),
    vcat(nPos_t./nVars_t, nNeg_t./nVars_t, nSat_t./nVars_t),
    stack=vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))),
    color=colorCycle.( vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))) )
)
ylims!(0,1)
# Legend(figTypeAge[1,2],
#     [PolyElement(polycolor = colorCycle(1)), PolyElement(polycolor = colorCycle(2)), PolyElement(polycolor = colorCycle(3))],
#     ["positive", "negative", "saturating"],
# )
axislegend(
    ax, [PolyElement(polycolor = colorCycle(1)), PolyElement(polycolor = colorCycle(2)), PolyElement(polycolor = colorCycle(3))],
    ["positive", "negative", "saturating"],
    position=:lb
)
display(figTypeAge)

figname="trajectoryTypesTimeBar.svg"
figloc="./Figures/ManuscriptDrafts/"
SAVEFIG && save(figloc*figname, figTypeAge)
#endregion

## --------------------------------------------------
#region - trajectory type per measurement age: Barplot 2

figTypeAge = Figure(fontsize=22)
ax = Axis(figTypeAge[1,1],
    # xlabel="mean age at measurement",
    xticks=(1:3, ["increasing", "decreasing", "saturating"]),
    ylabel="mean age at measurement"
)
# hist!(measureAge_vid[_posMask], normalization=:probability)
barplot!(
    _t,
    (nPos_t./nVars_t).*-0.8,
    # stack=vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))),
    # color=colorsTh.( vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))) ),
    direction=:x,
    # scale_to=-0.6,
    offset=1,
)
barplot!(
    _t,
    (nNeg_t./nVars_t).*-0.8,
    # stack=vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))),
    # color=colorsTh.( vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))) ),
    direction=:x,
    # scale_to=-0.6,
    offset=2,
)
barplot!(
    _t,
    (nSat_t./nVars_t).*-3,
    # stack=vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))),
    # color=colorsTh.( vcat(fill(1,length(_t)), fill(2,length(_t)), fill(3,length(_t))) ),
    direction=:x,
    # scale_to=-0.6,
    offset=3,
)
display(figTypeAge)
#endregion

##--------------------------------------------------
#region - trajectory type abundance with age

measureAge_vid = [mean(_t) for _t in dfVid[!,:_t]]
_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfVid[!,:vid],
    bins=20
)
_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nPos_t = [sum((vid->in(vid,vidT_id)).(dfVid[_posMask,:vid])) for vidT_id in vid_id_T]
nSat_t = [sum((vid->in(vid,vidT_id)).(dfVid[_satMask,:vid])) for vidT_id in vid_id_T]
nNeg_t = [sum((vid->in(vid,vidT_id)).(dfVid[_negMask,:vid])) for vidT_id in vid_id_T]
nVars_t = [sum((vid->in(vid,vid_t)).(dfVid[!,:vid])) for vid_t in vid_id_T]

nNeg_t ./ nVars_t
_tEdges[1]
_tEdges[2]
_t[1]
_t[end-2]

ρIncreasing = corspearman(_t[.!isnan.(nPos_t./nVars_t)], (nPos_t./nVars_t)[.!isnan.(nPos_t./nVars_t)])
ρDecreasing = corspearman(_t[.!isnan.(nNeg_t./nVars_t)], (nNeg_t./nVars_t)[.!isnan.(nNeg_t./nVars_t)])

figTrajTime = Figure(
    fontsize=26,
    resolution=(rescale*750,rescale*600)
)
Axis(figTrajTime[1,1],
    ylabel="fraction of trajectories increasing",
    xlabel="mean age at measurement",
)
scatterlines!(
    # _t, nPos_t./nVars_t,
    _t, nNeg_t ./ nVars_t,
    color=:black,
    # label="spearman correlation: "*string(round(ρIncreasing,digits=2))
    label="spearman correlation: "*string(round(ρDecreasing,digits=2))
)

axislegend(position=:rb)
# Axis(figTrajTime[1,2])
# scatter!(
#     _t, (nNeg_t.+nSat_t)./nVars_t,
#     color=colorCycle(2),
#     label="decreasing or saturating\ntrajectories"
# )
# axislegend(position=:rb)
# Label(figTrajTime[2,:], "mean age at measurement")
display(figTrajTime)

# figname="trajectoryTypesTime.png"
figname="trajectoryTypesTime.svg"
figloc="./Figures/Poster/"
SAVEFIG && save(figloc*figname, figTrajTime)

#endregion

## --------------------------------------------------
#region - Single trajectories

# rscale=0.7
# res = (rscale*600,rscale*500)
# fsize=16

# for (i,dfVCur) in enumerate(eachrow(dfVid[20:80,:]))
#     println("t0: ", dfVCur[:t0])
#     println("γ: ", dfVCur[:γ])
#     # println("x0: ", dfVCur[:t0])
#     fig1 = Figure(resolution=res, fontsize=fsize)
#     ax = Axis(
#         fig1[1,1],
#         title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene]*"; vid="*string(dfVCur[:vid]),
#     )
#     scatter!(dfVCur[:_t], dfVCur[:vaf_t],
#         linestyle=:solid,
#         label=dfVCur[:gene],
#     )
#     lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
#         linestyle=:dash,
#         label="s = "*string(round(dfVCur[:γ],digits=2)),
#     )
#     rangebars!(
#         dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
#         whiskerwidth = 20,
#         cycle = [:color],
#     )
#     legendpos = ( dfVCur[:ρ]>0 ? :lt : :rt )
#     axislegend(position=legendpos)
#     # figname = "dataFit_vid$i.png"
#     # save("Figures/DataFits/PerVariant/"*figname, fig1)
#     display(fig1)
# end
#endregion

## --------------------------------------------------
#region - multiple trajectories in single figure

rscale=0.8
res = (1300,1100)
fsize=22
nRows = 3
nCols = 3
_plotvid = [
    3,4,7,
    439,13,542,
    79,78,76,
]
figTraj = Figure(resolution=res, fontsize=fsize)
for (i,dfVCur) in enumerate(eachrow(dfVid[_plotvid,:]))
    col = (i-1) % nCols + 1
    row = floor((i-1)/nCols) + 1 |> Int
    color = colorType(dfVCur[:fitType])
    color_t = [(dfVCur[:fit_t][t] ? :black : :red) for t in eachindex(dfVCur[:fit_t])]
    marker_t = [(dfVCur[:fit_t][t] ? :diamond : :xcross) for t in eachindex(dfVCur[:fit_t])]
    ax = Axis(
        figTraj[row,col],
        title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene],
        xlabel=(row==nRows ? "age at measurement" : ""),
        ylabel=(col==1 ? "VAF" : ""),
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label=dfVCur[:fitType]*"\n"
             *"s = "*string(round(dfVCur[:γ],digits=2)),
        # color=color,
        color=:grey35,
        # color=colorsTheme[1],
        linewidth=2,
    )
    rangebars!(
        dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
        whiskerwidth = 20,
        # color=color,
        color=:grey35,
        # color=colorsTheme[1],
        linewidth=1.5,
    )
    for t in eachindex(dfVCur[:_t])
        scatter!(dfVCur[:_t][t], dfVCur[:vaf_t][t],
            # linestyle=:solid,
            # label=dfVCur[:gene],
            marker=marker_t[t],
            color=color_t[t],
            markersize=17,
        )
    end
    legendpos = ( dfVCur[:ρ]>0 ? :lt : :rt )
    axislegend(position=legendpos, bgcolor=(:white,0.), framevisible=false)
end
display(figTraj)

figname="fittedTrajectories.svg"
figloc="./Figures/Poster/"
# save(figloc*figname, figTraj)
#endregion

## --------------------------------------------------
#region - multiple trajectories in single plot
names(dfVid)
dfVid[(dfVid[!,:γ].<0) .&& (x->x[1]>0.1).(dfVid[!,:vaf_t]),[:vid,:vaf_t]]
dfVid[(x->x[1]<60).(dfVid[!,:_t]) .&& (x->x[1]>0.02).(dfVid[!,:vaf_t]),[:vid,:vaf_t]]

# _vidPlot = [6,10,16,34]
# _vidPlot = [10,16,34]
_vidPlot = [697, 107,105, 461]
rsc=0.5
res = (rsc*800,rsc*600)
fsize=20
fig0 = Figure(resolution=res, fontsize=fsize)
ax = Axis(
    fig0[1,1],
    # yticks=[0,0.02,0.04,0.06],
    xlabel="age at measurement",
    ylabel="variant allele frequency",
    xgridvisible=false,
    ygridvisible=false,
    rightspinevisible=false,
    topspinevisible=false,
)
# hidespines!(ax, :t, :r)
_scatMarker=[:circle, :diamond, :utriangle, :pentagon]
for (i,dfVCur) in enumerate(eachrow(dfVid[_vidPlot,:]))
    # colorCur=Makie.wong_colors()[i]
    colorCur=:black
    scatter!(
        dfVCur[:_t],
        # dfVCur[:_t] .- dfVCur[:_t][1],
        # dfVCur[:vaf_t] .- dfVCur[:vaf_t][1],
        dfVCur[:vaf_t],
        linestyle=:solid,
        label=dfVCur[:gene],
        color=colorCur,
        marker=_scatMarker[i],
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label=dfVCur[:fitType]*"\n"
            *"s = "*string(round(dfVCur[:γ],digits=2)),
        linewidth=2,
        color=(colorCur,0.4),
    )
    rangebars!(
        dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
        whiskerwidth = 20,
        color=(colorCur,0.4),
        linewidth=1.5,
    )
    fsize=15
    text!(
        dfVCur[:_t][end]+0.5, dfVCur[:vaf_t][end]+0.0;
        text="pid"*string(dfVCur[:pid]),
        fontsize=fsize,
        color=colorCur,
    )
    text!(
        dfVCur[:_t][end]+0.5, dfVCur[:vaf_t][end]-0.03;
        text=string(dfVCur[:gene]),
        fontsize=fsize,
        color=colorCur,
    )
end
xlims!(54,95)
ylims!(-0.01,0.43)
display(fig0)
figname="singleTrajectories.svg"
figloc="./Figures/Poster/"
# save(figloc*figname, fig0, px_per_unit=4)
#endregion

## --------------------------------------------------
#region - Fitted growth rates distribution

rescale = 0.9
figFitnessDist = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=24,
)
Axis(
    figFitnessDist[1,1],
    xlabel="logistic fitness",
    ylabel="density of variants"
)
hist!(
    dfVid[_fitMask,:γ],
    bins=54,
    normalization=:pdf,
    strokewidth=2,
    color=:grey45,
    strokecolor=:grey45,
)
display(figFitnessDist)
figname="variantFitnessDistribution.svg"
figloc="./Figures/Poster/"
# save(figloc*figname, figFitnessDist)
#endregion

## --------------------------------------------------
#region - Fitted Growth rates with time
rscale = 0.8
figGRTime = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=24,
)
grThreshold = 11
# DataStructuring.testFit!(dfVid, 0.99)
markersize = 7
Axis(
    figGRTime[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_posMask,:])],
    dfVid[_fitMask.&&_posMask, :γ];
    markersize,
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_negMask,:])],
    dfVid[_fitMask.&&_negMask, :γ];
    markersize,
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_satMask,:])],
    dfVid[_fitMask.&&_satMask, :γ];
    markersize,
)
# ylims!(-1,1)
display(figGRTime)

figname="dataGrowthRatesTime.png"
figloc="./Figures/ManuscriptDrafts/Figure1/"
SAVEFIG && save(figloc*figname, figGRTime, px_per_unit=4)
#endregion

## ----------------------------------------
#region - Sizes of clones for all categories

figSizeSat = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=24,
)
Axis(figSizeSat[1,1],
    xlabel="age of donor",
    ylabel="size of clone",
    yscale=log10
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_posMask,:])],
    [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_posMask,:])];
    markersize,
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_negMask,:])],
    [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_negMask,:])];
    markersize,
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask.&&_satMask,:])],
    [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_satMask,:])];
    markersize,
)
display(figSizeSat)

#endregion

## ----------------------------------------
#region - Sizes of clones boxplots

# xs = rand(1:3, 1000)
# ys = randn(1000)
# boxplot(xs, ys)



_tBinRange = [20,70,80,100]
_tBin = _tBinRange[1:end-1] + (_tBinRange[2:end].-_tBinRange[1:end-1])./2
vaf_vidPos_T = Vector{Vector{Float64}}(undef, length(_tBinRange)-1)
vaf_vidNeg_T = similar(vaf_vidPos_T)
vaf_vidSat_T = similar(vaf_vidPos_T)
for tid in 1:length(_tBinRange)-1
    _tMask = (_t->_tBinRange[tid] .< mean(_t) .<= _tBinRange[tid+1]).(dfVid[!,:_t])
    vaf_vidPos_T[tid] = [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_posMask.&&_tMask,:])]
    vaf_vidNeg_T[tid] = [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_negMask.&&_tMask,:])]
    vaf_vidSat_T[tid] = [mean(dfVidCur[:vaf_t]) for dfVidCur in  eachrow(dfVid[_fitMask.&&_satMask.&&_tMask,:])]
end

figSizeBox = Figure(fontsize=24)
Axis(figSizeBox[1,1],
    xticks=([4,8,12], 
            [
                "<"*string(_tBinRange[2]),
                string(_tBinRange[2])*"-"*string(_tBinRange[3]),
                ">"*string(_tBinRange[3])
            ]),
    yticks=[10^(-3),10^(-2),10^(-1),10^0],
    ylabel="Clone size",
    xlabel="Age (years)",
    yscale=log10,
    xgridvisible=false,
)
ylims!(0.003,1)
colorCycle(1)
for (tid,t) in enumerate(_tBin)
    xLoc = 4*tid
    boxplot!(
        (xLoc-1)*ones(length(vaf_vidPos_T[tid])),
        2*vaf_vidPos_T[tid],
        show_outliers=false,
        color=(colorCycle(1),0.5),
    )
    scatter!(
        (xLoc-1) .+ 0.1 .*randn(length(vaf_vidPos_T[tid])).*ones(length(vaf_vidPos_T[tid])), 2*vaf_vidPos_T[tid],
        color=colorCycle(1),
    )
    boxplot!(
        xLoc*ones(length(vaf_vidNeg_T[tid])),
        2*vaf_vidNeg_T[tid],
        show_outliers=false,
        color=(colorCycle(2),0.5),
    )
    scatter!(
        xLoc .+ 0.1 .*randn(length(vaf_vidNeg_T[tid])).*ones(length(vaf_vidNeg_T[tid])), 2*vaf_vidNeg_T[tid],
        color=colorCycle(2),
    )
    boxplot!(
        (xLoc+1)*ones(length(vaf_vidSat_T[tid])),
        2*vaf_vidSat_T[tid],
        show_outliers=false,
        color=(colorCycle(3),0.5),
    )
    scatter!(
        (xLoc+1) .+ 0.1 .*randn(length(vaf_vidSat_T[tid])).*ones(length(vaf_vidSat_T[tid])), 2*vaf_vidSat_T[tid],
        color=colorCycle(3),
    )
end

display(figSizeBox)
SAVEFIG && save("Figures/cloneSizesBoxplot.png",figSizeBox)

#endregion

## --------------------------------------------------
#region - Growth rates with time per gene
using Statistics

datapointThreshold = 20
#list of genes
_gene = DataStructuring.dropdupes(dfVid[!,:gene])
for gene in _gene
    _geneMask = dfVid[!,:gene] .== gene
    if sum(_geneMask .&& _fitMask) < datapointThreshold continue end
    _t = [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_fitMask.&&_geneMask,:])]
    _x = [γ for γ in dfVid[_fitMask.&&_geneMask, :γ]]
    cor = corspearman(_t, _x)
    rscale = 0.8
    figGRTimeGene = Figure(
        resolution=(rscale*700,rscale*600),
        fontsize=24,
    )
    grThreshold = 11
    # DataStructuring.testFit!(dfVid, 0.99)
    markersize = 13
    Axis(
        figGRTimeGene[1,1],
        xlabel="age at measurement",
        ylabel="fitted logistc growth rate",
        title="gene: "*gene*"; cor="*string(round(cor, digits=2)),
    )
    scatter!(
        _t, _x;
        markersize,
    )
    # ylims!(-1,1)
    # axislegend(position=:rt)
    display(figGRTimeGene)
    figname="gene_"*gene*".png"
    save("./Figures/DataFits/LogisticFitsPerGene/"*figname, figGRTimeGene)
end
#endregion

## --------------------------------------------------
#region - Fitted Growth rates with time Alternative plotstyle

rscale = 0.6
figGRTime = Figure(
    resolution=(rscale*800,rscale*600),
    fontsize=22,
)
grThreshold = 11
# DataStructuring.testFit!(dfVid, 0.99)

Axis(
    figGRTime[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
)
markersize=6
color=(:black,0)
strokewidth=1
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_fitMask.&&_posMask,:])],
    dfVid[_fitMask.&&_posMask, :γ];
    markersize,
    color,
    strokewidth,
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_fitMask.&&_negMask,:])],
    dfVid[_fitMask.&&_negMask, :γ];
    markersize,
    color,
    strokewidth,
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_fitMask.&&_satMask,:])],
    dfVid[_fitMask.&&_satMask, :γ];
    markersize,
    color,
    strokewidth,
)
# ylims!(-1,1)
display(figGRTime)

figname="dataGrowthRatesTime.png"
figloc="./Figures/Poster/"
# save(figloc*figname, figGRTime, px_per_unit=4)
#endregion

## --------------------------------------------------
#region - mean growth rates over time
bins = 19
_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask,:])], dfVid[_fitMask, :γ]; bins)
fig3a = Figure(
    # resolution=(800,600),
    # fontsize=22,
    resolution=(rescale*750,rescale*600),
    fontsize=26,
)
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate",
    # title="number of bins: "*string(bins)
)

ρMeanGR = corspearman(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
scatterlines!(_tBin, λAv_tBin,
    color=:black,
    label="spearman correlation: "*string(round(ρMeanGR,digits=2)),
)
axislegend(position=:rt)
display(fig3a)
figname="binnedAvFitnessLogit.png"
figloc="./Figures/ManuscriptDrafts/Figure1/"
# save(figloc*figname, fig3a)
#endregion

## --------------------------------------------------
#region - trajectory types + mean growth rate single plot

using HypothesisTests
using GLM

rscale = 1
bins = 19
measureAge_vid = [mean(_t) for _t in dfVid[_fitMask,:_t]]
_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfVid[_fitMask,:vid],
    bins=bins
)
_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nNeg_t = [sum((vid->in(vid,vidT_id)).(dfVid[_negMask .&& _fitMask,:vid])) for vidT_id in vid_id_T]
nVars_t = length.(vid_id_T)
sampleWeights_t = nVars_t ./ sum(nVars_t)
# ρDecreasing = corspearman(
#     _t[.!isnan.(nNeg_t./nVars_t)],
#     (nNeg_t./nVars_t)[.!isnan.(nNeg_t./nVars_t)]
# )
nonEmpty_t = .!isnan.(nNeg_t./nVars_t)
ρDecreasing = cor(
    _t[nonEmpty_t],
    (nNeg_t./nVars_t)[nonEmpty_t]
)

pDecreasing = pvalue(CorrelationTest(
    _t[nonEmpty_t],
    (nNeg_t./nVars_t)[nonEmpty_t]
    ))
ρDecreasingPearson = cor(
    _t[nonEmpty_t],
    (nNeg_t./nVars_t)[nonEmpty_t]
)

# function weightedCor(x_s, y_s, weights_s)
#     xMean = sum(weights_s .* x_s) / sum(weights_s)
#     yMean = sum(weights_s .* y_s) / sum(weights_s)
    
#     covWeighted = sum(weights_s .* (x_s .- xMean) .* (y_s .- yMean)) / sum(weights_s)
#     xVarWeighted = sum(weights_s .* (x_s .- xMean).^2) / sum(weights_s)
#     yVarWeighted = sum(weights_s .* (y_s .- yMean).^2) / sum(weights_s)
    
#     corrWeighted = covWeighted / sqrt(xVarWeighted * yVarWeighted)
    
#     return corrWeighted
# end

# weightedCor(
#     _t[nonEmpty_t],
#     (nNeg_t./nVars_t)[nonEmpty_t],
#     sampleWeights_t[nonEmpty_t],
# )

println("ρ: ", ρDecreasing)
println("p: ", pDecreasing)
# Perform linear regression
model = lm(@formula(y ~ x), DataFrame(x=_t[nonEmpty_t], y=(nNeg_t./nVars_t)[nonEmpty_t]))
decInt = coef(model)[1]
decSlope = coef(model)[2]

figTrajTime = Figure(
    fontsize=26,
    resolution=(rscale*750,rscale*600)
)
Axis(figTrajTime[2,1],
    ylabel="fraction of trajectories\ndecreasing",
    # xlabel="mean age at measurement",
    yticks=[0,.2,.4,.6,.8],
    xticklabelsvisible=false,
    xticksvisible=false,
)
ylims!(-0.03,0.7)
scatter!(
    # _t, nPos_t./nVars_t,
    _t, nNeg_t ./ nVars_t,
    color=:black,
)
lines!(
    _t, (t->decInt+decSlope*t).(_t),
    color=:black,
    label="correlation: "*string(round(ρDecreasing,digits=2))*" (p="*string(round(pDecreasing,digits=3))*")",
)
axislegend(position=:lt,labelsize=20, framevisible=false)

_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask,:])], dfVid[_fitMask, :γ]; bins)
Axis(
    figTrajTime[1,1],
    # yaxisposition=:right,
    ylabel="mean logistic\nfitness",
    xlabel="mean age at measurement",
    # title="number of bins: "*string(bins)
)
ylims!(0,0.15)

# _tAll = [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfVid[_fitMask,:])]
# gr_tAll = [γ for γ in dfVid[_fitMask, :γ]]
# corspearman(_tAll, gr_tAll)
# pvalue(CorrelationTest(_tAll, gr_tAll))

nonEmpty_t = .!isnan.(λAv_tBin)
# ρFit = corspearman(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
ρFit = cor(_tBin[nonEmpty_t], λAv_tBin[nonEmpty_t])
# pFit = pvalue(CorrelationTest(
#     _tBin[.!isnan.(λAv_tBin)],
#     λAv_tBin[.!isnan.(λAv_tBin)]
#     ))
pFit = pvalue(CorrelationTest(
            _tBin[nonEmpty_t],
            λAv_tBin[nonEmpty_t]
    )
)

# weightedCor(
#     _t[nonEmpty_t],
#     λAv_tBin[nonEmpty_t],
#     nVars_t[nonEmpty_t],
# )

println("ρ: ", ρFit)
println("p: ", pFit)
ρPearsonFit = cor(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
# Perform linear regression
model = lm(@formula(y ~ x), DataFrame(x=_tBin[.!isnan.(λAv_tBin)], y=λAv_tBin[.!isnan.(λAv_tBin)]))
fitInt = coef(model)[1]
fitSlope = coef(model)[2]
scatter!(_tBin, λAv_tBin,
    color=:black,
)
lines!(
    _t, (t->fitInt+fitSlope*t).(_t),
    color=:black,
    label="correlation: "*string(round(ρFit,digits=2))*" (p="*string(round(pFit,digits=2))*")",
)
axislegend(position=:rt, labelsize=20, framevisible=false)
display(figTrajTime)

figname="trajectoryTypes_MeanFitness_Time.svg"
figloc="./Figures/ManuscriptDrafts/"
# save(figloc*figname, figTrajTime)

#endregion

## ---------------------------------------------------------------------
#region - Spearman correlation of mean fitness over time per number of bins

mean(dfVid[_fitMask,:γ])
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, λAv_tBin = DataStructuring.meanVariantParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_fitMask,:])], dfVid[_fitMask, :γ]; bins)
    corS_nBins[i] = corspearman(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
end
fig3 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(
    fig3[1,1],
    xlabel="number of bins",
    ylabel="spearman correlation of mean fitness with time"
)
scatter!(nBins, corS_nBins, color=:black)
lines!(nBins, zeros(length(nBins)),color=:red)
display(fig3)

#endregion

## ----------------------------------------
#region - Number of variant types per time bin

_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]
_tEdges, fitType_vid_T = DataStructuring.binVariantParamTime(
    [mean(dfVCur[:_t]) for dfVCur in eachrow(dfVid[_timeMask,:])],
    dfVid[_timeMask,:fitType],
    bins=23
)

nConstant_t = [sum(fitType_vid.=="constant")/length(fitType_vid) for fitType_vid in fitType_vid_T]
nNegative_t = [sum(fitType_vid.=="negative logistic")/length(fitType_vid) for fitType_vid in fitType_vid_T]
nIncreasing_t = [sum(fitType_vid.=="positive logistic")/length(fitType_vid) for fitType_vid in fitType_vid_T]
nConstOrNegative_t = nConstant_t + nNegative_t

_t = (_tEdges[1:end-1].+_tEdges[2:end])./2
fig3 = Figure()
Axis(fig3[1,1],
    xlabel="mean age at measurement",
    ylabel="proportion of trajectories"
)
# scatterlines!(_t,nConstant_t,label="saturating")
# scatterlines!(_t,nNegative_t,label="decreasing")
scatterlines!(_t,nIncreasing_t,label="")
# scatterlines!(_t,nIncreasing_t,label="increasing")
axislegend(position=:lt)
display(fig3)
#endregion

## ----------------------------------------
#region - Spearman correlation of variant types with time

nBins = 10:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tEdges, fitType_vid_T = DataStructuring.binVariantParamTime(
    [mean(dfVCur[:_t]) for dfVCur in eachrow(dfVid[!,:])],
    dfVid[!,:fitType],
    bins=bins
    )
    _t = (_tEdges[1:end-1].+_tEdges[2:end])./2
    nConstant_t = [sum(fitType_vid.=="constant")/length(fitType_vid) for fitType_vid in fitType_vid_T]
    corS_nBins[i] = corspearman(_t[.!isnan.(nConstant_t)], nConstant_t[.!isnan.(nConstant_t)])
end

fig3 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(
    fig3[1,1],
    xlabel="number of bins",
    ylabel="spearman correlation of mean fitness with time"
)
scatter!(nBins, corS_nBins, color=:black)
lines!(nBins, zeros(length(nBins)),color=:red)
display(fig3)
#endregion

## ========================================
#region = Intermezzo: gamma dists

using Distributions


parA = (
    s = 0.1,
    σ = 0.03,
)
distA = Gamma(parA.s^2/parA.σ^2 , parA.σ^2/parA.s)

distA_x = [pdf(distA, x) for x in _x]

# lines(_x, distA_x)
println(mean(distA))
println(std(distA))
##

_x = 0:0.001:0.4
s = 0.1
_par = [0.02, 0.05, 0.1, 0.2,]
density_x_par = Array{Float64, 2}(undef, length(_x), length(_par))
for (i,σ) in enumerate(_par)
    density_x_par[:,i] = (x -> pdf(Gamma(s^2/σ^2, σ^2/s),x)).(_x)
end

figGammaDists = Figure(resolution=(680,550),fontsize=26)
Axis(figGammaDists[1,1],
    title="s = "*string(s)
)
for (i,σ) in enumerate(_par)
    lines!(_x, density_x_par[:, i],
        label="σ = "*string(σ)
    )
end
axislegend()
display(figGammaDists)

#endregion

## ----------------------------------------
#region - Correlation Negative fitness vs initial frequency
using HypothesisTests

vafIn_vidNeg = [dfVidCur[:vaf_t][1] for dfVidCur in eachrow(dfVid[_negMask,:])]
s_vidNeg = [dfVidCur[:γ] for dfVidCur in eachrow(dfVid[_negMask,:])]

figNegFitCor = Figure(fotnsize=24)
Axis(figNegFitCor[1,1],
    yscale=log10,
    ylabel="Clone size",
    xlabel="(negative) fitness strength"
)
# xlims!(-.01,0.8)
xlims!(0.001,0.8)
# ylims!(0.002,0.35)
# ylims!(-0.01,0.15)
# scatter!(vafIn_vidNeg, -s_vidNeg, color=colorCycle(2))
scatter!(-s_vidNeg, vafIn_vidNeg, color=colorCycle(2))

modelFitFreq = lm(@formula(y ~ x), DataFrame(x=-s_vidNeg, y=log.(vafIn_vidNeg)))
decInt = coef(modelFitFreq)[1]
decSlope = coef(modelFitFreq)[2]
_s = 0:0.01:0.8
linReg_s = (s->decInt+decSlope*s).(_s)
# lines!(
#     _s[linReg_s.>0], linReg_s[linReg_s.>0],
#     color=:black,
# )
lines!(
    _s, exp.(linReg_s),
    color=:black,
)

display(figNegFitCor)

cor(-s_vidNeg, log.(vafIn_vidNeg))
pvalue(CorrelationTest(
    -s_vidNeg,
    log.(vafIn_vidNeg),
    ))


cor(vafIn_vidNeg,-s_vidNeg)
corspearman(vafIn_vidNeg,-s_vidNeg)
pvalue(CorrelationTest(
    -s_vidNeg,
    vafIn_vidNeg,
    ))

#endregion


