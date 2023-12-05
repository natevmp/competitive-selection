include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using CSV, DataFrames
using Statistics
include("../src/mullerPlot.jl")
using .MullerPlot

## ========================================
#region = parameters + Run sims

params = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 100000,
    :α => 0.6,
    :sType => "gamma",
    :q => 1,
    :μ => 8,
    :s => 0.1,
    :σ => 0.04,
)
ctrlParams = Dict{}(
    :errorQuantile => 0.98,
)

@time solEns, simArgs = CompetitiveSelection.evolvePopSim(
    params;
    runs=1000,
    noDiffusion=false
)

solId = 1

#endregion

## ========================================
#region = Figures

using CairoMakie
rscale = 0.8
linewidth = 1.5

#endregion

## ----------------------------------------
#region - Russian doll statistics
using CairoMakie

sizeThreshold = 0.001
_t = range(40,90)
meanSizeDolls_t = Vector{Float64}(undef, length(_t))
meanNumberDolls_t = Vector{Float64}(undef, length(_t))
propDollCarrying_t = Vector{Float64}(undef, length(_t))
for (tid,t) in enumerate(_t)
    nDollFree = 0
    nDollCarrying = 0
    sizeDolls_vid1 = []
    nDolls_vid1 = []
    for sid in eachindex(solEns)
        sol = solEns[sid]
        n_vid = sol(t)
        parentVid_vid = simArgs.parentId_vid[sid]
        visibleT_vid = n_vid .> sizeThreshold
        russianDoll_vid = parentVid_vid .!== 0
        nDolls = sum(visibleT_vid .&& russianDoll_vid)
        if nDolls > 0
            nDollCarrying += 1
        else
            nDollFree += 1
        end
        push!(sizeDolls_vid1, sum(n_vid[visibleT_vid .&& russianDoll_vid]))
        push!(nDolls_vid1, nDolls)
    end
    meanSizeDolls = isempty(sizeDolls_vid1) ? 0 : mean(sizeDolls_vid1)
    propDollCarrying_t[tid] = nDollCarrying / (nDollCarrying + nDollFree)
    meanSizeDolls_t[tid] = meanSizeDolls
    meanNumberDolls_t[tid] = mean(nDolls_vid1)
end

fig = Figure(fontsize=28)
ax1 = Axis(
    fig[1,1],
    ylabel="total size of clones"
)
linAvCl = lines!(
    _t, meanSizeDolls_t,
    color=:black,
    label="average size of clones"
)
ax2 = Axis(
    fig[1,1],
    ylabel="number of clones",
    xlabel="age",
    yaxisposition=:right
)
linNCl = lines!(
    _t, meanNumberDolls_t,
    color=:black,
    linestyle=:dash,
    label="av number of clones"
)
ax2 = Axis(
    fig[1,1],
    ylabel="number of clones",
    xlabel="age",
    yaxisposition=:right
)
lines!(_t, propDollCarrying_t, color=Cycled(2)) 
axislegend(ax1, [linAvCl,linNCl], ["mean size of\nmulti-mutant clones", "mean number of\nmulti-mutant clones"], position=:lt)
display(fig)

figname="multimutantCloneStats.svg"
# save("Figures/ManuscriptDrafts/"*figname, fig)
#endregion

## ========================================
#region = Muller plot

include("../src/mullerPlot.jl")
using .MullerPlot
solId = 6
sol = solEns[solId]
println("total variant fraction at end: ", sum(sol.u[end]))

lines([sum(sol(t)) for t in 1:100], figure=(fontsize=24,), axis=(xlabel="age",ylabel="proportion of HSC population\nwith at least 1 driver"))

simArgs.parentId_vid[solId]
parentVid_vid = simArgs.parentId_vid[solId]
russianDoll_vid = parentVid_vid .!== 0
vid_child_Vid0 = MullerPlot.buildFamilyArray(parentVid_vid)
sizeThreshold = 0.001
n_vid = sol.u[end]
visible_vid = n_vid .> sizeThreshold

sum(visible_vid)
visible_vid0 = prepend!(copy(visible_vid), [true,])
vid_child_Vid0[visible_vid0]
println(
    "total size of russian doll clones: ", sum(n_vid[russianDoll_vid])
)
println("number of russian doll clones: ", sum(.!isempty.(vid_child_Vid0[visible_vid0]))-1)

_t = range(0,100,length=length(sol))

n_t_vid0 = Array{Float64,2}(undef, length(sol), length(sol[1])+1)
for i in 1:length(sol[1])
    n_t_vid0[:,1+i] = [solT[i] for solT in sol.u]
end
n_t_vid0[:,1] = [1-sum(n_t_vid0[t,2:end]) for t in 1:size(n_t_vid0,1)]
n_t_vid0[findall(x->x<0, n_t_vid0)] .= 0

xL_t_vid, xU_t_vid = MullerPlot.mullerBounds(n_t_vid0, vid_child_Vid0, visible_vid0)

size(xL_t_vid)

tBirth_vid = MullerPlot.variantBirthTimes(n_t_vid0[:,visible_vid0])[2:end]

using CairoMakie
using ColorSchemes
colorMap = ColorSchemes.leonardo
colorMap = ColorSchemes.tableau_20
scale = 0.8
lw = 1.5
# set_theme!()
fig1 = Figure(resolution=(scale*900,scale*800), fontsize=25)
axtop = Axis(
    fig1[1,1],
    backgroundcolor=:grey65,
    # palette=(patchcolor=cgrad(:tab20),),
    xgridvisible=false,
    ygridvisible=false,
    yticksvisible=false,
    yticklabelsvisible=false,
    xticks=[0,25,50,75,100],
    xlabel="time",
    ylabel="size of clone"
)
# hidedecorations!(axtop)
for i in 1:(sum(visible_vid0)-1)
    band!(
        _t, xL_t_vid[:,1+i], xU_t_vid[:,1+i],
        show_grid=false,
        # color=colorMap[i],
    )
    scatter!([tBirth_vid[i],], [xU_t_vid[Int(tBirth_vid[i]), 1+i]+0.015,],
    # color=colorMap[i],
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
# save("Figures/ManuscriptDrafts/exampleSim.svg", fig1)

#endregion

## ----------------------------------------
#region - Average total variant fraction

detectThresh=0.01
_t = 0:(length(solEns[1])-1)
nTotAv_t = Vector{Float64}(undef, length(_t))
nClones_t = Vector{Float64}(undef, length(_t))
nClonesThresh_t = Vector{Float64}(undef, length(_t))
for tInd in eachindex(_t)
    nTot = 0
    nClones = 0
    nClonesThresh = 0
    for sol in solEns
        nTot += sum(sol[tInd])
        nClones += sum(sol[tInd].>0)
        nClonesThresh += sum(sol[tInd].>detectThresh)
    end
    nTotAv_t[tInd] = nTot / length(solEns)
    nClones_t[tInd] = nClones / length(solEns)
    nClonesThresh_t[tInd] = nClonesThresh / length(solEns)
end


nTotTheory_t = (t->Theory.compCoverage(t, params[:α], params[:s], params[:μ], params[:N])).(_t)

fig = Figure(fontsize=30)
Axis(fig[1,1],
    xlabel="time",
    # ylabel="number of clones",
    ylabel="clonal fraction",
)

# lines!(nTotTheory_t)
# lines!(nClones_t, label="clones alive")
# lines!(nClonesThresh_t, label="clones above threshold")
lines!(nTotAv_t,
    # label="sMean="*string(params[:s])*"; sStd="*string(params[:σ])
    # label="N="*string(params[:N])*"; α="*string(params[:α])
    label="α="*string(params[:α])*"; sMean="*string(params[:s]),
)
ylims!(0,1)
axislegend(position=:lt)
display(fig)

#endregion

## ----------------------------------------
#region - Barplot variant sizes after sampling over time
using StatsBase

function cellSampler(sol, S::Int, t::Real)
    weight_vid0 = prepend!(sol(t), 1-sum(sol(t)))
    vid_sid = sample(range(0,length(weight_vid0)-1), Weights(weight_vid0), S)
    _vidS = sort!(unique(vid_sid))
    f_vidS = Vector{Float64}(undef, length(_vidS))
    for (i,vid) in enumerate(_vidS)
        f_vidS[i] = count(el->el==vid, vid_sid) / S
    end
    return _vidS, f_vidS
end

# _tt = [40,50,60,70,80]
_tt = 0:10:100

simID = solId
sol = solEns[simID]
f_vidS_Tt = Vector{Vector{Float64}}(undef, length(_tt))
for (tid,t) in enumerate(_tt)
    _vidS, f_vidS = cellSampler(sol, 300, t)
    vidIn = findfirst(_vidS.!=0)
    if isnothing(vidIn)
        f_vidS_Tt[tid] = Int[]
        continue
    end
    f_vidS_Tt[tid] = f_vidS[vidIn:end]
end

##


detected_vidS_Tt = (f_vidS->f_vidS.>detectThresh).(f_vidS_Tt)
println(length(f_vidS_Tt[end][detected_vidS_Tt[end]]))
f_vidS_Tt
rscale=0.8
fig = Figure(fontsize=30, resolution=(rscale*900,rscale*800))
Axis(fig[1,1],
    xticks=[0,25,50,75,100],
    xlabel="time",
    ylabel="size of clone"
)
xlims!(0,105)
# position = cat(1:nSims, 2*nSims+1:(2*nSims+nSims), 4*nSims+1:(4*nSims+nSims), 6*nSims+1:(6*nSims+nSims), 8*nSims+1:(8*nSims+nSims), dims=1)

positions = cat([_tt[tid]*ones(length(f_vidS_Tt[tid][detected_vidS_Tt[tid]])) for tid in eachindex(_tt)]..., dims=1)
heights = cat([f_vidS_Tt[tid][detected_vidS_Tt[tid]] for tid in eachindex(_tt)]..., dims=1)
stack = cat( [range(1,length(f_vidS_Tt[tid][detected_vidS_Tt[tid]])) for tid in eachindex(_tt)]..., dims=1 )
barplot!(
    positions, heights,
    stack=stack,
    color=stack,
)

display(fig)

#endregion

## ----------------------------------------
#region - example trajectories

# solID = 6
sol = solEns[solID]
size_vid_T = solEns[solID].u
fig1 = Figure(resolution=(rscale*800, rscale*600), fontsize=24)
axbottom = Axis(
    fig1[1,1], 
    xlabel="time",
    ylabel="size of clone",
    xticks=[0,25,50,75,100],
    # xticksvisible=false,
    # xticklabelsvisible=false,
    # ygridvisible=true,
    # xgridvisible=true,
)
ylims!(0,0.4)
xlims!(0,100)
# hidedecorations!(axbottom, grid=false)
for i in 1:length(size_vid_T[1])
    sizeV_t = [size_vid[i] for size_vid in size_vid_T]
    lines!(
        sizeV_t,
        # color=Cycled(i),
        linewidth=linewidth,
    )
end

# ylims!(0,0.5)
# xlims!(0,90)
display(fig1)

# save("./Figures/ManuscriptDrafts/Figure 2/simulationExampleTrajectories.png", fig1, px_per_unit=2)
#endregion

## ========================================
#region = fit simulations to logistic growth with sampling

# @time dfSimsVid = CompetitiveSelection.fitSamplesGrowth(solEns, params, tMeasure=(60,85); errorQuantile=ctrlParams[:errorQuantile])


# run `simulationsExamplePlots.jl` first to run this
bins = 19
measureAge_vid = [_t[1] for _t in dfVid[_fitMask,:_t]]
_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfVid[_fitMask,:vid],
    bins=bins
)
length(_tEdges)
length(nVars_t)
_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nVars_t = length.(vid_id_T)
@time dfSimsVid = CompetitiveSelection.fitSamplesGrowth(solEns, params, tMeasure=(_tEdges, nVars_t); errorQuantile=ctrlParams[:errorQuantile])
dfSimsVid.vid = 1:size(dfSimsVid,1)
_posMaskSims = dfSimsVid[!,:fitType] .== "positive logistic"
_negMaskSims = dfSimsVid[!,:fitType] .== "negative logistic"
_satMaskSims = dfSimsVid[!,:fitType] .== "constant"
_timeMaskSims = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfSimsVid)]
grThreshold = 1
_λMaskSims = [dfVCur[:γ]<grThreshold for dfVCur in eachrow(dfSimsVid)]
_fitMaskSims = dfSimsVid[!,:goodFit]


# dfSimsVid[!,:γ]
# size(dfSimsVid[!,:])

#endregion

## ----------------------------------------
#region - fitted logistic fitnesses with time
sum(_negMaskSims)
# sum(_satMaskSims)
rscale = 0.8
figGRTime = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=24,
)
# DataStructuring.testFit!(dfVid, 0.99)
markersize=5
Axis(
    figGRTime[1,1],
    # title="Simulations",
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_posMaskSims,:])],
    dfSimsVid[_fitMaskSims.&&_posMaskSims, :γ];
    markersize,
)
scatter!(
    [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_negMaskSims,:])],
    dfSimsVid[_fitMaskSims.&&_negMaskSims, :γ];
    markersize,
)
# scatter!(
#     [mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_satMaskSims,:])],
#     dfSimsVid[_fitMaskSims.&&_satMaskSims, :γ];
#     markersize,
# )
ylims!(-0.4,0.8)
# xlims!()
display(figGRTime)
figname="SimsGrowthRatesTime.png"
figloc="./Figures/ManuscriptDrafts/"
# save(figloc*figname, figGRTime)

#endregion

## ----------------------------------------
#region - size of saturating clones

#endregion

## ----------------------------------------
#region - mean growth rates over time
include("../src/dataStructuring.jl")
using .DataStructuring

bins = 26
_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_λMaskSims,:])], dfSimsVid[_fitMaskSims.&&_λMaskSims, :γ]; bins)
fig3a = Figure(
    resolution=(800,600),
    fontsize=22,
)
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate",
    title="number of bins: "*string(bins)
)
scatterlines!(_tBin, λAv_tBin)
display(fig3a)

#endregion

## ----------------------------------------
#region - average fitness + decreasing trajectories w time

using GLM
using Statistics, StatsBase

dfSimsVid.vid = 1:nrow(dfSimsVid)
rscale = 1
bins = 19
measureAge_vid = [mean(_t) for _t in dfSimsVid[_fitMaskSims,:_t]]
_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfSimsVid[_fitMaskSims,:vid],
    bins=bins
)
_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nNeg_t = [sum((vid->in(vid,vidT_id)).(dfSimsVid[_negMaskSims .&& _fitMaskSims,:vid])) for vidT_id in vid_id_T]
nVars_t = [sum((vid->in(vid,vid_t)).(dfSimsVid[_fitMaskSims,:vid])) for vid_t in vid_id_T]
ρDecreasing = corspearman(
    _t[.!isnan.(nNeg_t./nVars_t)],
    (nNeg_t./nVars_t)[.!isnan.(nNeg_t./nVars_t)]
)
ρDecreasingPearson = cor(
    _t[.!isnan.(nNeg_t./nVars_t)],
    (nNeg_t./nVars_t)[.!isnan.(nNeg_t./nVars_t)]
)
# Perform linear regression
model = lm(@formula(y ~ x), DataFrame(x=_t[.!isnan.(nNeg_t./nVars_t)], y=(nNeg_t./nVars_t)[.!isnan.(nNeg_t./nVars_t)]))
decInt = coef(model)[1]
decSlope = coef(model)[2]

figTrajTime = Figure(
    fontsize=26,
    resolution=(rscale*750,rscale*600)
)
Axis(figTrajTime[1,1],
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
    label="Spearman correlation: "*string(round(ρDecreasing,digits=2)),
)
axislegend(position=:lt,labelsize=20, framevisible=false)

_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfSimsVid[_fitMaskSims,:])], dfSimsVid[_fitMaskSims, :γ]; bins)
Axis(
    figTrajTime[2,1],
    # yaxisposition=:right,
    ylabel="mean logistic\nfitness",
    xlabel="mean age at measurement",
    # title="number of bins: "*string(bins)
)
ylims!(0,0.15)
ρFit = corspearman(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
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
    label="Spearman correlation: "*string(round(ρFit,digits=2)),
)
axislegend(position=:rt, labelsize=20, framevisible=false)
display(figTrajTime)
#endregion

## ----------------------------------------
#region - Spearman correlation of mean fitness over time per number of bins

mean(dfSimsVid[!,:γ])
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, λAv_tBin = DataStructuring.meanVariantParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfSimsVid[!,:])], dfSimsVid[!, :γ]; bins)
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
#region - Fitted growt rates distribution

tRange = (70,75)
_ageAvMaskSims = (x-> x>tRange[1] && x<tRange[2]).(mean.(dfSimsVid[!,:_t]))
_ageAvMask = (x-> x>tRange[1] && x<tRange[2]).(mean.(dfSimsVid[!,:_t]))

rescale = 0.9
figFitnessDist = Figure(
    fontsize=26,
    resolution=(rescale*750,rescale*600)
)
Axis(
    figFitnessDist[1,1],
    title="simulations",
    xlabel="logistic fitness",
    ylabel="density of variants"
)
hist!(
    dfSimsVid[_fitMaskSims.&&_λMaskSims.&&_ageAvMaskSims,:γ],
    bins=30,
    normalization=:pdf,
    # strokewidth=2,
    # color=:grey45,
    # strokecolor=:grey45,
    label="sims",
)
hist!(
    dfSimsVid[_fitMaskSims.&&_ageAvMask,:γ],
    bins=30,
    normalization=:pdf,
    # strokewidth=2,
    # color=:grey45,
    # strokecolor=:grey45,
    label="data",
)
axislegend()
display(figFitnessDist)
figname="simsVariantFitnessDistribution.png"
figloc="./Figures/ManuscriptDrafts/"
# save(figloc*figname, figFitnessDist)
#endregion

## ========================================
#region = Sim Data masks

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

## --------------------------------------------------
#region - trajectory types + mean growth rate single plot

using HypothesisTests
using GLM

rscale = 1
bins = 15

measureAge_vid = [mean(_t) for _t in dfSimsVid[_fitMaskSims,:_t]]
_tEdges, vid_id_T = DataStructuring.binVariantParamTime(
    measureAge_vid,
    dfSimsVid[_fitMaskSims,:vid],
    bins=bins
)
_t = (_tEdges[2]-_tEdges[1])/2 .+_tEdges[1:end-1]
nNeg_t = [sum((vid->in(vid,vidT_id)).(dfSimsVid[_negMaskSims .&& _fitMaskSims,:vid])) for vidT_id in vid_id_T]
nVars_t = length.(vid_id_T)
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

_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([mean(dfVidCur[:_t]) for dfVidCur in eachrow(dfSimsVid[_fitMaskSims,:])], dfSimsVid[_fitMaskSims, :γ]; bins)
Axis(
    figTrajTime[1,1],
    # yaxisposition=:right,
    ylabel="mean logistic\nfitness",
    xlabel="mean age at measurement",
    # title="number of bins: "*string(bins)
)
ylims!(0,0.15)

nonEmpty_t = .!isnan.(λAv_tBin)
# ρFit = corspearman(_tBin[.!isnan.(λAv_tBin)], λAv_tBin[.!isnan.(λAv_tBin)])
ρFit = cor(_tBin[nonEmpty_t], λAv_tBin[nonEmpty_t])

pFit = pvalue(CorrelationTest(
            _tBin[nonEmpty_t],
            λAv_tBin[nonEmpty_t]
    )
)

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