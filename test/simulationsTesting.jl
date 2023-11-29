include("../src/competitiveSelection.jl")
using .CompetitiveSelection
using Statistics
using CairoMakie
using ElasticArrays, StatsBase

##

myTheme = Theme(fontsize = 20)
set_theme!(myTheme)

##

paramsPop = Dict(
    :N => Int(2E5),
    :α => 1,
    # :s => 0.12,
    :s => 0.15,
    :σ => 0.015,
    :T => 90,
    :μ => 3.5,
    # :μ => 0,
    :sType => "gamma",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 500

##
# @time sol = CompetitiveSelection.evolvePop(paramsPop)
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(paramsPop, runs=runs);
println("")

##
@time dfSimsVid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop, tMeasure=(50,90))
##

_posMaskSims = dfSimsVid[!,:fitType] .== "positive logistic"
_negMaskSims = dfSimsVid[!,:fitType] .== "negative logistic"
_satMaskSims = dfSimsVid[!,:fitType] .== "constant"
_timeMaskSims = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfSimsVid)]
grThreshold = 1
_λMaskSims = [dfVCur[:γ]<grThreshold for dfVCur in eachrow(dfSimsVid)]
_fitMaskSims = dfSimsVid[!,:goodFit]

##

rscale = 0.8
figGRTime = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=17,
)
# DataStructuring.testFit!(dfVid, 0.99)

Axis(
    figGRTime[1,1],
    title="Simulations",
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_posMaskSims,:])],
    dfSimsVid[_fitMaskSims.&&_posMaskSims, :γ],
    markersize=4,
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_negMaskSims,:])],
    dfSimsVid[_fitMaskSims.&&_negMaskSims, :γ],
    markersize=4,
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfSimsVid[_fitMaskSims.&&_satMaskSims,:])],
    dfSimsVid[_fitMaskSims.&&_satMaskSims, :γ],
    markersize=4,
)
ylims!(-0.4,1)
# xlims!()
display(figGRTime)
figname="SimsGrowthRatesTime.png"
figloc="./Figures/ManuscriptDrafts/"
# save(figloc*figname, figGRTime)

## --------- mean growth rates over time ----------
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

## -------- Spearman correlation of mean fitness over time per number of bins -----------

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

## ---------- Fitted growth rates distribution ---------

tRange = (70,75)
_ageAvMaskSims = (x-> x>tRange[1] && x<tRange[2]).(mean.(dfSimsVid[!,:_t]))
_ageAvMask = (x-> x>tRange[1] && x<tRange[2]).(mean.(dfVid[!,:_t]))

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
    dfVid[_fitMask.&&_ageAvMask,:γ],
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
