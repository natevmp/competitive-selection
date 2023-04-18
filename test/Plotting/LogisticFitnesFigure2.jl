using CSV, DataFrames, Statistics, StatsBase
# using Optim, Distributions
# include("../../src/competitiveSelection.jl")
# using .CompetitiveSelection
include("../../src/dataStructuring.jl")
using ProgressMeter
using .DataStructuring

## ====================== Data Analysis =====================

df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
Nf = 2E5
dfVid = DataStructuring.analyseVariants(df, Nf; adjustedLogit=false)


## ===================== Simulations =======================

include("../../src/competitiveSelection.jl")
using .CompetitiveSelection

##

paramsPop = Dict(
    :N => Nf,
    :α => 1,
    # :s => 0.12,
    :s => 0.05,
    :σ => 0.25,
    :T => 105,
    # :μ => 0.14,
    :μ => 0.24,
    :sType => "gamma",
    # :sType => "free",
    :q => 0.3,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 5000

##
@time solEns = CompetitiveSelection.evolvePopSim(paramsPop, runs=runs)
@time fit_tλ_cid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop, tMeasure=(50,90), timeLimit=0.4)


## ======================= Plotting ========================
using CairoMakie

## ----------- Growth curve per variant ----------
rscale=0.7
res = (rscale*600,rscale*500)
fsize=16

@showprogress for (i,dfVCur) in enumerate(eachrow(dfVid[1:10,:]))
    fig1 = Figure(resolution=res, fontsize=fsize)
    ax = Axis(
        fig1[1,1],
        title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene],
    )
    scatter!(dfVCur[:_t], dfVCur[:vaf_t],
        linestyle=:solid,
        label=dfVCur[:gene],
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label="s = "*string(round(dfVCur[:λ],digits=2))*"\nz = "*string(round(dfVCur[:z],digits=2)),
    )
    rangebars!(
        dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
        whiskerwidth = 20,
        cycle = [:color],
    )
    legendpos = ( dfVCur[:ρ]>0 ? :lt : :rt )
    axislegend(position=legendpos)
    # figname = "dataFit_vid$i.png"
    # save("Figures/DataFits/PerVariant/"*figname, fig1)
    display(fig1)
end

## ---------- All fitted fitness values -----------

fig2 = Figure(
    resolution=(700,600),
    fontsize=22,
)
λThresh=2
# λThresh=9
_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]
_λMask = [dfVCur[:λ]<λThresh for dfVCur in eachrow(dfVid)]
# ρS = corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_timeMask .& _λMask,:])], dfVid[_timeMask .& _λMask, :λ])
ρS = corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])], dfVid[!, :λ])
Axis(
    fig2[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
    title="Spearman correlation: ρ="*string(round(ρS,digits=3))
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])],
    dfVid[_λMask, :λ],
    color=:black,
    markersize=5,
)
# ylims!(-1,1)
display(fig2)

## ---------- Fitness distribution -----------

fig2b = Figure()
Axis(
    fig2b[1,1],
    xlabel="inferred growth rate",
    ylabel="density of variants",
    title="mean growth rate: "*string( round(mean(dfVid[_λMask,:λ]), digits=3) )
)
# hist!(dfVid[_λMask, :λ], normalization=:pdf, bins=25)
hist!(dfVid[_λMask, :λ], normalization=:pdf, bins=25)
display(fig2b)

## ---------- Mean fitness over time ---------

bins = 15
_tBin, λAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])], dfVid[_λMask, :λ]; bins)
fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate",
)
scatterlines!(_tBin, λAv_tBin)
display(fig3a)

## ---------- Fitness distribution over time -----------

bins = 8
_tBinEdge, λ_vid_TBin = DataStructuring.binVariantsTime(dfVid[_λMask,:]; bins)
_tBin = _tBinEdge[1:end-1] .+ (_tBinEdge[2] - _tBinEdge[1])

for (i,t) in enumerate(_tBin)
    fig = Figure()
    Axis(fig[1,1])
    hist!(λ_vid_TBin[i], normalization=:pdf, bins=15)
    xlims!(-1,1)
    display(fig)
end

## ---------- Binned Means ------------

mean(dfVid[_λMask,:λ])
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, λAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])], dfVid[_λMask, :λ]; bins)
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
scatter!(nBins, corS_nBins)
display(fig3)

## ----------- Fitted fitness values simulations ------------

fig4 = Figure()
_λmask = fit_tλ_cid[2,:] .<= 1
ρSpSim = corspearman(fit_tλ_cid[1,_λmask], fit_tλ_cid[2,_λmask])
Axis(
    fig4[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
    title="Spearman correlation: ρ="*string(round(ρSpSim,digits=3))
)
scatter!(
    fit_tλ_cid[1,_λmask], fit_tλ_cid[2,_λmask],
    color=:black,
    markersize=5,
)
xlims!(50,100)
ylims!(-1,1)
display(fig4)

## ---------- Mean fitness over time Sims ---------

bins = 15
_tBin, λAv_tBin = meanVarFitParamTime(fit_tλ_cid[1,_λmask], fit_tλ_cid[2,_λmask]; bins)
fig4a = Figure()
Axis(
    fig4a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate"
)
scatterlines!(_tBin, λAv_tBin)
display(fig4a)

## ---------- Binned Means Sims ----------

# mean(dfVid[!,:λ])
nBins = 6:30
corSpSims_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    @views _tBin, λAv_tBin = meanVarFitParamTime(fit_tλ_cid[1,_λmask],fit_tλ_cid[2,_λmask]; bins)
    corSpSims_nBins[i] = corspearman(_tBin, λAv_tBin)
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
scatter!(nBins, corSpSims_nBins)
display(fig3)

## ---------- All fitted Z-values -----------

fig2 = Figure(
    resolution=(700,600),
    fontsize=22,
)
_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]
_λMask = [dfVCur[:λ]<λThresh for dfVCur in eachrow(dfVid)]
# ρS = corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_timeMask .& _λMask,:])], dfVid[_timeMask .& _λMask, :λ])
ρS = corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])], dfVid[!, :z])
Axis(
    fig2[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc competition factor Z",
    title="Spearman correlation: ρ="*string(round(ρS,digits=3))
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])],
    dfVid[!, :z],
    color=:black,
    markersize=5,
)
# ylims!(-1,1)
display(fig2)

##
mean(dfVid[!,:z])
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, zAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)], dfVid[!, :z]; bins)
    corS_nBins[i] = corspearman(_tBin[.!isnan.(zAv_tBin)], zAv_tBin[.!isnan.(zAv_tBin)])
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
scatter!(nBins, corS_nBins)
display(fig3)

