using CSV, DataFrames, Statistics, StatsBase
include("../src/dataStructuring.jl")
using ProgressMeter
using .DataStructuring
using CairoMakie

## ====================== Data Analysis =====================

variableSaturation = false

df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
Nf = 2E5
dfVid = DataStructuring.analyseVariants(df, Nf; competitionSaturate=variableSaturation)



## ====================== Plotting ===========================


## ----------- Growth curve per variant ----------
# rscale=0.7
# res = (rscale*600,rscale*500)
# fsize=16

# for (i,dfVCur) in enumerate(eachrow(dfVid[1:1,:]))
#     println("t0: ", dfVCur[:t0])
#     println("γ: ", dfVCur[:γ])
#     # println("x0: ", dfVCur[:t0])
#     fig1 = Figure(resolution=res, fontsize=fsize)
#     ax = Axis(
#         fig1[1,1],
#         title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene],
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


## ----------- Multiple data curves per variant in single plot ----------

rscale=0.7
res = (1300,1100)
fsize=18
fig0 = Figure(resolution=res, fontsize=fsize)
@showprogress for (i,dfVCur) in enumerate(eachrow(dfVid[1:12,:]))
    col = (i-1) % 3 + 1
    row = floor((i-1)/3) + 1 |> Int
    ax = Axis(
        fig0[row,col],
        title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene],
    )
    scatter!(dfVCur[:_t], dfVCur[:vaf_t],
        linestyle=:solid,
        label=dfVCur[:gene],
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label="s = "*string(round(dfVCur[:γ],digits=2))*"\n"
             *"k = "*string(round(dfVCur[:xF],digits=2)),
    )
    rangebars!(
        dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
        whiskerwidth = 20,
        cycle = [:color],
    )
    legendpos = ( dfVCur[:ρ]>0 ? :lt : :rt )
    axislegend(position=legendpos)
end
display(fig0)
saveFolder = "Figures/Presentations/"
fileName = "individualTrajectories_LogAdjustedFitted.png"
# save(saveFolder*fileName, fig0)

## ---------- Fitted Growth rates -----------

rscale = 0.8
fig2 = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=17,
)
grThreshold = 11
_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]
_λMask = [dfVCur[:γ]<grThreshold for dfVCur in eachrow(dfVid)]
ρS = corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])], dfVid[_λMask, :γ])
Axis(
    fig2[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
    # title="Spearman correlation: ρ="*string(round(ρS,digits=3))
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])],
    dfVid[_λMask, :γ],
    color=:black,
    markersize=4,
)
# ylims!(-1,1)
display(fig2)
saveFolder = "Figures/Presentations/"
fileName = "fitnessFits_logAdj.png"
# save(saveFolder*fileName, fig2)


## --------- mean growth rates over time ----------
bins = 28
_tBin, λAv_tBin = DataStructuring.meanVariantParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])], dfVid[_λMask, :γ]; bins)
fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate",
    title="number of bins: "*string(bins)
)
scatterlines!(_tBin, λAv_tBin)
display(fig3a)
saveFolder = "Figures/Presentations/"
fileName = "binnedAvFitnessLogit.png"
# save(saveFolder*fileName, fig3a)

## -------- Spearman correlation per number of bins -----------

mean(dfVid[_λMask,:γ])
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, λAv_tBin = DataStructuring.meanVariantParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[_λMask,:])], dfVid[_λMask, :γ]; bins)
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
saveFolder = "Figures/Presentations/"
fileName = "correlationPerBins.png"
# save(saveFolder*fileName, fig3)



## ----------- Fitted saturation levels -------------

rscale = 0.8
fig2 = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=17,
)
Axis(
    fig2[1,1],
    xlabel="age at measurement",
    ylabel="fitted saturation frequency",
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)],
    dfVid[!, :xF],
    color=:black,
    markersize=4,
)
# ylims!(-1,1)
display(fig2)

## ----------- mean saturation over time ---------
bins = 15
_tBin, xAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])], dfVid[!, :xF]; bins)
fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average saturation size",
)
scatterlines!(_tBin, xAv_tBin)
display(fig3a)

## -------------------- Plot curves of high fitness variants ------------

rscale=0.7
res = (rscale*600,rscale*500)
fsize=16

γHighThresh = 1.5
_γHMask = [dfVCur[:γ]>γHighThresh for dfVCur in eachrow(dfVid)]
for (i,dfVCur) in enumerate(eachrow(dfVid[_γHMask,:]))
    println("t0: ", dfVCur[:t0])
    println("γ: ", dfVCur[:γ])
    # println("x0: ", dfVCur[:t0])
    fig1 = Figure(resolution=res, fontsize=fsize)
    ax = Axis(
        fig1[1,1],
        title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene]*"; ρ "*string(round(dfVCur[:ρ],digits=3)),
    )
    scatter!(dfVCur[:_t], dfVCur[:vaf_t],
        linestyle=:solid,
        label=dfVCur[:gene],
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label="s = "*string(round(dfVCur[:γ],digits=2)),
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

## ----------------- Rank correlation per variant -------------

rscale = 0.8
fig2 = Figure(
    resolution=(rscale*700,rscale*600),
    fontsize=17,
)
grThreshold = 3
_timeMask = [dfVCur[:_t][1]<90 for dfVCur in eachrow(dfVid)]

Axis(
    fig2[1,1],
    xlabel="age at measurement",
    ylabel="measured rank correlation",
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])],
    dfVid[!, :ρ],
    color=:black,
    markersize=4,
)
display(fig2)

## ------------- Mean rank correlation over time ---------------
bins = 25
_tBin, ρAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])], dfVid[!, :ρ]; bins)
fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average rank correlation",
    title="number of bins: "*string(bins)
)
scatterlines!(_tBin, ρAv_tBin)
display(fig3a)

## ------------ Rank correlation of mean rank correlation with time ------------

nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, ρAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])], dfVid[!, :ρ]; bins)
    corS_nBins[i] = corspearman(_tBin[.!isnan.(ρAv_tBin)], ρAv_tBin[.!isnan.(ρAv_tBin)])
end
fig3 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(
    fig3[1,1],
    xlabel="number of bins",
    ylabel="spearman cor of mean rank cor"
)
scatter!(nBins, corS_nBins, color=:black)
lines!(nBins, zeros(length(nBins)),color=:red)
display(fig3)
saveFolder = "Figures/Presentations/"
fileName = "correlationPerBins.png"
# save(saveFolder*fileName, fig3)



## -------------- Mean variant size over time ----------------
import Base.Iterators.flatten
bins = 24
tThresh = 90
x_sid = [x for x in flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
tMask_sid = t_sid .< tThresh
_tBin, xAv_tBin = DataStructuring.meanVariantParamTime(t_sid[tMask_sid], x_sid[tMask_sid]; bins)

fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average size",
)
scatterlines!(_tBin, xAv_tBin)
display(fig3a)

## ------------- Patients with multiple clones ---------------

nCloneThreshold = 5
_pid = dropdupes(dfVid[!,:pid])

dfVid[!,:pid] .== _pid[1]

for pid in _pid
    pMask = dfVid[!,:pid] .== pid
    sum(pMask) < nCloneThreshold && continue
    fig = Figure()
    Axis(fig[1,1])
    for dfVCur in eachrow(dfVid[pMask, :])
        scatter!(dfVCur[:_t], dfVCur[:vaf_t],
            linestyle=:solid,
            # label=dfVCur[:gene],
        )
        lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
            linestyle=:dash,
            label="s = "*string(round(dfVCur[:γ],digits=2))*"\n"
                *"k = "*string(round(dfVCur[:xF],digits=2)),
        )
        rangebars!(
            dfVCur[:_t], dfVCur[:vafThQl_t], dfVCur[:vafThQu_t],
            whiskerwidth = 20,
            cycle = [:color],
        )
    end
    axislegend(position=:lt)
    display(fig)
end


## --------------- Size distribution ------------------
# include("../src/dataStructuring.jl")
# using .DataStructuring
tThresh = 90
x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
tMask_sid = t_sid .< tThresh
xMask_sid = x_sid .> 0
xtMask_sid = (t_sid.<tThresh) .& (x_sid.>0)
_T, _f, densV_f_T = DataStructuring.sizeDistBinTime(t_sid[xtMask_sid], x_sid[xtMask_sid]; tBins=4, fBins=50)

##
fig = Figure(
    resolution=(800,600),
    fontsize=22,
)
Axis(
    fig[1,1],
    xscale=log10,
    yscale=log10,
    xlabel="size of variant",
    ylabel="density of variants"
)
ylims!(0.00001, 0.008)
for tInd in eachindex(densV_f_T)
    scatterlines!(_f, densV_f_T[tInd], label="t = "*string(round(_T[tInd],digits=2)))
end
axislegend()
display(fig)