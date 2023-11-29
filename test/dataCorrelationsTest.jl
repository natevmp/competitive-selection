using CSV, DataFrames, Statistics, StatsBase
include("../src/dataStructuring.jl")
using .DataStructuring
using ProgressMeter
using CairoMakie

using JLD2


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

## ----------------------------------------
#region - Correlations

_pid = dropdupes(dfVid[!,:pid])

dfPatient = dfVid[dfVid.pid .== _pid[1], :]
dfPatient.dxdtx_t = [Float64[] for _ in eachrow(dfPatient)]

for dfVidCur in eachrow(dfPatient)
    dfVidCur.dxdtx_t = 2*(dfVidCur.vaf_t[2:end] .- dfVidCur.vaf_t[1:end-1]) ./ dfVidCur._t[1:end-1] ./ ( 2 .* dfVidCur.vaf_t[1:end-1].*(1 .- 2 .* dfVidCur.vaf_t[1:end-1]) )
end

# make pairs of dxdt1 - x2
dxdtx_tS_PairID = Vector{Float64}[]
x_tS_PairID = Vector{Float64}[]
_tS_PairID = Vector{Float64}[]
for (vid,dfVid1) in enumerate(eachrow(dfPatient))
    _t1 = dfVid1._t
    for dfVid2 in eachrow(dfPatient)[vid+1:end]
        _t2 = dfVid2._t
        _tS = [t for t in _t1 if t in _t2][1:end-1]
        if length(_tS)<3 continue end
        pick_t1 = [t in _t2 for t in _t1]
        pick_t2 = [t in _t1 for t in _t2]
        dxdtxA_tS = dfVid1.dxdtx_t[pick_t1[1:end-1]]
        xA_tS = dfVid2.vaf_t[pick_t2][1:end-1]
        push!(dxdtx_tS_PairID, dxdtxA_tS)
        push!(x_tS_PairID, xA_tS)
        push!(_tS_PairID, _tS)
        dxdtxB_tS = dfVid2.dxdtx_t[pick_t2[1:end-1]]
        xB_tS = dfVid1.vaf_t[pick_t1][1:end-1]
        push!(dxdtx_tS_PairID, dxdtxB_tS)
        push!(x_tS_PairID, xB_tS)
        push!(_tS_PairID, _tS)
    end
end

dxdtx_tS_PairID
x_tS_PairID
length.(x_tS_PairID)

cor_pairID = [corspearman(dxdtx_tS, x_tS) for (dxdtx_tS, x_tS) in zip(dxdtx_tS_PairID, x_tS_PairID)]

for cor in cor_pairID println(cor) end

for (dxdtx_tS, x_tS, _tS, cor) in zip(dxdtx_tS_PairID, x_tS_PairID, _tS_PairID, cor_pairID)
    println("variant: ")
    println("dxdt_t variant1:")
    println(dxdtx_tS)
    println("x_t variant2:")
    println(x_tS)
    println("")
    fig = Figure()
    Axis(fig[1,1], ylabel="x_1", title="correlation: "*string(round(cor, digits=2)))
    lines!(_tS, x_tS)
    Axis(fig[1,1], yaxisposition=:right, ylabel="dx_2/dt/x_2")
    lines!(_tS, dxdtx_tS, color=Cycled(2), linestyle=:dash)
    display(fig)
end

##
dfPatient.dxdtx_t
dfPatient.vaf_t

fig = Figure()
Axis(fig[1,1])
for dfVidCur in eachrow(dfPatient)
    lines!(dfVidCur._t, dfVidCur.vaf_t)
end
display(fig)
#endregion