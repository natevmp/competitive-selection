## ======================== Fabre Data analysis ======================

using CSV
using DataFrames
using Statistics
using LsqFit, StatsBase
include("../../src/competitiveSelection.jl")
using .CompetitiveSelection
using Optim, Distributions

##

df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
# remove 0 VAF "non-variant" calls (why do these exist?)
df = subset(df, :Gene => g -> g .!= "NA")

function fixArrTypeToType(arr, T::DataType)
    eltype(arr) == T && return arr
    newArr = Vector{T}(undef, length(arr))
    for (i,x) in enumerate(arr)
        try 
            newArr[i] = parse(T,x)
        catch
            newArr[i] = 0.
        end
    end
    return newArr
end
df[!,:VAF] .= fixArrTypeToType(df[!,:VAF], Float64)
df[!,:TOTALcount] .= fixArrTypeToType(df[!,:TOTALcount], Int64)
df[!,:WTcount] .= fixArrTypeToType(df[!,:WTcount], Int64)
df[!,:MUTcount_Xadj] .= fixArrTypeToType(df[!,:MUTcount_Xadj], Float64)
names(df)

##
function distbin(x_vid::Vector{T} where T<:Real; bins::Int=25, xMin::Real=0, xMax::Real=1., normalized::Bool=true)
    # println("any x<0: ", any(x_vid .< 0))
    # nVars = length(x_vid)
    _xEdges = range(xMin, xMax, length=bins+1)
    dx = _xEdges[2]-_xEdges[1]
    _f = _xEdges[1:end-1] .+ dx/2
    withinRangeBool_vid = (x -> x>=xMin && x<xMax).(x_vid)   #apply thresholds
    nV_f = zeros(Float64, bins)
    for x in x_vid[withinRangeBool_vid]
        k = ((x-xMin)/dx |> floor |> Int) + 1
        nV_f[k] += 1
    end
    # normalized && (nV_f ./= length(x_vid[withinRangeBool_vid])*dx)
    normalized && (nV_f ./= (length(x_vid)*dx))
    return _f, nV_f
end

function dropdupes(arr::Vector{<:Real})
    newArr = eltype(arr)[]
    for el in arr
        el in newArr ? continue : push!(newArr, el)
    end
    return newArr
end
function dropdupes(arr::AbstractVector)
    newArr = eltype(arr)[]
    for el in arr
        el in newArr ? continue : push!(newArr, el)
    end
    return newArr
end
_pid = dropdupes(df[!,:SardID])

nVars_pid = [
    subset(df, :SardID => id -> id .== pid, :Phase => p -> p .== 1, :VAF => f->f.>0)[!, :VAF] |> length
    for pid in _pid
]
_gene = dropdupes(df[!,:Gene])

##


Nf = 1E5
λ_vid = Float64[]
t0_vid = Float64[]
age_vid = Float64[]
λ_Gene_Pid = Vector{Vector{Float64}}(undef, length(_gene))
Age_Gene_Pid = Vector{Vector{Float64}}(undef, length(_gene))
vaf_Gene_Pid_Age = Vector{Vector{Vector{Float64}}}(undef, length(_gene))
_age_Gene_Pid = Vector{Vector{Vector{Float64}}}(undef, length(_gene))
@time for (gInd,gene) in enumerate(_gene)
    λ_Gene_Pid[gInd] = Float64[]
    Age_Gene_Pid[gInd] = Float64[]
    vaf_Gene_Pid_Age[gInd] = []
    _age_Gene_Pid[gInd] = []
    _gPid = subset(df,:Gene => g -> g.== gene,)[!,:SardID] |> dropdupes
    for pid in _gPid
        _tG = subset(df,:Gene => g -> g.== gene,:SardID => id -> id.== pid,)[!,:Age]
        if length(_tG)>=3
            vaf_tG = subset(df,:Gene => g -> g.== gene,:SardID => id -> id.== pid,)[!,:VAF]
            fit = CompetitiveSelection.fitLogGrowth(_tG, vaf_tG, Nf)
            append!(t0_vid, fit.param[1])
            append!(λ_Gene_Pid[gInd], fit.param[2])
            append!(Age_Gene_Pid[gInd], _tG[1])
            push!(vaf_Gene_Pid_Age[gInd], vaf_tG)
            push!(_age_Gene_Pid[gInd], _tG)
        end
    end
end
age_vid = append!(Float64[],[Age_Gene_Pid...]...)
λ_vid = append!(Float64[],[λ_Gene_Pid...]...)
ρS = corspearman(age_vid, λ_vid)

function allVariantsVafTime(df::DataFrame)
    vaf_Age_Vid = Vector{Float64}[]
    _Age_Vid = Vector{Float64}[]
    for (j,pid) in enumerate(_pid)
        _vid = 
        dropdupes(
            subset(df, :SardID => s->s.==pid)[!,"AAChange.refGene"]
        )
        for (i,vid) in enumerate(_vid)
            _age = subset(
                df,
                "AAChange.refGene" => ref->ref.==vid,
                "SardID" => s->s.==pid,
            )[!,:Age]
            vaf_age = subset(
                df,
                "AAChange.refGene" => ref->ref.==vid,
                "SardID" => s->s.==pid,
            )[!,:VAF]
            push!(vaf_Age_Vid, vaf_age)
            push!(_Age_Vid, _age)
        end
    end
    return _Age_Vid, vaf_Age_Vid
end

##
subset(
    df[!, [:SardID, :Phase, :Age, :Gene, :VAF, :START]],
    :SardID => id -> id .== 260
)
_vid = dropdupes(df[!,:START])
subset(df, :START => st->st.==_vid[1])[!,:Age]
subset(df, :START => st->st.==_vid[1])[!,:VAF]
##
bins = 27
function binVariantsTime(age_vid, λ_vid; bins::Int=21)
    _tBinEdge = range(minimum(age_vid),maximum(age_vid)+1, length=bins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    vaf_TBin_Vid = [Float64[] for i in eachindex(_tBin)]
    for (i,t) in enumerate(age_vid)
        iBin = 1 + floor((t-_tBinEdge[1])/Float64(_tBinEdge.step)) |> Int
        append!(vaf_TBin_Vid[iBin], λ_vid[i])
    end
    return _tBinEdge, vaf_TBin_Vid
end

function meanVarSizeTime(age_vid, λ_vid; bins::Int=21)
    _tBinEdge, vaf_TBin_Vid = binVariantsTime(age_vid, λ_vid; bins=bins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    vafAv_tBin = Vector{Float64}(undef, length(_tBin))
    for i in eachindex(vaf_TBin_Vid)
        vafAv_tBin[i] = mean(vaf_TBin_Vid[i])
    end
    return _tBin, vafAv_tBin
end

##
include("../../src/competitiveSelection.jl")
using .CompetitiveSelection

function analyseVariants(df::DataFrame, Nf::Real)
    stdVaf(n,p) = sqrt(n*p*(1-p))/n
    # _age_Vid, vaf_age_Vid = allVariantsVafTime(df)
    _pid = dropdupes(df[!,:SardID])
    dfVid = DataFrame(
        _t=Vector{Float64}[],
        vaf_t=Vector{Float64}[],
        pid=Int64[],
        gene=String7[],
        cov_t=Vector{Int64}[], # coverage
        n_t=Vector{Float64}[], # mutation counts
    )
    for (j,pid) in enumerate(_pid)
        _vid = 
        dropdupes(
            subset(df, :SardID => s->s.==pid)[!,"AAChange.refGene"]
        )
        for (i,vid) in enumerate(_vid)
            subDf = subset(
                df,
                "AAChange.refGene" => ref->ref.==vid,
                "SardID" => s->s.==pid,
            )
            length(subDf[:,:Age]) < 3 && continue # skip variant if less than 3 datapoints
            length( findall((v->v>0).(subDf[:,:VAF]) )) < 2 && continue
            _t = subDf[:,:Age]
            vaf_t = subDf[:,:VAF]
            gene = subDf[1,:Gene]
            cov_t = subDf[:, :TOTALcount] # coverage
            n_t = subDf[:, :MUTcount_Xadj] # mutation counts

            push!(dfVid, (_t, vaf_t, pid, gene, cov_t, n_t))
        end
    end
    nVars = size(dfVid,1)
    dfVid.t0 = Vector{Float64}(undef, nVars)
    dfVid.λ = Vector{Float64}(undef, nVars)
    dfVid.ρ = Vector{Float64}(undef, nVars)
    dfVid.x0 = Vector{Float64}(undef, nVars)
    dfVid._tt = Vector{
        StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    }(undef, nVars)
    dfVid.vafTh_tt = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafTh_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.stdVafTh_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.fitConverge = Vector{Bool}(undef, nVars)
    for i in 1:nVars
        # dfVid[i, :ρ] = corspearman(dfVid[i,:_t], dfVid[i,:vaf_t])
        # fit = CompetitiveSelection.fitLogGrowth(dfVid[i,:_t], dfVid[i,:vaf_t], Nf)
        # dfVid[i, :t0] = fit.param[1]
        # dfVid[i, :λ] = fit.param[2]
        # dfVid[i, :x0] = (dfVid[i, :ρ]>=0 ? 1/Nf : fit.param[3])
        fit = CompetitiveSelection.fitLogGrowthMLE(dfVid[i,:_t], dfVid[i,:vaf_t], dfVid[i,:cov_t], Nf)
        println("finished variant ", i)
        dfVid[i, :fitConverge] = Optim.converged(fit)
        _β = Optim.minimizer(fit)
        dfVid[i, :ρ] = corspearman(dfVid[i,:_t], dfVid[i,:vaf_t])
        dfVid[i, :t0] = _β[1]
        dfVid[i, :λ] = _β[2]
        dfVid[i, :x0] = (dfVid[i, :ρ]>=0 ? 1/Nf : _β[3])
        dfVid[i, :_tt] = range(dfVid[i,:_t][1],dfVid[i,:_t][end], length=25)
        dfVid[i, :vafTh_tt] = [
            CompetitiveSelection.fLogistic(t, dfVid[i,:t0], dfVid[i,:λ], dfVid[i,:x0], 0.5)
            for t in dfVid[i,:_tt]
        ]
        dfVid[i, :vafTh_t] = [
            CompetitiveSelection.fLogistic(t, dfVid[i,:t0], dfVid[i,:λ], dfVid[i,:x0], 0.5)
            for t in dfVid[i,:_t]
        ]
        dfVid[i, :stdVafTh_t] = Vector{Float64}(undef, length(dfVid[i,:_t]))
        for (j,t) in enumerate(dfVid[i,:_t])
            vafTh = CompetitiveSelection.fLogistic(t, dfVid[i,:t0], dfVid[i,:λ], dfVid[i,:x0], 1/2)
            dfVid[i,:stdVafTh_t][j] = stdVaf(dfVid[i,:cov_t][j], vafTh)
        end
    end
    return dfVid
end

##
using CairoMakie
##
paramsPop = Dict(
    :N => Int(2E5),
    :α => 1,
    # :s => 0.12,
    :s => 0.15,
    :σ => 0.08,
    :T => 105,
    :μ => 0.2,
    # :μ => 0,
    :sType => "gamma",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
@time dfVid = analyseVariants(df, paramsPop[:N])

##


_noConverge = findall((b->!b).([dfCur[:fitConverge] for dfCur in eachrow(dfVid)]))

##
varId=482
dfVid[varId,:][:cov_t]
dfVid[varId,:][:n_t]
dfVid[varId,:][:n_t] ./ dfVid[varId,:][:cov_t]

dfVid[varId,:][:vaf_t] .* dfVid[varId,:][:cov_t]
##


# dfVid[varId,:_t]
dfVid[varId,:vaf_t]
lines(dfVid[varId,:vaf_t])

dfVid[varId,:cov_t]

##
include("../../src/competitiveSelection.jl")
using .CompetitiveSelection

# for varId in 1:50
varId=31
@time res = CompetitiveSelection.fitLogGrowthMLE(dfVid[varId,:_t], dfVid[varId,:vaf_t], dfVid[varId,:cov_t], paramsPop[:N])
_β = Optim.minimizer(res)
# CompetitiveSelection.fitLogGrowthMLE(dfVid[i,:_t], dfVid[i,:vaf_t], dfVid[i,:cov_t], Nf)

vafFit_t = let 
    if _β[2]>0
        (t->CompetitiveSelection.fLogistic(t, _β[1], _β[2], 1/paramsPop[:N], 1/2)).(dfVid[varId,:_t])
    else
        (t->CompetitiveSelection.fLogistic(t, _β[1], _β[2], _β[3], 1/2)).(dfVid[varId,:_t])
    end
end
quantError_lu_T = [
    quantile.(Binomial(dfVid[varId,:cov_t][i], vafFit_t[i]), [0.01,0.99]) / dfVid[varId,:cov_t][i] 
    for i in eachindex(dfVid[varId,:_t])
]
#
fig0a = Figure()
Axis(fig0a[1,1])
scatter!(dfVid[varId,:_t], dfVid[varId,:vaf_t])
rangebars!(dfVid[varId,:_t], [qe_lu[1] for qe_lu in quantError_lu_T], [qe_lu[2] for qe_lu in quantError_lu_T],
    whiskerwidth = 20,
    cycle = [:color],
)
lines!(dfVid[varId,:_t], vafFit_t)
display(fig0a)
# end
##

## ======== plot per patient ===========
rscale=0.7
res = (rscale*700,rscale*500)
fsize=16

for pid in _pid
    dfPid = subset(dfVid, :pid => id->id.==pid)
    # size(dfCur,1) < 2 && continue
    fig = Figure(resolution=res, fontsize=fsize)
    ax = Axis(
        fig[1,1],
        title="patient "*string(pid),
    )
    # for i in 1:size(dfPid,1)
    for dfPidVid in eachrow(dfPid)
        scatter!(dfPidVid[:_t], dfPidVid[:vaf_t],
            linestyle=:solid,
            label=dfPidVid[:gene],
        )
        lines!(dfPidVid[:_tt], dfPidVid[:vafTh_tt],
            linestyle=:dash,
            label="s = "*string(round(dfPidVid[:λ],digits=2)),
        )
        quantError_lu_T = [
            quantile.(Binomial(dfPidVid[:cov_t][i], dfPidVid[:vafTh_t][i]), [0.01,0.99]) / dfPidVid[:cov_t][i]
            for i in eachindex(dfPidVid[:_t])
        ]
        rangebars!(
            dfPidVid[:_t], [qe_lu[1] for qe_lu in quantError_lu_T], [qe_lu[2] for qe_lu in quantError_lu_T],
            whiskerwidth = 20,
            cycle = [:color],
        )
    end
    Legend(fig[1,2], ax, framevisible = false)
    # figname = "dataFit_pid$(pid).png"
    # save("Figures/DataFits/PerPatient/"*figname, fig)
    # display(fig)
end


## ----------- Plot per variant ----------
rscale=0.7
res = (rscale*600,rscale*500)
fsize=16
# for vid in 1:size(dfVid,1)
# for (i,dfVCur) in enumerate(eachrow(dfVid[_noConverge, :]))
@showprogress for (i,dfVCur) in enumerate(eachrow(dfVid))
    fig = Figure(resolution=res, fontsize=fsize)
    ax = Axis(
        fig[1,1],
        title="patient "*string(dfVCur[:pid])*"; gene "*dfVCur[:gene],
    )
    scatter!(dfVCur[:_t], dfVCur[:vaf_t],
        linestyle=:solid,
        label=dfVCur[:gene],
    )
    lines!(dfVCur[:_tt], dfVCur[:vafTh_tt],
        linestyle=:dash,
        label="s = "*string(round(dfVCur[:λ],digits=2)),
    )
    quantError_lu_T = [
        quantile.(Binomial(dfVCur[:cov_t][i], dfVCur[:vafTh_t][i]), [0.01,0.99]) / dfVCur[:cov_t][i]
        for i in eachindex(dfVCur[:_t])
    ]
    rangebars!(
        dfVCur[:_t], [qe_lu[1] for qe_lu in quantError_lu_T], [qe_lu[2] for qe_lu in quantError_lu_T],
        whiskerwidth = 20,
        cycle = [:color],
    )
    legendpos = ( dfVCur[:ρ]>0 ? :lt : :rt )
    axislegend(position=legendpos)
    figname = "dataFit_vid$i.png"
    save("Figures/DataFits/PerVariant/"*figname, fig)
    # display(fig)
end


## ----------- Data Figures -------------

sum(length.(_age_Vid).>=4)
##
_age_Vid, vaf_age_Vid = allVariantsVafTime(df)
vids = [1,2,3,4,5]
fig0 = Figure()
Axis(
    fig0[1,1],
    xlabel="time",
    ylabel="variant allele frequency",
)
for i in eachindex(_age_Vid)
    lines!(_age_Vid[i], vaf_age_Vid[i])
end
display(fig0)

##


fig1 = Figure()
Axis(
    fig1[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
    title="Spearman correlation: ρ="*string(round(ρS,digits=3))
)

for i in eachindex(_gene[1:end.!=30])
    scatter!(
        Age_Gene_Pid[i],λ_Gene_Pid[i],
        label=_gene[i],
    )
end
# end
# axislegend()
# println(_gene[ind])
display(fig1)

_tBin, vafAv_tBin = meanVarSizeTime(age_vid, λ_vid; bins)
println(corspearman(_tBin[.!isnan.(vafAv_tBin)], vafAv_tBin[.!isnan.(vafAv_tBin)]))
fig2 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(fig2[1,1],
    xlabel="age at first measurement", 
    ylabel="mean fitness",
)
scatter!(_tBin, vafAv_tBin)
display(fig2)

nBins = 7:35
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, vafAv_tBin = meanVarSizeTime(age_vid, λ_vid; bins=bins)
    corS_nBins[i] = corspearman(_tBin[.!isnan.(vafAv_tBin)], vafAv_tBin[.!isnan.(vafAv_tBin)])
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

## -------- retest scatter of all fitnesses ----------------

fig4 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(
    fig4[1,1],
    xlabel="age at measurement",
    ylabel="fitted logistc growth rate",
    # title="Spearman correlation: ρ="*string(round(ρS,digits=3))
)
scatter!(
    [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)],
    dfVid[!, :λ],
    color=:black,
    markersize=5,
)
ylims!(-1,1)
display(fig4)
corspearman([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[dfVid[!,:λ].<1,:])], dfVid[dfVid[!,:λ].>1, :λ])
##

function binVariantsTime(age_vid, λ_vid; bins::Int=21)
    _tBinEdge = range(minimum(age_vid),maximum(age_vid)+1, length=bins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    vaf_TBin_Vid = [Float64[] for i in eachindex(_tBin)]
    for (i,t) in enumerate(age_vid)
        iBin = 1 + floor((t-_tBinEdge[1])/Float64(_tBinEdge.step)) |> Int
        append!(vaf_TBin_Vid[iBin], λ_vid[i])
    end
    return _tBinEdge, vaf_TBin_Vid
end

function meanVarSizeTime(age_vid, λ_vid; bins::Int=21)
    _tBinEdge, vaf_TBin_Vid = binVariantsTime(age_vid, λ_vid; bins=bins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    vafAv_tBin = Vector{Float64}(undef, length(_tBin))
    for i in eachindex(vaf_TBin_Vid)
        vafAv_tBin[i] = mean(vaf_TBin_Vid[i])
    end
    return _tBin, vafAv_tBin
end

##
bins = 25
_tBin, vafAv_tBin = meanVarSizeTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)], dfVid[!, :λ]; bins)
println(corspearman(_tBin[.!isnan.(vafAv_tBin)], vafAv_tBin[.!isnan.(vafAv_tBin)]))
fig1 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(fig1[1,1],
    xlabel="age at first measurement", 
    ylabel="mean fitness",
)
scatter!(_tBin, vafAv_tBin)
display(fig1)

##
nBins = 6:30
corS_nBins = Vector{Float64}(undef, length(nBins))
for (i,bins) in enumerate(nBins)
    _tBin, vafAv_tBin = meanVarSizeTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)], dfVid[!, :λ]; bins)
    corS_nBins[i] = corspearman(_tBin[.!isnan.(vafAv_tBin)], vafAv_tBin[.!isnan.(vafAv_tBin)])
end
fig = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(
    fig[1,1],
    xlabel="number of bins",
    ylabel="spearman correlation of mean fitness with time"
)
scatter!(nBins, corS_nBins)
display(fig)

















## ===================== Simulations ==========================
include("../../src/competitiveSelection.jl")
using .CompetitiveSelection
# using Statistics
using CairoMakie
# using ElasticArrays, StatsBase, LsqFit

##
myTheme = Theme(fontsize = 20)
set_theme!(myTheme)
##
paramsPop = Dict(
    :N => Int(2E4),
    :α => 1,
    # :s => 0.12,
    :s => 0.15,
    :σ => 0.08,
    :T => 105,
    :μ => 0.2,
    # :μ => 0,
    :sType => "gamma",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 10000
##
@time solEns = CompetitiveSelection.evolvePop(paramsPop, runs=runs)

##
measureTime=50
# @time λ_tVaf_cid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop; t=measureTime)
@time λ_tVaf_cid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop)
##
fig = scatter(λ_tVaf_cid[1,:], λ_tVaf_cid[2,:])
xlims!(50,100)
ylims!(-1,5)
display(fig)
##
# thresh = 2
# fig4 = Figure()
# Axis(
#     fig4[1,1],
#     xlabel="measured fitness",
#     ylable="density of variants",
#     title="t=$(measureTime); s=$(paramsPop[:s]); σ=$(paramsPop[:σ]); q=$(paramsPop[:q])",
# )
# hist!(λ_tVaf_cid[2, λ_tVaf_cid[2,:].<=thresh]; bins=50, normalization=:pdf)
# xlims!(-1,thresh)
# display(fig4)

##
_tSimBin, vafSimAv_tSimBin = meanVarSizeTime(λ_tVaf_cid[1,:], λ_tVaf_cid[2,:]; bins)

println(corspearman(_tSimBin[.!isnan.(vafSimAv_tSimBin)], vafSimAv_tSimBin[.!isnan.(vafSimAv_tSimBin)]))

fig5 = Figure(
    resolution=(700,600),
    fontsize=22,
)
Axis(fig5[1,1],
    xlabel="age at first measurement", 
    ylabel="mean fitness",
)
scatter!(_tSimBin, vafSimAv_tSimBin)
display(fig5)
