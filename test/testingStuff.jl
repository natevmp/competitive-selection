using CSV, DataFrames, Statistics, StatsBase
# using Optim, Distributions
# include("../../src/competitiveSelection.jl")
# using .CompetitiveSelection
include("../src/dataStructuring.jl")
using ProgressMeter
using .DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools
using Optim

##
df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
Nf = 2E5
##
# dfVid = DataStructuring.analyseVariants(df, Nf)
@time dfVid = DataStructuring.analyseVariants(df, Nf; competition=false)


##
for vid in 1:10
    fig = Figure(resolution=(600,600), fontsize=30)
    Axis(fig[1,1],
        # title="s=$(round(dfVid[vid, :s],digits=2)), a=$(round(dfVid[vid, :a],digits=2)), μ=$(round(dfVid[vid, :μ],digits=2))"
    )
    lines!(dfVid[vid,:_tt], dfVid[vid,:vafTh_tt])
    scatter!(dfVid[vid,:_t], dfVid[vid, :vaf_t])
    display(fig)
end
##

fig = Figure(resolution=(800,600), fontsize=24)
Axis(fig[1,1],xlabel="time", ylabel="inferred fitness")
scatter!([dfCur[:_t][1] for dfCur in eachrow(dfVid)], dfVid[!,:s],label="individual fitness")
scatter!([dfCur[:_t][1] for dfCur in eachrow(dfVid)], dfVid[!,:a],label="average fitness")
axislegend(position=:rt)
display(fig)
##
fig = Figure(resolution=(800,600), fontsize=24)
Axis(fig[1,1], xlabel="time", ylabel="inferred mutation rate")
scatter!([dfCur[:_t][1] for dfCur in eachrow(dfVid)], dfVid[!,:μ],label="μ")
axislegend()
display(fig)

##
println(dfVid[!,:t0])
##

@showprogress for dfVCur in eachrow(dfVid[!,:])
    DataStructuring.fitLogGrowthMLE(dfVCur[:_t], dfVCur[:vaf_t], dfVCur[:cov_t], Nf; constAdjustedLogit=true)
end

##
vid = 9
dfit = DataStructuring.fitLogGrowthMLE(dfVid[vid,:_t], dfVid[vid,:vaf_t], dfVid[vid,:cov_t], Nf; constAdjustedLogit=true)
ρ = corspearman(dfVid[vid,:_t], dfVid[vid,:vaf_t])
_β = Optim.minimizer(dfit)
t0 = _β[1]
λ = _β[2]
z = _β[3]
x0 = (ρ>=0 ? 1/Nf : _β[4])
_tt = range(dfVid[vid,:_t][1], dfVid[vid,:_t][end], length=25)
vafTh_tt = [
            AnalysisTools.fLogistic(t, t0, λ, x0, 0.5, z)
            for t in _tt
        ]

fig = Figure()
Axis(fig[1,1])
scatter!(dfVid[vid,:_t], dfVid[vid,:vaf_t])
lines!(_tt, vafTh_tt)
display(fig)
println("t0: ", t0)
println("x0: ", x0)
println("λ: ", λ)
println("z: ", z)

##
trackerVariant=3.2
dfTest = DataFrame()
dfTest.k = fill(10,10)
dfTest.t₀_vid = (k->(10*rand(k) |> sort!)).(dfTest[!,:k])
if !isnothing(trackerVariant)
    for row in eachrow(dfTest)
        row.k += 1
        push!(row.t₀_vid, trackerVariant)
        sort!(row.t₀_vid)
    end
end

##

dfVid = DataStructuring.analyseVariants(df, Nf)
##

using LinearAlgebra

function linearRegression(x::AbstractVector{T} where T<:Real, y::Vector{T} where T<:Real)
    ỹ = @. !isnan(y)
    n = sum(ỹ)
    X = [ones(n) x[ỹ]]
    β̂ = X \ y[ỹ]
    return β̂
end

##

bins = 15
_tBin, λAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)], dfVid[!, :λ]; bins)
fig3a = Figure()
Axis(
    fig3a[1,1],
    xlabel="age at measurement",
    ylabel="average growth rate",
    title="bins = "*string(bins)
)
scatter!(_tBin, λAv_tBin)
λ0, slope = linearRegression(_tBin, λAv_tBin)
relSlope = round(slope/λ0, digits=4)
lines!(_tBin, (t-> λ0+slope*t).(_tBin),label="relative slope = "*string(relSlope))
axislegend()
display(fig3a)
# saveFolder = "../../Figures/Presentations/"
# save(saveFolder*"binnedAvFitnessLogit.png", fig3a, px_per_unit = 2)


##

_bins = 6:30
λ0_bins = Vector{Float64}(undef, length(_bins))
λslope_bins = Vector{Float64}(undef, length(_bins))
for (i, bins) in enumerate(_bins)
    _tBin, λAv_tBin = DataStructuring.meanVarFitParamTime([dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid)], dfVid[!, :λ]; bins)
    λ0, slope = linearRegression(_tBin, λAv_tBin)
    λ0_bins[i] = λ0
    λslope_bins[i] = slope / λ0
end

##

fig4 = Figure()
Axis(fig4[1,1],
    xlabel="number of bins",
    ylabel="relative linear slop (lsq regression)"
)
scatter!(_bins, λslope_bins)
display(fig4)

##

symbolList = [:a, :b]

dfTest = DataFrame()

dfTest.:aa = [1,2]
dfTest[:, symbolList[1]] = [2.,3.]
dfTest

##

struct TestType
    a
    b
end

test = TestType(4, 5)

fieldnames(TestType)
##
testD = Dict(
    :a => 1,
    :b => 2.
)
keyL = keys(testD)
df = DataFrame()
for key in keyL
    df[:,key] = Vector{Float64}(undef, 4)
end
df

##
using DataFrames
using Optimization, OptimizationOptimJL, Distributions
include("../src/analysisTools.jl")
using .AnalysisTools
##

abstract type ModelShape end
mutable struct LogisticModel <: ModelShape 
    γ::Float64
    t0::Float64
    x0::Float64
    xF::Float64
    LogisticModel() = new(0.1, 1.,0.1, 0.5)
end
mutable struct MaxAdjustedLogisticModel <: ModelShape
    γ::Float64
    t0::Float64
    x0::Float64
    xF::Float64
    MaxAdjustedLogisticModel() = new(0.1, 1., 0.1, 0.25)
end

abstract type ExpansionType end
struct PositiveExpansion <: ExpansionType end
struct NegativeExpansion <: ExpansionType end

struct ModelFit{S<:ModelShape, E<:ExpansionType}
    shape::S
    expansion::E
end
ModelFit(shape, exp) = ModelFit{typeof(shape), typeof(exp)}(shape, exp)
##

shapeA = MaxAdjustedLogisticModel()
expA = PositiveExpansion()
modelA = ModelFit(shapeA, expA)
display(modelA)

##

# b0[1]=t0, b0[2]=γ, b0[3]=xF
function setInitialParamsMLE!(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 1.
    model.shape.γ = 0.05
    model.shape.x0 = 1/(2Nf)
    model.shape.xF = 0.25
    return [model.shape.t0, model.shape.γ, model.shape.xF]
end
function boundsMLE(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _t, vaf_t)
    ([-1E6, 0, 1E-4],
    [_t[1], 10, 0.49])
end
# b0[1]=t0, b0[2]=γ, b0[3]=x0, b0[4] = xF
function setInitialParamsMLE!(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _t, vaf_t, Nf)
    model.shape.t0 = 1.
    model.shape.γ = -0.1
    model.shape.x0 = vaf_t[1]
    model.shape.xF = 0.45
    return [model.shape.t0, model.shape.γ, model.shape.x0, model.shape.xF]
end
function boundsMLE(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _t, vaf_t)
    ([-1E6, -Inf, 0, 1E-4],
    [_t[1], 0, 0.45, 0.49])
end

function fLogModel(t,β,model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion})
    fLogistic(t, β[1], β[2], model.shape.x0, β[3])
end

function fLogModel(t,β,model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion})
    fLogistic(t, β[1], β[2], β[3], β[4])
end 

function addFitParams!(model::ModelFit{MaxAdjustedLogisticModel, PositiveExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
    model.shape.xF = _β[3]
end
function addFitParams!(model::ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}, _β)
    model.shape.t0 = _β[1]
    model.shape.γ = _β[2]
    model.shape.x0 = _β[3]
    model.shape.xF = _β[4]
end

##
function fitLogGrowthMLE!(_t, vaf_t, nCov_t, Nf, timeLimit=1.0)
    corr = corspearman(_t, vaf_t)
    # println(corr)
    if corr < 0
        return fitLogGrowthMLENeg!(_t, vaf_t, nCov_t, Nf; timeLimit)
    else
        return fitLogGrowthMLEPos!(_t, vaf_t, nCov_t, Nf; timeLimit)
    end
end

function fitLogGrowthMLENeg!(_t, vaf_t, nCov_t, Nf; timeLimit=1.0)
    X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on
    # println("negative growth")
    # 1: t0
    # 2: ρ
    # 3: x0
    # 4: xF
    function logprobX(X, t, n, β)
        # println("\nβ: ", β)
        p = AnalysisTools.fLogistic(t, β[1], β[2], β[3], β[4])
        #! debug
        if p < 0 || p > 1 || isnan(p)
            # println("\nβ: ", β)
            # println("p: ", p)
            # println("t: ", t)
            p = 1E-30
        end
        #! debug
        logpdf(Binomial(n, p), X)
    end
    logllh(β, p) = -([
        logprobX(X_t[i], _t[i], nCov_t[i], β) for i in eachindex(_t)
    ] |> sum)
    x0 = [20, -0.5, 0.2, 0.3]
    cons(res, x, p) = (res .= [x[1], x[2], x[3], x[4], x[3] - x[4]])
    lcons = [0, -10, 0, 0, -0.5]
    ucons = [_t[1], 0, 0.48, 0.49, -0.01]
    optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff(), cons = cons)
    # optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optprob, x0, lcons=lcons, ucons=ucons)
    # sol = solve(prob, IPNewton())
    sol = solve(prob, Ipopt.Optimizer(); print_level = 0)
    
    # println("[t0, ρ, x0, xF]")
    return sol
end

function fitLogGrowthMLEPos!(_t, vaf_t, nCov_t, Nf; timeLimit=1.0)
    X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on
    # println("positive growth")
    function logprobX(X, t, n, β)
        # println("\nβ: ", β)
        p = AnalysisTools.fLogistic(t, β[1], β[2], 1/(2Nf), β[3])
        #! debug
        if p < 0 || p > 1 || isnan(p)
            # println("\nβ: ", β)
            # println("p: ", p)
            # println("t: ", t)
            # p = 1E-50
        end
        #! debug
        logpdf(Binomial(n, p), X)
    end
    logllh(β, p) = -([
        logprobX(X_t[i], _t[i], nCov_t[i], β) for i in eachindex(_t)
    ] |> sum)
    x0 = [20, 0.5, 0.25]
    lcons = [0, 0, 0]
    ucons = [_t[1], 10, 0.49]
    optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optprob, x0, lb=lcons, ub=ucons)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
    # println("[t0, ρ, xF]")
    return sol
end
##
df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
Nf = 2E5
DataStructuring.structureData!(df)
dfVidTest = DataStructuring.sortDataPerVariant(df);

##
using Optimization, OptimizationBBO, OptimizationMOI, Ipopt

vid = 20
sol = fitLogGrowthMLE!(dfVidTest[vid,:_t], dfVidTest[vid,:vaf_t], dfVidTest[vid,:cov_t], Nf)
# @time [fitLogGrowthMLE!(dfVidTest[vid,:_t], dfVidTest[vid,:vaf_t], dfVidTest[vid,:cov_t], Nf) for vid in 1:200]

4+4


##



rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2

x0 = zeros(2)
_p = [1.0, 1.0]

cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])
# cons = (res, x, p) -> res .= [x[1]^2 + x[2]^2]

optprob = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(); cons = cons)
prob = OptimizationProblem(optprob, x0, _p, lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
sol = solve(prob, Ipopt.Optimizer())




##

tSingleSim = 1
nParams=3
20^nParams
20^nParams * tSingleSim / 3600