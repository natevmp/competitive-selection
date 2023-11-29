include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using CSV, DataFrames
# using Statistics
# include("../src/mullerPlot.jl")
# using .MullerPlot

## ========================================
#region = parameters + Run sims

params = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 1.,
    :μ => 3,
    :s => 0.09,
    :σ => 0.025,
)
ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 500,
    :fBins => 25,
    :sBins => 50,
    :sFitBounds => (-0.5, 1),
    :sQError => 0.95,
    :nSamples => 2, # number of variants to sample per patient (sim) when performing logistic fit
    :nTimeSteps => 4, # number of time points per variant to sample when performing logistic fit
    :tStep => 2.5, # time between sample time points when performing logistic fit 
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :coverage => 1000, # sequencing coverage to use in sampling
    :fixPar => Symbol[],
    :sBounds => (0.01, 0.4), # ABC particle parameter bounds
    :σBounds => (0.005, 0.1), # ABC particle parameter bounds
    :μBounds => (0.1, 15), # ABC particle parameter bounds
    :params => modelParams, # model parameters to use for simulations
    # -------- rejection parameters --------
    :nParticles => 1,
    :maxIter => 20,
    :metric => "chisquare",
    :cumulativeDist => false,
    :thresholds => (0.1,(0.1,0.1)),
    :nVarsBoundsLate => (1,25),
    :nVarsBoundsEarly => (0,2),
)

@time solEns, simArgs = CompetitiveSelection.evolvePopSim(
    params;
    runs=500,
    noDiffusion=false
)
sol = solEns[1]
#endregion

##

tMeasure = 72
parentVid_vid_Sim = simArgs[!, :parentId_vid]
vaf_t = CompetitiveSelection.sizeDistSims(tMeasure, solEns, parentVid_vid_Sim, ctrlParams)

_tEdges = range(ctrlParams[:fMin], ctrlParams[:fMax], length=ctrlParams[:fBins]+1)
_t = _tEdges[1:end-1] .+ (_tEdges[2]-_tEdges[1])/2
fig = Figure()
Axis(fig[1,1],xscale=log10, yscale=log10)
# xlims!(1E-1, 0.5)
ylims!(1E-5, 0.02)
scatterlines!(_t, vaf_t)
display(fig)

##
using Distributions

vec = rand(Uniform{Float64}(0.001, 0.1), 100000)
hist(vec, bins=100)

