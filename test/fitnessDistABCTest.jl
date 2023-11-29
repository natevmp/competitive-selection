include("../src/dataStructuring.jl")
using .DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools
include("../src/abcSims.jl")
using .ABCSims
include("../src/abcSimsBeta.jl")
using .ABCSimsBeta
using JLD2, CSV, DataFrames
using CairoMakie

##
# include(ARGS[1])
# include("../HPC/params/paramsChisquareSFixed6571.jl")
tBounds = (65,75)

modelParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => round((tBounds[1]+tBounds[2])/2) |> Int,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
)

ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 100,
    :fBins => 25,
    :sBins => 50,
    :sFitBounds => (-0.5, 1),
    :sQError => 0.95,
    :nSamples => 3,
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :coverage => 1000,
    :sBounds => (0.01, 1.),
    :σBounds => (0.005, 0.1),
    :μBounds => (0.1, 8),
    :params => modelParams,
    :cumulativeDist => false,
    :metric => "chisquare",
    # :threshold => [1E-2, 4E-3, 1E-3, 7E-4, 4E-4],
    :threshold => [4E-3,],
    :nParticles => 10,
    :maxIter => 20,
    :nTimeSteps => 4,
    :tStep => 2.5,
)

##



# @time n_s, nVarsThreshAv = ABCSimsBeta.modelSim([0.5, 0.04, 7], ctrlParams)
n_s = vec(n_s)
println(nVarsThreshAv)
##
using CairoMakie

_s = range(ctrlParams[:sFitBounds][1], ctrlParams[:sFitBounds][2], length=ctrlParams[:sBins])

barplot(_s, n_s)
