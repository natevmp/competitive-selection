include("../src/dataStructuring.jl")
using .DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools
include("../src/abcSims.jl")
using .ABCSims
using JLD2, CSV, DataFrames
using CairoMakie

##
# include(ARGS[1])
# include("../HPC/params/paramsChisquareSFixed6571.jl")
tBounds = (70,76)

modelParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => round((tBounds[1]+tBounds[2])/2) |> Int,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
)

ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 500,
    :fBins => 25,
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
)


##

println("Available threads: ", Threads.nthreads())

##
# ----- Data --------
df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
DataStructuring.structureData!(df)
dfVid = DataStructuring.sortDataPerVariant(df)
x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
xMask_sid = x_sid.>0
@time _fEdges, nData_f = AnalysisTools.sizeDistDens(
    x_sid[xMask_sid], t_sid[xMask_sid], tBounds; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax]
)
nCumData_f = [sum(nData_f[1:i]) for i in eachindex(nData_f)]
data_f = ctrlParams[:cumulativeDist] ? nCumData_f : nData_f

## ----- plot sizeDist data -----
_f = (_fEdges[2]-_fEdges[1])/2 .+ _fEdges[1:end-1]
fig = Figure()
Axis(fig[1,1], xscale=log10, yscale=log10)
ylims!(1E-5, 2E-2)
# ylims!()
lines!(_f, data_f)
display(fig)

## ----- ABC --------
@time abcRes = ABCSims.runABC(data_f, ctrlParams)

##
propertynames(abcRes)
abcRes.population
abcRes.n_tries
abcRes.n_accepted


## --------- save data -------------
println("saving data")
filename = begin
    "abcResult_"*string(ctrlParams[:metric])*
    "_ages"*string(tBounds[1])*string(tBounds[2])*
    # "_cumulative"*string(ctrlParams[:cumulativeDist])*
    (ctrlParams[:cumulativeDist] ? "_cumulative" : "")*
    (!haskey(ctrlParams, :sBounds) ? "_SFixed" : "")*
    (!haskey(ctrlParams, :σBounds) ? "_SFixed" : "")*
    "_particles"*string(ctrlParams[:nParticles])*
    ".jld2"
end
    # jldsave(filename; abcRes, ctrlParams)

