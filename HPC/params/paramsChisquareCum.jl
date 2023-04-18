
const tBounds = (70,80)
const fixedParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => (tBounds[1]+tBounds[2])/2 |> Int,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
)

const ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 1000,
    :fBins => 25,
    :nSamples => 3,
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :sBounds => (0.01, 1.5),
    :σBounds => (0.01, 0.5),
    :μBounds => (0.1, 6),
    :params => fixedParams,
    :cumulativeDist => true,
    :metric => "chisquare",
    :threshold => [0.001, 0.0005, 0.00001],
    :nParticles => 300,
    :maxIter => 10000,
)