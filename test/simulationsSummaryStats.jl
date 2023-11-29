include("../src/dataStructuring.jl")
using .DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools
include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using CSV, DataFrames
using Statistics
include("../src/abcSims.jl")
using .ABCSims


## ----------------- Load data -------------

df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
DataStructuring.structureData!(df)
dfVid = DataStructuring.sortDataPerVariant(df)

x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
xMask_sid = x_sid.>0

## ---------------- ABC params --------------
tBounds = (70,80)
fixedParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => (tBounds[1]+tBounds[2])/2 |> Int,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
)
ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 1000,
    :fBins => 25,
    :nSamples => 3,
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :sBounds => (0.01, 1.5),
    :σBounds => (0.01, 0.2),
    :μBounds => (0.1, 5),
    :params => fixedParams,
    :cumulativeDist => false,
    :metric => "mean",
    :threshold => [0.1, 0.02, 0.005],
    :nParticles => 2,
    :maxIter => 10,
)

##

_fEdges, nData_f = AnalysisTools.sizeDistDens(
    x_sid[xMask_sid], t_sid[xMask_sid], tBounds; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax]
)
nCumData_f = [sum(nData_f[1:i]) for i in eachindex(nData_f)]


# ## ------------ Test sims -----------

# params = Dict{Symbol,Union{Int64,Float64,String}}(
#     :T => 75,
#     :N => 200000,
#     :α => 1,
#     :sType => "gamma",
#     :q => 0,
#     :μ => 1,
#     :s => 0.1,
#     :σ => 0.05,
# )
# @time solEns, _ = CompetitiveSelection.evolvePopSim(params; runs=500, noDiffusion=false)
# println("")


## ============================ Testing ABC ======================

##

# @time nSim_f = vec(ABCSims.modelSim([0.114, 0.018, 2], ctrlParams))
@time nSim_f = vec(ABCSims.modelSim([0.114, 0.022, 3], ctrlParams))
nSimCum_f = [sum(nSim_f[1:i]) for i in eachindex(nSim_f)]

##
using CairoMakie
fig = Figure()
Axis(
    fig[1,1],
    xscale=log10, yscale=log10,
    xlabel="f",
    ylabel="density"
)
ylims!(1E-5, 0.61E-1)
# ylims!(1E-5, 0.61E-1)
lines!(nData_f, label="data")
# lines!(nCumData_f, label="data")
lines!(nSim_f, label="sims")
# lines!(nCumSim_f, label="sims", linestyle=:dash)
axislegend(position=:rt)
display(fig)

##
using Distances
euclidean(mean(nCumData_f), mean(nSim_f))
evaluate(Distances.chisq_dist, nCumData_f, nSim_f)
evaluate(Distances.chisq_dist, nData_f, nSim_f)
ABCSims.kolmogorovSmirnovDist(nCumData_f, nSimCum_f)
##

sum(nCumData_f)
sum(nCumSim_f)

mean(nData_f)
mean(nSim_f)



##
include("../src/abcSims.jl")
using .ABCSims

##

tBounds = (70,80)
fixedParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => (tBounds[1]+tBounds[2])/2 |> Int,
    :N => 200000,
    :α => 1,
    :s => 0.15,
    :sType => "gamma",
    :q => 0,
)
ctrlParams = Dict{}(
    :tBounds => tBounds,
    :simRuns => 1000,
    :fBins => 25,
    :nSamples => 3,
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :sBounds => (0.01, 1.5),
    :σBounds => (0.01, 0.2),
    :μBounds => (0.1, 5),
    :params => fixedParams,
    :cumulativeDist => false,
    :metric => "chisquare",
    :threshold => [0.5, 0.1, 0.1],
    :nParticles => 2,
    :maxIter => 4,
)

## ----- Data --------
tBounds = (70,80)
x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
xMask_sid = x_sid.>0
_fEdges, nData_f = AnalysisTools.sizeDistDens(
    x_sid[xMask_sid], t_sid[xMask_sid], tBounds; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax]
)
nCumData_f = [sum(nData_f[1:i]) for i in eachindex(nData_f)]

## ----- ABC --------
@time abcRes = ABCSims.runABC(nData_f, ctrlParams)


##

# @profview ABCSims.runABC(nCumData_f, ctrlParams)

# @code_warntype ABCSims.runABC(nCumData_f, ctrlParams)




##

abcRes.n_accepted
abcRes.n_tries

abcRes.population

fig1 = Figure(resolution=(1300,1000))
Axis(
    fig1[1,1],
    xlabel="s",
)
hist!(abcRes.population[:,1])
Axis(
    fig1[1,2],
    xlabel="σ",
)
hist!(abcRes.population[:,2])
Axis(
    fig1[2,2],
    xlabel="μ",
)
hist!(abcRes.population[:,3])
display(fig1)




##

fig2 = Figure(resolution=(1300,1000))
Axis(
    fig1[1,1],
    xlabel="s",
    ylabel="σ"
)
scatter!(abcRes.population[:,1], abcRes.population[:,2])
Axis(
    fig1[1,2],
    xlabel="s",
    ylabel="μ"
)
scatter!(abcRes.population[:,1], abcRes.population[:,3])
Axis(
    fig1[2,2],
    xlabel="σ",
    ylabel="μ"
)
scatter!(abcRes.population[:,2], abcRes.population[:,3])
display(fig1)











##

# runs = 1
# @time sol, simArgs = CompetitiveSelection.evolvePopSim(params; runs, getArgs=true, noDiffusion=false)

## ------- Individual trajectories --------
using CairoMakie
scale=0.9
makieTheme = Theme(
    fontsize=16,
    resolution=(scale*800,scale*600),
)
set_theme!(makieTheme)
##
sid = 4
sol = solEns[sid]
f_t_vid = [[sol[t][vid] for t in eachindex(sol)] for vid in eachindex(sol[1])]



fig = Figure()
Axis(fig[1,1],
    xlabel="time",
    ylabel="clone frequency",
)
for vid in eachindex(sol[1])
    lines!(f_t_vid[vid])
end
display(fig)

# figname = "individualTrajectories.png"
# save("Figures/Presentations/"*figname, fig)

##



##

# sid = 12
# f_t_Vid = [[solEns[sid][t][i] for t in eachindex(solEns[sid])] for i in eachindex(solEns[sid][1])]

# fig = Figure()
# Axis(fig[1,1])
# for f_t in f_t_Vid
#     lines!(f_t)
# end
# display(fig)

##

@time _f, nV_t_f = CompetitiveSelection.variantsSizeDistribution(solEns)

##

fig = Figure()
Axis(
    fig[1,1],
    xscale=log10,
    yscale=log10,
)
# xlims!(1E-5,1E3)
ylims!(1E-5,1)
for t in [50,60,70,80,90,100]
    lines!(_f, nV_t_f[t, :])
end
display(fig)
