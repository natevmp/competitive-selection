include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using Statistics
using CairoMakie

##
params = Dict(
    :N => 200000,
    :α => 1,
    :s => 0.15,
    :σ => 0.03,
    :T => 100,
    :μ => 2,
    :sType => "fixed",
    # :sType => "gamma",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
##
runs = 1000
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(params; runs, noDiffusion=false)
println("")

## ----------------- Total variant fraction above threshold -------------------
_t = 0:params[:T]
# sizeThreshold = 0.01
sizeThreshold = 0
vCovAv_t = Vector{Float64}(undef, length(_t))
for tInd in eachindex(_t)
    vCovAv_t[tInd] = mean([ sum(sol[tInd][sol[tInd].>sizeThreshold]) for sol in solEns ])
end


# for sol in solEns
#     println(sol[tInd])
# end

vCovTh_t = [Theory.compCoverage(t, params[:α], params[:s], params[:μ], params[:N]) for t in _t]

rscale=0.6
fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="total variant fraction"
)
lines!(_t, vCovTh_t, label="theory")
lines!(_t, vCovAv_t, label="simulations average", linestyle=:dash)
axislegend(position=:lt)
display(fig)
# saveFolder="Figures/Presentations/"
# save(saveFolder*"totalVariantFraction_ThVSims.png", fig, px_per_unit = 2)

## --------------- Size distribution ----------------

tInd = 91
sizeThreshold = 0.01

f_vid = Float64[]
for sol in solEns
    append!(f_vid,sol[tInd][sol[tInd].>sizeThreshold])
end
bins = 20
_f, n_f = CompetitiveSelection.distbin(f_vid; xMin=0, xMax=0.5, bins=bins)

fig = Figure()
Axis(fig[1,1],
    xscale=log10,
    yscale=log10,
)
xlims!(0.01, 0.5)
ylims!(1E-5,1E-1)
lines!(_f, n_f)
display(fig)

sizeDistDF = DataFrame()
sizeDistDF._f = _f
sizeDistDF.n_f = n_f

filename = "variantSizeDist"*
    "_t"*string(solEns[1].t[tInd])*
    "_bins"*string(bins)*``
    ".csv"

sizeDistDF |> CSV.write(filename)
## ================== Export result to CSV file ===================

using CSV, DataFrames

filename = "variantFractionTime"*
    "_s"*string(params[:s])*
    "_sigma"*string(params[:σ])*
    "_mu"*string(params[:μ])*
    ".csv"

# vCovMat_t = Array{Float64,2}(undef, (length(_t),3))
vCovMat_t = DataFrame()
vCovMat_t[:,:_t] = _t
vCovMat_t[:,:theory] = vCovTh_t
vCovMat_t[:,:simsAverage] = vCovAv_t
vCovMat_t |> CSV.write(filename)

## ==================================================

fig = Figure()
Axis(fig[1,1])
lines!(
    _t, (t->Theory.fCompetition(t, 0, params[:α], params[:s], params[:μ], params[:N])).(_t)
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 0, params[:α], params[:s]-.01, params[:s], params[:μ], params[:N]))
)
# t->fCompetition(t, 0, params[:α], params[:s], params[:μ], params[:N])
# t->fIndFitnessWithAvCompetition(t, 10, params[:α], params[:s]-.02, params[:s], params[:μ], params[:N])
display(fig)

##

params = Dict(
    :N => 200000,
    :α => 1,
    :s => 0.3,
    # :σ => 0.25,
    :T => 100,
    :μ => 1,
    :sType => "fixed",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)

## -------------Theory plots ---------------

fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="variant size",
    title="fixed fitness",
)
lines!(
    _t, (t->Theory.compCoverage(t, params[:α], params[:s], params[:μ], params[:N])).(_t),
    label="total variant fraction"
)
lines!(
    _t, (t->Theory.fCompetition(t, 1, params[:α], params[:s], params[:μ], params[:N])).(_t),
    label="t₀ = 1"
)
lines!(
    _t, (t->Theory.fCompetition(t, 5, params[:α], params[:s], params[:μ], params[:N])).(_t),
    label="t₀ = 5"
)
lines!(
    _t, (t->Theory.fCompetition(t, 10, params[:α], params[:s], params[:μ], params[:N])).(_t),
    label="t₀ = 10"
)
axislegend(position=:lt)
display(fig)
# save(saveFolder*"theoryFixedFitness.png", fig, px_per_unit=2)

##

trackerFitness = 0.19
fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="variant size",
    title="tracker variant in average background a = "*string(params[:s]),
)

lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 5, params[:α], 0.29, params[:s], params[:μ], params[:N])).(_t),
    label="s < a"
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 5, params[:α], 0.3, params[:s], params[:μ], params[:N])).(_t),
    label="s = a"
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 5, params[:α], 0.31, params[:s], params[:μ], params[:N])).(_t),
    label="s > a"
)

axislegend(position=:lt)
display(fig)
# save(saveFolder*"theoryTrackerFitness.png", fig, px_per_unit=1)

##
runs = 3000
trackerBirthTime = 10
# trackerFitness = 0.33
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(params; runs,trackerVariant=trackerBirthTime, getArgs=true, noDiffusion=false)
sizeTrackedAverage_t = CompetitiveSelection.averageTrackedVariant(solEns, simArgs)
##
fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="variant size",
    title="t₀ = "*string(trackerBirthTime)
)
lines!(
    _t, (t->Theory.fCompetition(t, trackerBirthTime, params[:α], params[:s], params[:μ], params[:N])).(_t),
    label="theory"
)
lines!(_t, sizeTrackedAverage_t, label="sims average")
axislegend(position=:lt)
display(fig)
# save(saveFolder*"trackedSimVTheory_u"*string(trackerBirthTime)*".png", fig, px_per_unit = 1)

##
##
ensID = 47
trackId = simArgs[ensID, :trackerID]

any([solEns[ensID][t][trackId+1] for t in eachindex(solEns[ensID].t)].!==0.0)
##
params = Dict(
    :N => 200000,
    :α => 1,
    :s => 0.2,
    # :σ => 0.25,
    :T => 100,
    :μ => 1,
    :sType => "fixed",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 500
trackerFitness = 0.2

##
trackerBirthTime = 5
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(params; runs,trackerVariant=(trackerBirthTime,trackerFitness), getArgs=true, noDiffusion=true)
sizeTrackedAverage5_t = CompetitiveSelection.averageTrackedVariant(solEns, simArgs)
##
trackerBirthTime = 10
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(params; runs,trackerVariant=(trackerBirthTime,trackerFitness), getArgs=true, noDiffusion=true)
sizeTrackedAverage10_t = CompetitiveSelection.averageTrackedVariant(solEns, simArgs)
##
trackerBirthTime = 15
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(params; runs,trackerVariant=(trackerBirthTime,trackerFitness), getArgs=true, noDiffusion=true)
sizeTrackedAverage15_t = CompetitiveSelection.averageTrackedVariant(solEns, simArgs)
##
rscale=0.9
fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="variant size",
    title="sAv = "*string(params[:s])
    # title="t₀ = "*string(trackerBirthTime)*
    # ", sAv = "*string(params[:s])*
    # ", sTracker = "*string(trackerFitness)
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 5, params[:α], trackerFitness, params[:s], params[:μ], params[:N])).(_t),
    label="theory ev",
    linestyle=:dash,
    color=Cycled(1),
)
lines!(_t, sizeTrackedAverage5_t, 
    label="sims average",
    # linestyle=:dash,
    color=Cycled(2),
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 10, params[:α], trackerFitness, params[:s], params[:μ], params[:N])).(_t),
    # label="theory ev" 
    linestyle=:dash,
    color=Cycled(1),
)
lines!(_t, sizeTrackedAverage10_t,
    # label="sims average"
    # linestyle=:dash,
    color=Cycled(2),
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, 15, params[:α], trackerFitness, params[:s], params[:μ], params[:N])).(_t),
    # label="theory ev" 
    linestyle=:dash,
    color=Cycled(1),
)
lines!(_t, sizeTrackedAverage15_t,
    # label="sims average"
    # linestyle=:dash,
    color=Cycled(2),
)
axislegend(position=:lt)
display(fig)
saveFolder="Figures/Presentations/"
save(saveFolder*"variantTrajectories_ThVSims.png", fig, px_per_unit = 2)


##
fig = Figure(resolution=(rscale*800,rscale*600), fontsize=18)
Axis(fig[1,1],
    xlabel="time",
    ylabel="variant size",
    title="t₀ = "*string(trackerBirthTime)*
    ", sAv = "*string(params[:s])*
    ", sTracker = "*string(trackerFitness)
)
lines!(
    _t, (t->Theory.fIndFitnessWithAvCompetition(t, trackerBirthTime, params[:α], trackerFitness, params[:s], params[:μ], params[:N])).(_t),
    label="theory ev"
)
# lines!(
#     _t, (t->Theory.fCompetition(t, trackerBirthTime, params[:α], params[:s], params[:μ], params[:N])).(_t),
#     label="theory fixed"
# )
lines!(_t, sizeTrackedAverage_t, label="sims average")
axislegend(position=:lt)
display(fig)
# save(saveFolder*"trackedSimVTheory_u"*string(trackerBirthTime)*"_s"*string(trackerFitness)*".png", fig, px_per_unit = 2)
##

ensID = 51
trackId = simArgs[ensID, :trackerID]
lines([solEns[ensID][t][trackId] for t in eachindex(solEns[ensID].t)])

##
include("../src/competitiveSelection.jl")
using .CompetitiveSelection
##
paramsPop = Dict(
    :N => 100000,
    :α => 1,
    # :s => 0.12,
    :s => 0.05,
    :σ => 0.2,
    :T => 105,
    # :μ => 0.14,
    :μ => 0.24,
    :sType => "gamma",
    # :sType => "free",
    :q => 0.3,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 5000;

@time solEns = CompetitiveSelection.evolvePopSim(paramsPop, runs=runs)
@time CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop, tMeasure=(50,90), timeLimit=0.4)