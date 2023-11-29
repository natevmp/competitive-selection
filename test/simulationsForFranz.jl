include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using Statistics
using CairoMakie


## ----------------------------------------
#region - Run sims
alpha=1
detectThresh = 0.01
params = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 100000,
    :α => alpha,
    :sType => "gamma",
    :q => 0.,
    :μ => 2,
    :s => 0.1/alpha,
    :σ => 0.04,
)
@time solEns, _ = CompetitiveSelection.evolvePopSim(
    params;
    runs=500,
    noDiffusion=false,
);
4
#endregion

## ----------------------------------------
#region - Average total variant fraction

_t = 0:(length(solEns[1])-1)
nTotAv_t = Vector{Float64}(undef, length(_t))
nClones_t = Vector{Float64}(undef, length(_t))
nClonesThresh_t = Vector{Float64}(undef, length(_t))
for tInd in eachindex(_t)
    nTot = 0
    nClones = 0
    nClonesThresh = 0
    for sol in solEns
        nTot += sum(sol[tInd])
        nClones += sum(sol[tInd].>0)
        nClonesThresh += sum(sol[tInd].>detectThresh)
    end
    nTotAv_t[tInd] = nTot / length(solEns)
    nClones_t[tInd] = nClones / length(solEns)
    nClonesThresh_t[tInd] = nClonesThresh / length(solEns)
end


nTotTheory_t = (t->Theory.compCoverage(t, params[:α], params[:s], params[:μ], params[:N])).(_t)

fig = Figure(fontsize=26)
Axis(fig[1,1],
    xlabel="time",
    ylabel="number of clones",
)
# lines!(nTotTheory_t)
# lines!(nClones_t, label="clones alive")
# hlines!(ones(length(_t)).*3, color=:black)
lines!(nClonesThresh_t,
    label="alpha = "*string(params[:α])
    *"μ = "*string(params[:μ])
    ,
    )
axislegend(position=:lt)
ylims!(0,20)
xlims!(20,80)
display(fig)

#endregion

## ----------------------------------------
#region - save data

using CSV, DataFrames

df = DataFrame(
    "t"=>_t,
    "Expected total variant fraction" => nTotTheory_t,
    "Average total variant fraction" => nTotAv_t,
    "Average number of existing clones" => nClones_t,
    "Average number of clones above threshold "*string(detectThresh) => nClonesThresh_t,
)

# CSV.write("totalVariantFracTime.csv", df)
##
sol = solEns[6]
sol.u[1]

n_t_vid = Array{Float64,2}(undef, length(sol), length(sol.u[1]))
for t in 1:size(n_t_vid,1)
    n_t_vid[t,:] .= sol.u[t]
end

time = 100
sum(n_t_vid[time,:])
sum(n_t_vid[time,:].>0)
sum(n_t_vid[time,:].>0.01)

fig = Figure()
Axis(fig[1,1])
lines!([sum(n_t_vid[t,:]) for t in 1:params[:T]])
display(fig)

#endregion

## ----------------------------------------
#region - Plot number of variants alive with time

nAlive_t = [sum(n_t_vid[t,:].>0) for t in 1:size(n_t_vid,1)]

fig = Figure()
Axis(fig[1,1], xlabel="age", ylabel="number of live clones")
scatterlines!(nAlive_t)
display(fig)

#endregion


