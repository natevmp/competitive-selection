include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/theory.jl")
using .Theory
using Statistics
using CairoMakie


## ----------------------------------------
#region - Run sims

alpha = 5
params = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 100000,
    :α => alpha,
    :sType => "gamma",
    :q => 1.,
    :μ => 3.,
    # :s => 0.1 / alpha,
    # :σ => 0.04 / alpha,
    :s => 2,
    :sMax => 50,
    :σ => 0.04,
)
@time solEns, _ = CompetitiveSelection.evolvePopSim(
    params;
    runs=100,
    noDiffusion=false,
);
3;
#endregion

## ----------------------------------------
#region - Average total variant fraction

_t = 0:(length(solEns[1])-1)
nTotAv_t = Vector{Float64}(undef, length(_t))
nClones_t = Vector{Float64}(undef, length(_t))
detectThresh = 0.005
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
# nTotTheory_t = (t->Theory.compCoverage(t, params[:α], params[:s], params[:μ], params[:N])).(_t)

fig = Figure(fontsize=24)
Axis(fig[1,1],
    xlabel="time",
    ylabel="number of clones above threshold",
)
# lines!(nTotTheory_t)
# lines!(nClones_t, label="clones alive")
lines!(nClonesThresh_t, label="α = "*string(alpha))
axislegend(position=:lt)
# ylims!(0,25)
display(fig)

#endregion