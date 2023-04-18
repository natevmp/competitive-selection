include("../src/competitiveSelection.jl")
using .CompetitiveSelection
using Statistics
using CairoMakie

##

paramsPop = Dict(
    :N => 200000,
    :α => 1,
    :s => 0.12,
    :T => 100,
    :μ => 10,
    :σ => 0.1,
    :sType => "fixed",
)
runs = 1000

# @time sol = CompetitiveSelection.evolvePop(paramsPop)
@time solEns = CompetitiveSelection.evolvePop(paramsPop, runs=runs)

##

nVars_sim_t = Array{Float64}(undef, runs, length(solEns[1]))
for simId in 1:runs
    for tInd in 1:paramsPop[:T]
        detectableVariants_vid = solEns[simId][tInd] .> 0.01
        nVars_sim_t[simId,tInd] = sum(detectableVariants_vid)
    end
end
nVarsAv_t = vec(mean(nVars_sim_t,dims=1))
nVarsStd_t = vec(std(nVars_sim_t,dims=1))

_t = 0:paramsPop[:T]
fig0 = Figure()
Axis(fig0[1,1], xlabel="time", ylabel="number of variants > 1%")
lines!(_t, nVarsAv_t)
band!(_t, nVarsAv_t + nVarsStd_t, nVarsAv_t - nVarsStd_t)
display(fig0)

##

_f, nV_t_f = CompetitiveSelection.variantsSizeDistribution(solEns, xMin=0.01, xMax=0.26, bins=50)

##

fig1 = Figure()
Axis(
    fig1[1,1],
    # yscale=Makie.pseudolog10,
    xlabel="variant frequency",
    ylabel="density of variants",
    title="Variants above 1%",
)
for tInd in [20,40,60,80]
    lines!(_f, nV_t_f[tInd,:], label="age: "*string(tInd))
end
axislegend()
display(fig1)








##  
xTot_t = [sum(sol[t][sol[t] .> 0.01]) for t in 1:length(sol)]

fig = Figure()
Axis(
    fig[1,1],
    # yscale=Makie.pseudolog10,
    xlabel="age",
    ylabel="variant size",
    title="selection type = "*paramsPop[:sType]*"; s = "*string(paramsPop[:s])*"; μ = "*string(paramsPop[:μ]),
    # xscale=identity,
)
# ylims!(1E-6, 1)
for i in 1:20
    lines!([sol[t][i] for t in 1:length(sol)])
end
lines!(xTot_t, color=:black, linestyle=:dash)
display(fig)

##
lines(xTot_t)




##

function distbin(x_vid::Vector{T} where T<:Real; bins::Int=25, x0::Real=0, xMax=1.)
    _xEdges = range(x0, xMax, length=bins+1)
    dx = _xEdges[2]-_xEdges[1]
    _f = _xEdges[1:end-1] .+ dx
    nV_f = zeros(Int64, bins)
    for x in x_vid
        k = ((x-x0)/dx |> floor |> Int) + 1
        nV_f[k] += 1
    end
    return nV_f
end


##

distbin([0.1, 0.05, 0.6,0.6,0.61, 0.8, 0.9], bins=5)

