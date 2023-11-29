include("../src/competitiveSelection.jl")
using .CompetitiveSelection
using Statistics
using CairoMakie
using ElasticArrays, StatsBase, LsqFit

##

myTheme = Theme(fontsize = 20)
set_theme!(myTheme)

##

paramsPop = Dict(
    :N => Int(2E4),
    :α => 1,
    # :s => 0.12,
    :s => 0.2,
    :σ => 0.15,
    :T => 105,
    :μ => 0.5,
    # :μ => 0,
    :sType => "gamma",
    # :sType => "free",
    :q => 0,    #double hit probability: 0 = no double hits; 1 = double hits succeed with probability 1.
)
runs = 5000

##
# @time sol = CompetitiveSelection.evolvePop(paramsPop)
@time solEns, simArgs = CompetitiveSelection.evolvePopSim(paramsPop, runs=runs);
4+4

##
measureTime=50
# @time λ_tVaf_cid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop; t=measureTime)
@time λ_tVaf_cid = CompetitiveSelection.fitSamplesGrowth(solEns, paramsPop)
fig = scatter(λ_tVaf_cid[1,:], λ_tVaf_cid[2,:])
xlims!(0,100)
display(fig)
##
thresh = 2
fig2 = Figure()
Axis(
    fig2[1,1],
    xlabel="measured fitness",
    ylable="density of variants",
    title="t=$(measureTime); s=$(paramsPop[:s]); σ=$(paramsPop[:σ]); q=$(paramsPop[:q])",
)
hist!(λ_tVaf_cid[2, λ_tVaf_cid[2,:].<=thresh]; bins=50, normalization=:pdf)
xlims!(-1,thresh)
display(fig2)
##

simID = 1
t = 80
_t, _cid, vaf_t_cid = CompetitiveSelection.sampleSim(solEns[simID], t=t)
println(vaf_t_cid)

##
Nf = paramsPop[:N]
vaf_t = vaf_t_cid[:,3]
println(_t)
println(vaf_t)
fLogistic(t, t0, r) = 1 / ( 1 + (1 - 1/Nf)/(1/Nf) * exp(-r*(t-t0)) )
fAntiLog(t, t0, r, b0) = b0*Nf / ( Nf + exp(r*(t-t0)) - 1 )
@. modelLogistic(t, p) = fLogistic(t, p[1], p[2])
@. modelAntiLog(t, p) = fAntiLog(t, p[1], p[2], p[3])
corr = corspearman(_t, vaf_t)
fit = let 
    if corr >= 0
        curve_fit(modelLogistic, _t, vaf_t, [1., 0.2], lower=[-Inf,0], upper=[_t[1], +Inf])
    else
        curve_fit(modelAntiLog, _t, vaf_t, [1.,-0.1, 0.05], lower=[-Inf, -Inf, 0], upper=[_t[1], 0, 1])
    end
end
println(fit.param)
##


λ_tVaf_cid = ElasticArray{Float64}(undef,2,0)
for i in eachindex(_cid)
    fit = CompetitiveSelection.fitLogGrowth(_t, vaf_t_cid[:, i], paramsPop[:N])
    append!(λ_tVaf_cid, [t, fit.param[2]])
end


## ---------------------- Plot individual trajectories -------------------

simID = 9
xTot_t = [sum(solEns[simID][t][solEns[simID][t] .> 0.01]) for t in 1:length(solEns[simID])]
fig = Figure()
Axis(
    fig[1,1],
    # yscale=Makie.pseudolog10,
    xlabel="patient age",
    ylabel="variant size",
    title="selection type = "*paramsPop[:sType]*"; s = "*string(paramsPop[:s])*"; σ = "*string(paramsPop[:σ])*"; μ = "*string(paramsPop[:μ]),
    # xscale=identity,
)
# ylims!(1E-6, 1)
for i in 1:length(solEns[simID][end])
    lines!([solEns[simID][t][i] for t in 1:length(solEns[simID])])
end
# lines!(xTot_t, color=:black, linestyle=:dash)
display(fig)

save("trajectoryFitness.png", fig, pt_per_unit=1)

##


##
_f, nV_t_f = CompetitiveSelection.variantsSizeDistribution(solEns, xMin=0.01, xMax=1., bins=100, normalized=true)
sum(nV_t_f[paramsPop[:T], :])
# println("fraction of fixated clones: "*string(round(nV_t_f[end][end],digits=3)))

##
vfs = CompetitiveSelection.VFreqspace(paramsPop[:N], 501, 0)
@time CompetitiveSelection.evolveSizeDist(vfs, paramsPop, paramsPop[:T])
# nSpl = Spline1D(vfs._f, vfs.n_f; k=2)
# quadgk(nSpl, 0,1)[1]
# @time CompetitiveSelection.evolveDriftPDE(vfs, paramsPop, paramsPop[:T])
##
fig2 = Figure(
    resolution=(0.8*800,0.8*700),
)
fig2[1,1] = Axis(fig2[1,1],
    # xscale=log10,
    # yscale=log10,
    xlabel="f",
    ylabel="density of variants",
)
# for tt in [55,70,85,100]
# lines!(_f, nV_t_f[tt,:], label="age "*string(tt))
# end
lines!(_f[1:20], nV_t_f[end,1:20] / nV_t_f[end,1], label="simulations")
lines!(vfs._f[2:end-1], vfs.n_f[2:end-1] / vfs.n_f[findfirst(vfs._f.>_f[1])], label="pde")
xlims!(0.01,0.2)
ylims!(0.01,0.9)
axislegend()

display(fig2)

##

# paramsPop = Dict(
#     :N => Int(2E4),
#     :α => 10,
#     :s => 0.1,
#     :T => 90,
#     :μ => 1,
#     :σ => 0.1,
#     :sType => "fixed",
#     :n₀ => 100,
# )

# vfs = CompetitiveSelection.VFreqspace(paramsPop[:N], 100, 0)
# @time CompetitiveSelection.evolveExistingVariants(vfs, paramsPop, paramsPop[:T])

# ##

# fig3 = Figure()
# fig3[1,1] = Axis(fig3[1,1], xscale=identity, yscale=log10)
# xlims!(1/paramsPop[:N],0.9)
# ylims!(1E-10, 0.04)
# lines!(vfs._f[2:end-1], vfs.n_f[2:end-1])
# # lines!(vfs._f, vfs.n_f)
# display(fig3)




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
    dx = _xEdges[   2]-_xEdges[1]
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


##
nVars_sim_t = Array{Float64}(undef, runs, length(solEns[1]))
for simId in 1:runs
    for tInd in 1:paramsPop[:T]
        detectableVariants_vid = solEns[simId][tInd] .> 0.02
        nVars_sim_t[simId,tInd] = sum(detectableVariants_vid)
    end
end
nVarsAv_t = vec(mean(nVars_sim_t,dims=1))

##

function variantsList(solEns)
    tLen = length(solEns[1])
    nSims = length(solEns)
    x_T_Vid = Vector{Vector{Float64}}(undef, tLen)
    for tInd in 1:tLen
        x_Sim_Vid = Array{Array{Float64}}(undef, nSims)
        for simId in 1:nSims
            x_vid = solEns[simId][tInd]
            x_Sim_Vid[simId] = x_vid
        end
        x_T_Vid[tInd] = vcat(x_Sim_Vid...)
    end
    return x_T_Vid
end

x_T_Vid = variantsList(solEns)

##
using StatsBase, LinearAlgebra
histData = LinearAlgebra.normalize(fit(Histogram, x_T_Vid[90], range(0.01,0.26, length=50)), mode=:pdf)
lines(histData.weights)

##

hist(x_T_Vid[60][x_T_Vid[90] .> 0.02], normalization=:pdf, bins=25)

##

testf(a,b,c) = (a,b)

a,b,_ = testf(1,2,3)

##

include("../src/competitiveSelection.jl")
using .CompetitiveSelection


##
# @time nVarsAv_t = CompetitiveSelection.variantsAboveThreshold(solEns, 0.01; nSamples=500)
@time nVarsAv_t = CompetitiveSelection.detectableVariantsSampled(solEns, 500, threshold=0.01)
lines(nVarsAv_t)

##

particle = load("particleTest.jld2", "particle")

propertynames(particle)
particle.simResults[3]
particle.paramSet