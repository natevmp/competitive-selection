include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/abcRejection.jl")
using .ABCRejection
include("../src/dataStructuring.jl")
using.DataStructuring
include("../src/analysisTools.jl")
using .AnalysisTools

using JLD2, CSV, DataFrames
using CairoMakie
using Distributions

## ======== Parameter values ========

tBounds = (65,71)
modelParams = Dict{Symbol,Union{Int64,Float64,String}}(
    :T => 100,
    :N => 200000,
    :α => 1,
    :sType => "gamma",
    :q => 0,
)
hiddenParams = (
    s=0.1,
    σ=0.07,
    μ=3.,
)
ctrlParams = Dict(
    :tBounds => tBounds,
    :simRuns => 500,
    :fBins => 25,
    :sBins => 50,
    :sFitBounds => (-0.5, 1),
    :sQError => 0.98,
    :nSamples => 2, # number of variants to sample per patient (sim) when performing logistic fit
    :fMin => 0.0002, #the minimum size appearing in the data
    :fMax => 0.5,
    :coverage => 1000, # sequencing coverage to use in sampling
    :sBounds => (0.01, 1.), # ABC particle parameter bounds
    :σBounds => (0.005, 0.1), # ABC particle parameter bounds
    :μBounds => (0.1, 8), # ABC particle parameter bounds
    :params => modelParams, # model paraters to use for simulations
    :cumulativeDist => false,
    :metric => "chisquare",
    # :threshold => [4E-3,],
    :thresholds => (0.1,0.08),
    :nVarsBounds => (1,25),
    :nParticles => 1,
    :maxIter => 20,
    :nTimeSteps => 4, # number of time points per variant to sample when performing logistic fit
    :tStep => 2.5,
    :fixPar => Symbol[],
)

## -------- Get Simulated Data metrics --------

@time simResults = CompetitiveSelection.runModelSim(hiddenParams, ctrlParams)
nData_f, nData_s, _ = simResults
dataMetrics = (nData_f, nData_s)

## ================ Get simulation particle results ================

particle_pid=load("./Data/SimResults/Particles6571/particle_pid.jld2", "particle_pid")

## ----------- acceptance thresholds ------------

include("../src/abcRejection.jl")
using .ABCRejection


# ctrlParams[:thresholds] = (1, (0.015, 0.01))
ctrlParams[:thresholds] = (0.0008, (0.0012, 0.011))
ctrlParams[:nVarsBounds] = (1,25)

_s = range(ctrlParams[:sFitBounds][1], ctrlParams[:sFitBounds][2], length=ctrlParams[:sBins])
accepted_pid = ABCRejection.testABCParticles(
    (dataMetrics, simResults, thresholds)->CompetitiveSelection.compareDataVSim2(dataMetrics, simResults, thresholds, _s),
    particle_pid,
    dataMetrics,
    ctrlParams[:thresholds],
    # simRes->CompetitiveSelection.checkConstraints(simRes, ctrlParams[:nVarsBounds])
)
println("number of accepted particles: ", sum(accepted_pid))

params_aid = ABCRejection.acceptedParams(particle_pid, accepted_pid)
val_aid_param = Array{Float64,2}(undef, length(params_aid), length(params_aid[1]))
for (aid, params) in enumerate(params_aid)
    for parid in 1:length(params)
        val_aid_param[aid, parid] = params[parid]
    end
end

## ======================== Plots ========================
set_theme!(Theme(fontsize = 20))

## -------- histograms + 2d accepted particle scatterplots --------

figABCRes = Figure(resolution=(1400,800))
bounds = (s=ctrlParams[:sBounds], σ=ctrlParams[:σBounds], μ=ctrlParams[:μBounds])
for (i, par) in enumerate([:s, :σ, :μ])
    Axis(figABCRes[1,i],
        xlabel=string(par),
        ylabel=(i==1 ? "counts" : ""),
    )
    hist!(val_aid_param[:,i])
    xlims!(bounds[par][1], bounds[par][2])
end
# μ -> s
Axis(figABCRes[2,3], 
    xlabel="μ",
    ylabel="s",
    title="accepted particles"
)
xlims!(bounds[:μ]...)
ylims!(bounds[:s]...)
scatter!(val_aid_param[:,3], val_aid_param[:,1])
# ylims!(0.03, 0.17)

# s -> σ
Axis(figABCRes[2,2], 
    xlabel="σ",
    ylabel="s",
    title="accepted particles"
)
xlims!(bounds[:σ]...)
ylims!(bounds[:s]...)
scatter!(val_aid_param[:,2], val_aid_param[:,1])

# μ -> σ
Axis(figABCRes[2,1], 
    xlabel="σ",
    ylabel="μ",
    title="accepted particles"
)
xlims!(bounds[:σ]...)
ylims!(bounds[:μ]...)
scatter!(val_aid_param[:,2], val_aid_param[:,3])
Label(
    figABCRes[0, :],
    "true parameters: s="*string(hiddenParams[:s])*", σ="*string(hiddenParams[:σ])*", μ="*string(hiddenParams[:μ]),
    fontsize = 30
)
display(figABCRes)

save("./Figures/DataFits/ABC/abcResults_SimulatedData.png", figABCRes)

## ======================== accepted size distributions ========================

dists_tid = ABCRejection.measureABCParticles(
    (dataMetrics, simResults)->CompetitiveSelection.compareDataVSimError(dataMetrics, simResults, _s),
    particle_pid,
    dataMetrics
)
sizeDistEr_tid = [dists[1] for dists in dists_tid]

using ColorSchemes

minSDistError = minimum(sizeDistEr_tid[.!accepted_pid][1:end-1])
maxSDistError = maximum(sizeDistEr_tid[.!accepted_pid])


(maxSDistError - minSDistError)

distanceToColor(e) = (e - minSDistError)*1/(maxSDistError-minSDistError)
ColorSchemes.Blues_3[2]

# _f = _fEdges[1:end-1] .+ Float64(_fEdges.step)/2
_f = range(ctrlParams[:fMin], ctrlParams[:fMax], length=ctrlParams[:fBins])
figSizeDist = Figure(fontsize=20)
Axis(figSizeDist[1,1],
    xlabel="variant allele frequency",
    ylabel="density of variants",
    xscale=log10,
    yscale=log10,
    title="true parameters: s="*string(hiddenParams[:s])*", σ="*string(hiddenParams[:σ])*", μ="*string(hiddenParams[:μ])
)
minimum(particle_pid[4].simResults[1])
ylims!(2E-5,0.21E-1)
xlims!(_f[1], _f[end])
lines!(_f, particle_pid[1].simResults[1], color=:grey65, label="rejected particle")
for pid in findall(.!accepted_pid)[2:end]
    lines!(_f, particle_pid[pid].simResults[1],
    color=:grey65,
    )
end
lines!(_f, particle_pid[accepted_pid][1].simResults[1], color=ColorSchemes.Blues_3[2], label="accepted particle")
for pid in findall(accepted_pid)[2:end]
    lines!(_f, particle_pid[pid].simResults[1], color=ColorSchemes.Blues_3[2])
end
scatterlines!(_f, nData_f, label="simulated data")
axislegend(position=:rt)
display(figSizeDist)

figname="sizeDistABCParticles_SimVSim.png"
save("./Figures/DataFits/ABC/"*figname, figSizeDist)

##
#region -------- accepted innate fitness distributions --------

gammaMeanStdToShapeScale(mean, std) = ((mean/std)^2, std^2/mean)
_s = range(0,0.25,length=300)
prob_s_pid = Array{Float64,2}(undef, length(_s), length(particle_pid))
for (i,particle) in enumerate(particle_pid)
    (α, θ) = gammaMeanStdToShapeScale(particle.paramSet[:s], particle.paramSet[:σ])
    prob_s_pid[:,i] = pdf.(Gamma(α, θ), _s)
end

rsc=0.8
figFitDist = Figure(fontsize=24, resolution=(rsc*800,rsc*600))
Axis(
    figFitDist[1,1],
    xlabel="s",
    ylabel="pdf"
    # limits=((_s[1],_s[end]), nothing),
)
xlims!(_s[1],_s[end])
ylims!(0, nothing)
for pid in findall(.!accepted_pid)
    lines!(_s, prob_s_pid[:,pid],
    color=(:grey60,0.2),
    )
end
for pid in findall(accepted_pid)
    lines!(_s, prob_s_pid[:,pid],
    color=(Makie.wong_colors()[1],0.2),
    )
end
display(figFitDist)

figname="innateFitnessDistABCParticles_SimVSim.png"
save("./Figures/DataFits/ABC/"*figname, figFitDist)

#endregion

##
#region -------- logistic fitness Distributions --------

_s = range(ctrlParams[:sFitBounds][1], ctrlParams[:sFitBounds][2], length=ctrlParams[:sBins])
rsc = 0.8
figFitness = Figure(fontsize=24, resolution=(rsc*800,rsc*600))
Axis(figFitness[1,1],
    xlabel="logistic fitness s*",
    ylabel="density",
    # title="Ages "*string(tBounds[1])*"-"*string(tBounds[2])*": particle no. "*string(pid),
)
lines!(
    _s, particle_pid[.!accepted_pid][1].simResults[2],
    color=:grey65,
    label="rejected particles",
)
for (i,particle) in enumerate(particle_pid[.!accepted_pid])
    lines!(
        _s, particle.simResults[2],
        color=:grey65,
    )
end
lines!(
    _s, particle_pid[accepted_pid][1].simResults[2],
    color=ColorSchemes.Blues_3[2],
    label="accepted particles",
)
for (i,particle) in enumerate(particle_pid[accepted_pid])
    lines!(
        _s, particle.simResults[2],
        color=ColorSchemes.Blues_3[2],
    )
end
scatterlines!(_s, dataMetrics[2], label="data")
axislegend(position=:rt)
xlims!(-0.5,0.8)
ylims!(0,0.012)
display(figFitness)

figname="logisticFitnessDistABCParticles_SimVSim.png"
save("./Figures/DataFits/ABC/"*figname, figFitness)

#endregion
