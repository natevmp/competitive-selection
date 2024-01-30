
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
using Glob, Statistics, ProgressMeter

## ========================================
#region = Parameter values

# ctrlParams = load("./Data/SimResults/Particles5D_7076/parameters.jld2", "ctrlParams")
dataDir = "C:/Users/monper02/Documents/Datasets/CompetitiveSelectionSimResults/Particles5D_7076/"
ctrlParams = load(dataDir*"parameters.jld2", "ctrlParams")

#endregion

## ----------------------------------------
#region - Get Data metrics

xTotData_t, nData_t, _t = load("./Data/variantNumbersSizesMitchell.jld2", "xTot_t", "n_t", "_t")

#endregion

## ========================================
#region = Get simulation particle results

SAVEPARTICLEVECTOR = false
if SAVEPARTICLEVECTOR
    fname_pid = glob("particle*", dataDir)
    typeParticle = typeof(load(fname_pid[1], "particle"))
    particle_pid = Vector{typeParticle}(undef, length(fname_pid))
    @showprogress for (pid, fname) in enumerate(fname_pid)
        particle_pid[pid] = load(fname, "particle")
    end
    save(
        simDir*"/simParticlesVariantNumberSizes.jld2",
        "particle_pid", particle_pid
    )
end

particle_pid=load(
    dataDir*
    "simParticlesVariantNumberSizes.jld2",
    "particle_pid",
)

#endregion

## ----------------------------------------
#region - Test parameter prior distributions


for par in keys(particle_pid[1].paramSet)
    display(hist([particle.paramSet[par] for particle in particle_pid], bins=100, axis=(xlabel=string(par),)))
end


#endregion


## ----------------------------------------
#region - acceptance thresholds
include("../src/abcRejection.jl")
using .ABCRejection
include("../src/competitiveSelection.jl")
using .CompetitiveSelection
##
propertynames(particle_pid[1])
particle_pid[2].simResults[4]

length(particle_pid)
eltype(particle_pid)

particle_pid[3].simResults[4]
xTotData_t

function constraintsF(simResults)
    _, _, n_t, _ = simResults
    return n_t[3] <= 1
end

ctrlParams[:thresholds] = (n=57.315,x=1.816) # strong fits logistic and size dist
acceptedN_pid = ABCRejection.testABCParticles(
    CompetitiveSelection.compareVariantNumbers,
    particle_pid,
    nData_t,
    ctrlParams[:thresholds][:n],
    # constraintsF
)
acceptedX_pid = ABCRejection.testABCParticles(
    CompetitiveSelection.compareVariantSizes,
    particle_pid,
    xTotData_t,
    ctrlParams[:thresholds][:x],
    # constraintsF
)
println("number of accepted particles variant numbers: ", sum(acceptedN_pid))
println("number of accepted particles variant sizes: ", sum(acceptedX_pid))
##
params_aNid = ABCRejection.acceptedParams(particle_pid, acceptedN_pid)
val_aNid_param = Array{Float64,2}(undef, length(params_aNid), length(params_aNid[1]))
for (aid, params) in enumerate(params_aNid)
    for parid in 1:length(params)
        val_aNid_param[aid, parid] = params[parid]
    end
end
params_aXid = ABCRejection.acceptedParams(particle_pid, acceptedX_pid)
val_aXid_param = Array{Float64,2}(undef, length(params_aXid), length(params_aXid[1]))
for (aid, params) in enumerate(params_aXid)
    for parid in 1:length(params)
        val_aXid_param[aid, parid] = params[parid]
    end
end

#endregion

## ----------------------------------------
#region - Plots
using ColorSchemes
using Makie.Colors
set_theme!(Theme(fontsize = 24))
#endregion

## ----------------------------------------
#region - 2d accepted particle scatterplots

SAVEFIGURES = true

bounds = (s=ctrlParams[:sBounds], σ=ctrlParams[:σBounds], μ=ctrlParams[:μBounds])

labels = (
    s="s",
    σ="σ",
    μ="μ",
    N="N",
    τ="τ",
)

ctrlParams
parlimits = (
    s=ctrlParams[:sBounds],
    σ=ctrlParams[:σBounds],
    μ=ctrlParams[:μBounds],
    N=ctrlParams[:NBounds],
    τ=ctrlParams[:τBounds],
)

parPairs = [
    (1, 2),
    (1, 3),
    (1, 4),
    (1, 5),
    (2, 3),
    (2, 4),
    (2, 5),
    (3, 4),
    (3, 5),
    (4, 5),
]

figN = Figure(size=(2200,1000))
for (i, pair) in enumerate(parPairs)
    row = (i-1)/5+1 |>floor|>Int
    col = (i-1)%5+1
    Axis(figN[row,col],
        xlabel=labels[pair[1]],
        ylabel=labels[pair[2]]
    )
    scatter!(val_aNid_param[:,pair[1]], val_aNid_param[:,pair[2]])
end
for (i, par) in enumerate([:s, :σ, :μ, :N, :τ])
    Axis(figN[3,i],
        xlabel=string(par),
        ylabel="counts",
    )
    hist!(val_aNid_param[:,i])
    xlims!(parlimits[par]...)
end
Label(figN[0, :], "ABC - variant numbers", fontsize = 40)
display(figN)
SAVEFIGURES && save("Figures/abc5dVariantNumbers.svg", figN)

figX = Figure(size=(2200,1000))
for (i, pair) in enumerate(parPairs)
    row = (i-1)/5+1 |>floor|>Int
    col = (i-1)%5+1
    Axis(figX[row,col],
        xlabel=labels[pair[1]],
        ylabel=labels[pair[2]]
    )
    scatter!(val_aXid_param[:,pair[1]], val_aXid_param[:,pair[2]])
end
for (i, par) in enumerate([:s, :σ, :μ, :N, :τ])
    Axis(figX[3,i],
        xlabel=string(par),
        ylabel="counts",
    )
    hist!(val_aXid_param[:,i])
    xlims!(parlimits[par]...)
end
Label(figX[0, :], "ABC - total variant fractions", fontsize = 40)
display(figX)
SAVEFIGURES && save("Figures/abc5dVariantFractions.svg", figX)

#endregion

## ----------------------------------------
#region detectable clones with age

figDetectableClones = Figure(fontsize=28)
Axis(
    figDetectableClones[1,1],
    xlabel="age",
    ylabel="average number\nof detected clones"
)
scatterlines!(_t, particle_pid[.!acceptedN_pid][1].simResults[3],
    color=(:grey75,0.2),
    label="accepted particle"
)
for (pid, particle) in enumerate(particle_pid[.!acceptedN_pid])
    nVars_t = particle.simResults[3]
    lines!(_t, nVars_t, 
        color=(:grey75,0.2),
    )
end
scatterlines!(_t, particle_pid[acceptedN_pid][1].simResults[3],
    color=(Makie.wong_colors()[1],0.55),
    label="accepted particle"
)
for (pid, particle) in enumerate(particle_pid[acceptedN_pid][2:end])
    nVars_t = particle.simResults[3]
    scatterlines!(_t, nVars_t,
        color=(Makie.wong_colors()[1],0.55),
    )
end
scatterlines!(
    _t, nData_t,
    color=(Makie.wong_colors()[2],1),
    linewidth=2,
    label="data",
)
axislegend(position=:lt)
display(figDetectableClones)

figname = "abcDetectableClonesAccepted.png"
save("./Figures/"*figname, figDetectableClones, px_per_unit=2)

#endregion

## ----------------------------------------
#region total fraction of expanded clones with age

figTotalFrac = Figure(fontsize=28)
Axis(
    figTotalFrac[1,1],
    xlabel="age",
    ylabel="average variant fraction\nin stem cell pool"
)
scatterlines!(_t, particle_pid[.!acceptedX_pid][1].simResults[4],
    color=(:grey75,0.2),
    label="rejected particle"
)
for (pid, particle) in enumerate(particle_pid[.!acceptedX_pid][2:end])
    nVars_t = particle.simResults[4]
    lines!(_t, nVars_t, 
        color=(:grey75,0.2),
    )
end
scatterlines!(_t, particle_pid[acceptedX_pid][1].simResults[4],
color=(Makie.wong_colors()[1],0.55),
label="accepted particle"
)
for (pid, particle) in enumerate(particle_pid[acceptedX_pid][2:end])
    nVars_t = particle.simResults[4]
    scatterlines!(_t, nVars_t,
        color=(Makie.wong_colors()[1],0.55),
    )
end
scatterlines!(
    _t, xTotData_t,
    color=(Makie.wong_colors()[2],1),
    linewidth=2,
    label="data"
)
ylims!(0,0.4)
axislegend(position=:lt)
display(figTotalFrac)

figname = "abcVariantFracAccepted.png"
save("./Figures/"*figname, figTotalFrac, px_per_unit=2)

#endregion