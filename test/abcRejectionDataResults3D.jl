
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

simDir = "C:/Users/monper02/Documents/Datasets/CompetitiveSelectionSimResults/Particles3D_7076/"
ctrlParams = load(simDir*"parameters.jld2", "ctrlParams")

#endregion

## ----------------------------------------
#region - Get Data metrics

CALCDATAMETRICS = false

if CALCDATAMETRICS
    df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
    DataStructuring.structureData!(df)
    dfVid = DataStructuring.sortDataPerVariant(df)
    DataStructuring.fitAllModels!(dfVid, ctrlParams[:params][:N], ["positive logistic", "negative logistic", "constant"]; errorQuantile=0.99)
    save(
        "./Data/dfDataFittedVid"*agesString*".jld2", Dict(
            "dfVid"=>dfVid
        )
    )
end
dfVid = load("./Data/dfDataFittedVid.jld2", "dfVid")
tBounds = ctrlParams[:tBounds]
agesString = string(tBounds[1])*string(tBounds[2])
x_sid = [x for x in Iterators.flatten([dfVidCur[:vaf_t] for dfVidCur in eachrow(dfVid)])]
t_sid = [t for t in Iterators.flatten([dfVidCur[:_t] for dfVidCur in eachrow(dfVid)])]
xMask_sid = x_sid.>0
_fEdges, nData_f = AnalysisTools.sizeDistDens(
    x_sid[xMask_sid], t_sid[xMask_sid], tBounds; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax]
)
sFit_vid = dfVid[dfVid[!,:goodFit],:γ]
_sEdges, nData_s = AnalysisTools.distbin(sFit_vid, 
    bins=ctrlParams[:sBins], 
    xMin=ctrlParams[:sFitBounds][1],
    xMax=ctrlParams[:sFitBounds][2],
)

dataMetrics = (nData_f, nData_s)

#endregion

## ========================================
#region = Get simulation particle results

SAVEPARTICLEVECTOR = false
if SAVEPARTICLEVECTOR
    fname_pid = glob("particle*", simDir)

    typeParticle = typeof(load(fname_pid[1], "particle"))
    particle_pid = Vector{typeParticle}(undef, length(fname_pid))
    @showprogress for (pid, fname) in enumerate(fname_pid)
        particle_pid[pid] = load(fname, "particle")
    end

    save(
        simDir*"/simParticlesSizeFitnessDist.jld2",
        "particle_pid", particle_pid
    )
end
particle_pid=load(
    simDir*
    "simParticlesSizeFitnessDist.jld2",
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

propertynames(particle_pid[1])
particle_pid[1].simResults

length(particle_pid)
eltype(particle_pid)

thresholdAllParticles = (Inf, (Inf, Inf)) # all particles
thresholdIntermediate = (0.0009, (0.0015, 0.018)) # intermediate fit logistic and size dist 
thresholdFitnessDist = (0.1, (0.001, 0.0175)) # maximize fit to fitness distribution
thresholdSizeDist = (0.0005, (1, 1)) # maximize fit to size distribution
ctrlParams[:thresholds] = thresholdSizeDist

_s = _sEdges[1:end-1] .+ Float64(_sEdges.step)/2
accepted_pid = ABCRejection.testABCParticles(
    (dataMetrics, simResults, thresholds)->CompetitiveSelection.compareDataVSim2(dataMetrics, simResults, thresholds, _s),
    particle_pid,
    dataMetrics,
    ctrlParams[:thresholds],
)
println("number of accepted particles: ", sum(accepted_pid))

params_aid = ABCRejection.acceptedParams(particle_pid, accepted_pid)
val_aid_param = Array{Float64,2}(undef, length(params_aid), length(params_aid[1]))
for (aid, params) in enumerate(params_aid)
    for parid in 1:length(params)
        val_aid_param[aid, parid] = params[parid]
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

SAVEFIG = false

labels = (
    σ="σ",
    μ="μ",
    τ="τ",
)
parlimits = (
    σ=ctrlParams[:σBounds],
    μ=ctrlParams[:μBounds],
    τ=ctrlParams[:τBounds],
)

parPairs = [
    (1, 2),
    (1, 3),
    (2, 3),
]

fig = Figure(size=(1300,800))
for (i, pair) in enumerate(parPairs)
    row = (i-1)/3+1 |>floor|>Int
    col = (i-1)%3+1
    Axis(fig[row,col],
        xlabel=labels[pair[1]],
        ylabel=labels[pair[2]]
    )
    scatter!(val_aid_param[:,pair[1]], val_aid_param[:,pair[2]])
    xlims!(parlimits[pair[1]]...)
    ylims!(parlimits[pair[2]]...)
end
for (i, par) in enumerate([:σ, :μ, :τ])
    Axis(fig[2,i],
        xlabel=string(par),
        ylabel="counts",
    )
    hist!(val_aid_param[:,i])
    xlims!(parlimits[par]...)
end
fitName, figureTitle = 
    if ctrlParams[:thresholds]==thresholdFitnessDist
        "FitnessDist", "Single ABC - fitness dist"
    elseif ctrlParams[:thresholds]==thresholdSizeDist
        "SizeDist", "Single ABC - size dist"
    elseif ctrlParams[:thresholds]==thresholdIntermediate
        "Intermediate", "Single ABC - size dist + fitness dist"
    else
        error()
    end
Label(fig[0, :], figureTitle, fontsize = 40)
display(fig)
filename = "Figures/abc3d"*fitName*".png"
SAVEFIG && save(filename, fig)

#endregion

## ----------------------------------------
#region - s-sigma relationship posterior

particle = particle_pid[1]
particle.paramSet[:τ]
particle.paramSet[:τ]*ctrlParams[:sFixed]

figSigmaS = Figure(fontsize=24)
Axis(figSigmaS[1,1],
    xlabel="η",
    ylabel="σ",
)
scatter!(val_aid_param[:,3].*ctrlParams[:sFixed], val_aid_param[:,1])
# xlims!(parlimits[pair[1]]...)
# ylims!(parlimits[pair[2]]...)
display(figSigmaS)

figSigmaSHist = Figure(fontsize=24)
Axis(figSigmaSHist[1,1],
    xlabel="σ/η",
)
hist!(val_aid_param[:,1] ./ (val_aid_param[:,3].*ctrlParams[:sFixed]))
# xlims!(parlimits[pair[1]]...)
# ylims!(parlimits[pair[2]]...)
display(figSigmaSHist)

figSigmaTauHist = Figure(fontsize=24)
Axis(figSigmaTauHist[1,1],
    xlabel="σ/τ",
)
hist!(val_aid_param[:,1] ./ (val_aid_param[:,3]))
# xlims!(parlimits[pair[1]]...)
# ylims!(parlimits[pair[2]]...)
display(figSigmaTauHist)

#endregion

## ----------------------------------------
#region - Comparing accepted particles: logistic fitness distribution
using ColorSchemes

rsc = 0.8
figFitness = Figure(fontsize=24, size=(rsc*800,rsc*600))
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
filename = "Figures/fitnessDist3D_"*fitName*".png"
SAVEFIG && save(filename, figFitness)

#endregion

## ----------------------------------------
#region - size distribution accepted/rejected
rsc=1.
ColorSchemes.Blues_3[2]
_f = _fEdges[1:end-1] .+ Float64(_fEdges.step)/2
figSizeDist = Figure(fontsize=24, size=(rsc*800,rsc*600))
Axis(figSizeDist[1,1],
    xlabel="variant allele frequency",
    ylabel="density of variants",
    xscale=log10,
    yscale=log10,
    # title="Ages "*string(tBounds[1])*"-"*string(tBounds[2])
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
scatterlines!(_f, nData_f, label="Fabre data")
axislegend(position=:rt)
display(figSizeDist)

filename = "Figures/SizeDist3D_"*fitName*".png"
SAVEFIG && save(filename, figSizeDist)

#endregion

## ----------------------------------------
#region - true fitness distribution

gammaMeanStdToShapeScale(mean, std) = ((mean/std)^2, std^2/mean)
_s = range(0,0.2,length=300)
prob_s_pid = Array{Float64,2}(undef, length(_s), length(particle_pid))
for (i,particle) in enumerate(particle_pid)
    (α, θ) = gammaMeanStdToShapeScale(ctrlParams[:sFixed], particle.paramSet[:σ] / particle.paramSet[:τ])
    prob_s_pid[:,i] = pdf.(Gamma(α, θ), _s)
end

#

figFitDist = Figure(fontsize=24, size=(rsc*800,rsc*600))
Axis(
    figFitDist[1,1],
    xlabel="s",
    ylabel="pdf"
    # limits=((_s[1],_s[end]), nothing),
)
xlims!(_s[1],_s[end])
ylims!(0, 100)
# for pid in findall(.!accepted_pid)
#     lines!(_s, prob_s_pid[:,pid],
#     color=(:grey60,0.2),
#     )
# end
for pid in findall(accepted_pid)
    lines!(_s, prob_s_pid[:,pid],
    color=(Makie.wong_colors()[1],0.55),
    )
end
display(figFitDist)
#
figname = "fitnessDistributionABCResults.png"
SAVEFIG && save("./Figures/ManuscriptDrafts/"*figname, figFitDist, px_per_unit=2)

#endregion

## ----------------------------------------
#region detectable clones with age

particle_pid[1].simResults[3]

figDetectableClones = Figure(fontsize=28)
Axis(
    figDetectableClones[1,1],
    xlabel="age",
    ylabel="average number of detected clones"
)
for (pid, particle) in enumerate(particle_pid[.!accepted_pid])
    nVars_t = particle.simResults[3]
    lines!(nVars_t, 
        color=(:grey60,0.2),
    )
end
for (pid, particle) in enumerate(particle_pid[accepted_pid])
    nVars_t = particle.simResults[3]
    lines!(nVars_t,
        color=(Makie.wong_colors()[1],0.55),
    )
end
display(figDetectableClones)

figname = "abcDetectableClonesAccepted.png"
SAVEFIG && save("./Figures/ManuscriptDrafts/"*figname, figDetectableClones, px_per_unit=2)

#endregion

