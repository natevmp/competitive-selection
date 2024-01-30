
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

simDir = "C:/Users/monper02/Documents/Datasets/CompetitiveSelectionSimResults/Particles5D_7076/"
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
ctrlParams
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

##
particle_pid[1].paramSet

sOverTau_pid = [particle.paramSet[:s]/particle.paramSet[:τ] for particle in particle_pid]
fig = Figure()
Axis(fig[1,1],
    xlabel="s/tau",
    ylabel="counts"
)
hist!(sOverTau_pid, bins=100)

# xlims!(0,0.5)
display(fig)



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
thresholdIntermediate = (0.001, (0.0015, 0.02)) # intermediate fit logistic and size dist 
thresholdIntermediate = (0.0014, (0.001, 0.02)) # intermediate fit logistic and size dist 
thresholdFitnessDist = (0.1, (0.0006, 0.017)) # maximize fit to fitness distribution
thresholdSizeDist = (0.0005, (1, 1)) # maximize fit to size distribution
ctrlParams[:thresholds] = thresholdFitnessDist
# ctrlParams[:tEarly] = 30
# ctrlParams[:tLate] = 70
# ctrlParams[:nVarsBoundEarly] = (0,3)
# ctrlParams[:nVarsBoundLate] = (14,40)
# ctrlParams[:nVarsBoundEarly] = (0,100)
# ctrlParams[:nVarsBoundLate] = (0,100)

_s = _sEdges[1:end-1] .+ Float64(_sEdges.step)/2
accepted_pid = ABCRejection.testABCParticles(
    (dataMetrics, simResults, thresholds)->CompetitiveSelection.compareDataVSim2(dataMetrics, simResults, thresholds, _s),
    particle_pid,
    dataMetrics,
    ctrlParams[:thresholds],
    # simRes->CompetitiveSelection.checkConstraintsDetectableVariants(
    #     simRes,
    #     ctrlParams[:tEarly],
    #     ctrlParams[:nVarsBoundEarly],
    #     ctrlParams[:tLate],
    #     ctrlParams[:nVarsBoundLate]
    # )
)
println("number of accepted particles: ", sum(accepted_pid))

params_aid = ABCRejection.acceptedParams(particle_pid, accepted_pid)
val_aid_param = Array{Float64,2}(undef, length(params_aid), length(params_aid[1]))
for (aid, params) in enumerate(params_aid)
    for (parid, par) in enumerate(keys(params))
        # if par==:s||par==:σ
        if false
            val_aid_param[aid, parid] = params[par]/params[:τ]
        else
            val_aid_param[aid, parid] = params[par]
        end
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

bounds = (s=ctrlParams[:sBounds], σ=ctrlParams[:σBounds], μ=ctrlParams[:μBounds])
labels = (
    s="s/τ",
    σ="σ/τ",
    μ="μ",
    N="N",
    τ="τ",
)
ctrlParams
parlimitsFull = (
    # s=ctrlParams[:sBounds] ./ 
        # (ctrlParams[:τBounds][2],ctrlParams[:τBounds][1]),
    s=ctrlParams[:sBounds],
    # σ=ctrlParams[:σBounds] ./ 
    #     (ctrlParams[:τBounds][2],ctrlParams[:τBounds][1]),
    σ=ctrlParams[:σBounds],
    μ=ctrlParams[:μBounds],
    N=ctrlParams[:NBounds],
    τ=ctrlParams[:τBounds],
)
parLimitsPlot = (
    s=(0,0.3),
    σ=(0,0.1),
    μ=(0,40),
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

fig = Figure(resolution=(2200,1000))
for (i, pair) in enumerate(parPairs)
    row = (i-1)/5+1 |>floor|>Int
    col = (i-1)%5+1
    Axis(fig[row,col],
        xlabel=labels[pair[1]],
        ylabel=labels[pair[2]]
    )
    scatter!(val_aid_param[:,pair[1]], val_aid_param[:,pair[2]])
    xlims!(parlimitsFull[pair[1]][1], parlimitsFull[pair[1]][2])
    ylims!(parlimitsFull[pair[2]][1], parlimitsFull[pair[2]][2])
end
for (i, par) in enumerate([:s, :σ, :μ, :N, :τ])
    Axis(fig[3,i],
        xlabel=labels[par],
        ylabel="counts",
    )
    hist!(val_aid_param[:,i], bins=40)
    xlims!(parlimitsFull[par]...)
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
figureTitle = figureTitle*"; "*string(size(val_aid_param,1))*" accepted particles"
Label(fig[0, :], figureTitle, fontsize = 40)
display(fig)
filename = "Figures/abc5d"*fitName*".png"
SAVEFIG && save(filename, fig)

#endregion

## ----------------------------------------
#region - Posterior distributions

# tau
figTauPost = Figure()
Axis(figTauPost[1,1])
hist!(val_aid_param[:,5], bins=15)
vlines!([mean(val_aid_param[:,5])], color=:red, label="mean: "*string(round(mean(val_aid_param[:,5]),digits=2)))
vlines!([median(val_aid_param[:,5])], color=:blue, label="median: "*string(round(median(val_aid_param[:,5]),digits=2)))
xlims!(0, 10)
axislegend(position=:rt)
display(figTauPost)


#endregion

## ----------------------------------------
#region - division rate normalized fitness parameter

figFitnessNormed = Figure(fontsize=24)

# for (i, par) in enumerate([:s, :σ, :μ, :N, :τ])
Axis(figFitnessNormed[1,1],
    xlabel="s/τ",
    ylabel="counts",
)
hist!(val_aid_param[:,1] ./ val_aid_param[:,5], bins=40)
xlims!(0,0.2)
display(figFitnessNormed)

#endregion

## ----------------------------------------
#region - s*tau inference:

hist(val_aid_param[:,1]./val_aid_param[:,5])

#endregion

## ----------------------------------------
#region - Comparing accepted particles: logistic fitness distribution
using ColorSchemes

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
filename = "Figures/fitnessDist_"*fitName*".png"
SAVEFIG && save(filename, figFitness)

#endregion

## ----------------------------------------
#region - size distribution accepted/rejected
rsc=1.
ColorSchemes.Blues_3[2]
_f = _fEdges[1:end-1] .+ Float64(_fEdges.step)/2
figSizeDist = Figure(fontsize=24, resolution=(rsc*800,rsc*600))
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

filename = "Figures/SizeDist_"*fitName*".png"
SAVEFIG && save(filename, figSizeDist)

#endregion

## ----------------------------------------
#region - true fitness distribution

gammaMeanStdToShapeScale(mean, std) = ((mean/std)^2, std^2/mean)
_s = range(0,0.25,length=300)
prob_s_pid = Array{Float64,2}(undef, length(_s), length(particle_pid))
for (i,particle) in enumerate(particle_pid)
    (α, θ) = gammaMeanStdToShapeScale(particle.paramSet[:s], particle.paramSet[:σ])
    prob_s_pid[:,i] = pdf.(Gamma(α, θ), _s)
end

#
rsc = 1
figFitDist = Figure(fontsize=24, size=(rsc*800,rsc*600))
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
    color=(Makie.wong_colors()[1],0.55),
    )
end
display(figFitDist)
#
figname = "fitnessDistributionABCResults.png"
SAVEFIG && save("./Figures/ManuscriptDrafts/"*figname, figFitDist, px_per_unit=2)

#endregion

## ----------------------------------------
#region - true fitness distribution - Normalized

gammaMeanStdToShapeScale(mean, std) = ((mean/std)^2, std^2/mean)
_s = range(0,0.25,length=300)
prob_s_pid = Array{Float64,2}(undef, length(_s), length(particle_pid))
for (i,particle) in enumerate(particle_pid)
    (α, θ) = gammaMeanStdToShapeScale(particle.paramSet[:s]/particle.paramSet[:τ], particle.paramSet[:σ]/particle.paramSet[:τ])
    prob_s_pid[:,i] = pdf.(Gamma(α, θ), _s)
end

#

figFitDist = Figure(fontsize=28, size=(rsc*800,rsc*600))
Axis(
    figFitDist[1,1],
    xlabel="s/τ",
    ylabel="pdf"
    # limits=((_s[1],_s[end]), nothing),
)
xlims!(_s[1],_s[end])
ylims!(0, 150)
for pid in findall(.!accepted_pid)
    lines!(_s, prob_s_pid[:,pid],
    color=(:grey60,0.2),
    )
end
for pid in findall(accepted_pid)
    lines!(_s, prob_s_pid[:,pid],
    color=(Makie.wong_colors()[1],0.4),
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

