
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

ctrlParams = load("./Data/SimResults/Particles5D_7076/parameters.jld2", "ctrlParams")

# tBounds = (70,76)
# modelParams = Dict{Symbol,Union{Int64,Float64,String}}(
#     :T => 100,
#     :N => 200000,
#     :α => 1,
#     :sType => "gamma",
#     :q => 1,
# )
# ctrlParams = Dict(
#     :tBounds => tBounds,
#     :simRuns => 500,
#     :fBins => 25,
#     :sBins => 50,
#     :sFitBounds => (-0.5, 1),
#     :sQError => 0.95,
#     :nSamples => 2, # number of variants to sample per patient (sim) when performing logistic fit
#     :nTimeSteps => 4, # number of time points per variant to sample when performing logistic fit
#     :tStep => 2.5, # time between sample time points when performing logistic fit 
#     :fMin => 0.0002, #the minimum size appearing in the data
#     :fMax => 0.5,
#     :coverage => 1000, # sequencing coverage to use in sampling
#     :fixPar => Symbol[],
#     :sBounds => (0.01, 1.), # ABC particle parameter bounds
#     :σBounds => (0.005, 0.1), # ABC particle parameter bounds
#     :μBounds => (0.1, 15), # ABC particle parameter bounds
#     :params => modelParams, # model parameters to use for simulations
#     # -------- rejection parameters --------
#     :nParticles => 1,
#     :maxIter => 20,
#     :metric => "chisquare",
#     :cumulativeDist => false,
#     :thresholds => (0.1,(0.1,0.1)),
#     :tEarly => 40,
#     :tLate => 70,
#     :nVarsBoundEarly => (0,2),
#     :nVarsBoundLate => (1,25),
# )
# tEarly, (lBoundEarly, uBoundEarly), tLate, (lBoundLate, uBoundLate)

#endregion

## ----------------------------------------
#region - Get Data metrics

# df = DataFrame(CSV.File("./Data/ALLvariants_exclSynonymous_Xadj.txt"))
# DataStructuring.structureData!(df)
# dfVid = DataStructuring.sortDataPerVariant(df)
# DataStructuring.fitAllModels!(dfVid, ctrlParams[:params][:N], ["positive logistic", "negative logistic", "constant"]; errorQuantile=0.99)
# save(
#     "./Data/dfDataFittedVid"*agesString*".jld2", Dict(
#         "dfVid"=>dfVid
#     )
# )
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

folder = "./Data/simResults/Particles5D_7076/"

# fname_pid = glob("particle*", folder)

# typeParticle = typeof(load(fname_pid[1], "particle"))
# particle_pid = Vector{typeParticle}(undef, length(fname_pid))
# @showprogress for (pid, fname) in enumerate(fname_pid)
#     particle_pid[pid] = load(fname, "particle")
# end

# save(
#     folder*"/simParticles.jld2",
#     "particle_pid", particle_pid
# )

particle_pid=load(
    folder*
    "simParticles.jld2",
    "particle_pid",
)

#endregion



##

particle = particle_pid[1]
propertynames(particle)

hist([particle.paramSet[1] for particle in particle_pid], bins=100)

length(particle.simResults)
particle.simResults[1]


## ----------------------------------------
#region - acceptance thresholds
include("../src/abcRejection.jl")
using .ABCRejection

propertynames(particle_pid[1])
particle_pid[1].simResults

length(particle_pid)
eltype(particle_pid)

# ctrlParams[:thresholds] = (Inf, (Inf, Inf)) # all particles
# ctrlParams[:thresholds] = (0.0019, (1, 1))
# ctrlParams[:thresholds] = (0.004, (0.001, 0.0157)) # maximize fit to logistic fitness
# ctrlParams[:thresholds] = (0.0012, (0.003, 0.02)) # maximize fit to size dist
ctrlParams[:thresholds] = (0.0012, (0.001, 0.018)) # intermediate fit logistic and size dist
# ctrlParams[:thresholds] = (0.0015, (0.001, 0.03)) # intermediate fit logistic and size dist
ctrlParams[:tEarly] = 30
ctrlParams[:tLate] = 70
ctrlParams[:nVarsBoundEarly] = (0,3)
ctrlParams[:nVarsBoundLate] = (5,25)
ctrlParams[:nVarsBoundEarly] = (0,100)
ctrlParams[:nVarsBoundLate] = (0,100)

# ctrlParams[:thresholds] = (5, (1, 1))
# ctrlParams[:tEarly] = 50
# ctrlParams[:tLate] = 70

_s = _sEdges[1:end-1] .+ Float64(_sEdges.step)/2
accepted_pid = ABCRejection.testABCParticles(
    (dataMetrics, simResults, thresholds)->CompetitiveSelection.compareDataVSim2(dataMetrics, simResults, thresholds, _s),
    particle_pid,
    dataMetrics,
    ctrlParams[:thresholds],
    simRes->CompetitiveSelection.checkConstraintsDetectableVariants(
        simRes,
        ctrlParams[:tEarly],
        ctrlParams[:nVarsBoundEarly],
        ctrlParams[:tLate],
        ctrlParams[:nVarsBoundLate]
    )
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

bounds = (s=ctrlParams[:sBounds], σ=ctrlParams[:σBounds], μ=ctrlParams[:μBounds])
labels = (
    s="s",
    σ="σ",
    μ="μ",
    N="N",
    α="α",
)
ctrlParams
parlimits = (
    s=ctrlParams[:sBounds],
    σ=ctrlParams[:σBounds],
    μ=ctrlParams[:μBounds],
    N=ctrlParams[:NBounds],
    α=ctrlParams[:αBounds],
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
end
for (i, par) in enumerate([:s, :σ, :μ, :N, :α])
    Axis(fig[3,i],
        xlabel=string(par),
        ylabel="counts",
    )
    hist!(val_aid_param[:,i])
    xlims!(parlimits[par]...)
end
display(fig)
save("Figures/abc5d.png", fig)

##
# for (i, par) in enumerate([:s, :σ, :μ, :N, α])
#     Axis(figABCRes[i,1],
#         xlabel=labels[par],
#         # ylabel=(i==1 ? "counts" : ""),
#         title=(i==1 ? "marginal distribution" : ""),
#     )
#     hist!(val_aid_param[:,i],
#         # color=:gray,
#         color=Makie.wong_colors()[1],
#         strokewidth = 2,
#         strokecolor = Makie.wong_colors()[1],
#     )
#     xlims!(parlimits[par]...)
# end

# s -> σ
Axis(figABCRes[1,2],
    xlabel="s",
    ylabel="σ",
    title="accepted particles"
)
xlims!(parlimits[:s]...)
ylims!(parlimits[:σ]...)
scatter!(val_aid_param[:,1], val_aid_param[:,2])
# ylims!(0.03, 0.17)

# σ -> μ
Axis(figABCRes[2,2], 
    xlabel="σ",
    ylabel="μ",
)
xlims!(parlimits[:σ]...)
ylims!(parlimits[:μ]...)
scatter!(val_aid_param[:,2], val_aid_param[:,3])

# μ -> s
Axis(figABCRes[3,2], 
    xlabel="μ",
    ylabel="s",
)
xlims!(parlimits[:μ]...)
ylims!(parlimits[:s]...)
scatter!(val_aid_param[:,3], val_aid_param[:,1])

display(figABCRes)

# save("./Figures/ManuscriptDrafts/Figure 3/abcResults.svg", figABCRes)

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

# save("./Figures/ManuscriptDrafts/Figure 3/fitDistCompare.png", figFitness, px_per_unit=2)

#endregion

## ----------------------------------------
#region - size distribution accepted/rejected

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

figname="sizeDistABCParticles.png"
# save("./Figures/DataFits/ABC/"*figname, figSizeDist)
# save("./Figures/ManuscriptDrafts/Figure 3/"*figname, figSizeDist, px_per_unit=2)

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
    color=(Makie.wong_colors()[1],0.55),
    )
end
display(figFitDist)
#
figname = "fitnessDistributionABCResults.png"
# save("./Figures/ManuscriptDrafts/"*figname, figFitDist, px_per_unit=2)

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
# save("./Figures/ManuscriptDrafts/"*figname, figDetectableClones, px_per_unit=2)

#endregion
