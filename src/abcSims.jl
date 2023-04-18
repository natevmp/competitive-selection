module ABCSims

include("../src/competitiveSelection.jl")
using .CompetitiveSelection
include("../src/analysisTools.jl")
using .AnalysisTools
using DifferentialEquations, Random, Distributions, Distances
using GpABC
# using ApproxBayes


## ============================ ABC functions: GpABC ===========================

function sampleSimTimepoint(sim::RODESolution, t::Real, nSamples::Int; fMin::Float64=0., fMax::Float64=1.)
    _vid = findall((f -> f>=fMin && f<=fMax).(sim(t))) |> shuffle!
    # make sure number of samples to draw does not exceed number of available variants
    length(_vid) < nSamples && (nSamples=length(_vid))
    vaf_vid = sim(t)[_vid[1:nSamples]]
    return vaf_vid
end

function modelSim(sweepParams::AbstractVector, ctrlParams::Dict)
    #Set parameters for sim: pars[1]=s, pars[2]=σ, pars[3]=μ
    params = ctrlParams[:params]
    params[:s] = sweepParams[1]
    params[:σ] = sweepParams[2]
    params[:μ] = sweepParams[3]

    # run simulations
    # solEns, _ = CompetitiveSelection.evolvePopSim(params; runs=ctrlParams[:simRuns], noDiffusion=false)
    model = CompetitiveSelection.GammaSelectionModel(ctrlParams[:params][:s], ctrlParams[:params][:σ], ctrlParams[:params][:q])
    solEns, _ = CompetitiveSelection.evolvePopSim(params, model; runs=ctrlParams[:simRuns], noDiffusion=false)
    
    # sample simulations at desired time point
    vaf_vid = Vector{Float64}(undef, ctrlParams[:simRuns]*ctrlParams[:nSamples])
    vCur = 1
    vNext = 1
    for (i,sim) in enumerate(solEns)
        f_vidSim = sampleSimTimepoint(sim, params[:T], ctrlParams[:nSamples]; fMin=2*ctrlParams[:fMin], fMax=2*ctrlParams[:fMax])./2
        vNext += length(f_vidSim)
        vaf_vid[vCur:vNext-1] .= f_vidSim
        vCur = vNext
    end

    # construct size distribution
    _f, n_f = AnalysisTools.sizeDistDens(@view vaf_vid[1:vCur-1]; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax])

    if ctrlParams[:cumulativeDist]
        # construct cumulative size distribution
        return hcat([sum(@view n_f[1:i]) for i in eachindex(n_f)])
    else
        return hcat(n_f)
    end
end

function runABC(nData_f, ctrlParams)
    # set priors
    #pars[1]=s, pars[2]=σ, pars[3]=μ
    priors = [
        Uniform(ctrlParams[:sBounds][1], ctrlParams[:sBounds][2]),  # s
        Uniform(ctrlParams[:σBounds][1], ctrlParams[:σBounds][2]),  # σ
        Uniform(ctrlParams[:μBounds][1], ctrlParams[:μBounds][2]),  # μ
    ]

    if ctrlParams[:metric]=="chisquare"
        summary_statistic = "keep_all"
        # distance_function = Distances.chisq_dist
        distance_function = (x,y)->evaluate(Distances.chisq_dist, x, y)
    elseif ctrlParams[:metric]=="euclidean"
        summary_statistic = "keep_all"
        distance_function = Distances.euclidean
    elseif ctrlParams[:metric]=="kolmogorov-smirnov"
        summary_statistic = "keep_all"
        distance_function = kolmogorovSmirnovDist
    elseif ctrlParams[:metric]=="mean"
        summary_statistic = "mean"
        distance_function = Distances.euclidean
    else
        error("Incorrect distance specified in `ctrlParams`")
    end

    # perform ABC
    # SimulatedABCRejection(
    #     hcat(nData_f),
    #     (par -> modelSim(par, ctrlParams)),
    #     priors,
    #     ctrlParams[:threshold],
    #     ctrlParams[:nParticles];
    #     summary_statistic=summary_statistic,
    #     distance_function=distance_function,
    #     max_iter=ctrlParams[:maxIter]
    # )
    SimulatedABCSMC(
        hcat(nData_f),
        (par -> modelSim(par, ctrlParams)),
        priors,
        ctrlParams[:threshold],
        ctrlParams[:nParticles];
        summary_statistic=summary_statistic,
        distance_function=distance_function,
        max_iter=ctrlParams[:maxIter]
    )
end


# ## ======================== ABC function: ApproxBayes =====================




# """
# """
# function runABC_MJW(nData_f, ctrlParams)

#     if ctrlParams[:metric]=="chisquare"
#         summary_statistic = "keep_all"
#         # distance_function = Distances.chisq_dist
#         distance_function = (x,y)->evaluate(Distances.chisq_dist, x, y)
#     elseif ctrlParams[:metric]=="euclidean"
#         summary_statistic = "keep_all"
#         distance_function = Distances.euclidean
#     elseif ctrlParams[:metric]=="kolmogorov-smirnov"
#         summary_statistic = "keep_all"
#         distance_function = kolmogorovSmirnovDist
#     elseif ctrlParams[:metric]=="mean"
#         summary_statistic = "mean"
#         distance_function = Distances.euclidean
#     else
#         error("Incorrect distance specified in `ctrlParams`")
#     end

#     # set priors
#     #pars[1]=s, pars[2]=σ, pars[3]=μ
#     priors = [
#         Uniform(ctrlParams[:sBounds][1], ctrlParams[:sBounds][2]),  # s
#         Uniform(ctrlParams[:σBounds][1], ctrlParams[:σBounds][2]),  # σ
#         Uniform(ctrlParams[:μBounds][1], ctrlParams[:μBounds][2]),  # μ
#     ]

#     simFunc(sweepParams, constants, nData_f) = 
#         distance_function(nData_f, modelSim(sweepParams, ctrlParams)), 1

#     setup = ABCRejection(
#         simFunc, #simulation function
#         3, # number of parameters
#         0.1, #target ϵ
#         Prior(priors); # Prior for each of the parameters
#         nparticles = ctrlParams[:nParticles],
#         maxiterations = ctrlParams[:maxIter], #Maximum number of iterations before the algorithm terminates
#     )

#     runabc(setup, nData_f)
# end

## ========================== Distance Functions ==========================

function kolmogorovSmirnovDist(X_t, Y_t)
    maxDist = 0.
    for i in eachindex(X_t)
        abs(X_t[i] - Y_t[i]) > maxDist && (maxDist = abs(X_t[i] - Y_t[i]))
    end
    return maxDist
end







end