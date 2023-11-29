

"""
get size binned histogram from list of variant sizes.
"""
function distbin(x_vid; bins::Int=25, xMin::Real=0, xMax::Union{Nothing,Real}=nothing, normalized::Bool=true)
    isnothing(xMax) && (xMax = maximum(x_vid)+0.01)
    _xEdges = range(xMin, xMax, bins+1)
    dx = _xEdges[2]-_xEdges[1]
    _x = _xEdges[1:end-1] .+ dx/2
    inRange_vid = (x -> x>=xMin && x<xMax).(x_vid)   #apply thresholds
    n_x = zeros(Float64, bins)
    for x in x_vid[inRange_vid]
        k = ((x-xMin)/dx |> floor |> Int) + 1
        n_x[k] += 1
    end
    normalized && (n_x .*= dx/sum(inRange_vid))
    return _x, n_x
end

"""
Get sizes distribution of variants at multiple time points determined by `bins`, from ensemble of simulations.
"""
function variantsSizeDistributionTimeBinned(solEns::EnsembleSolution; bins=25, xMin::Real=0.01, xMax::Union{Nothing,Real}=nothing, normalized::Bool=true)
    tLen = length(solEns[1])
    nSims = length(solEns)
    nV_t_f = Array{Float64,2}(undef, tLen, bins)
    _f = Vector{Float64}(undef, bins)
    for tInd in 1:tLen
        x_Sim_Vid = Array{Array{Float64}}(undef, nSims)
        for simId in 1:nSims
            x_vid = solEns[simId][tInd]
            x_Sim_Vid[simId] = x_vid
        end
        xGen_vid = vcat(x_Sim_Vid...)
        _f, nV_t_f[tInd,:] = distbin(xGen_vid; bins, xMin, xMax, normalized)

    end
    return _f, nV_t_f
end

"""
Get size distribution of variants at multiple time points determined by `bins`, from ensemble of simulations.
This method uses the StatsBase Histogram type, but is much slower than the main implementation.
"""
function variantsSizeDistributionAlt(solEns::EnsembleSolution; bins=25, xMin::Real=0.01, xMax::Real=1, normalized::Bool=true)
    normMode = normalized ? :pdf : :none
    tLen = length(solEns[1])
    nSims = length(solEns)
    nV_t_f = Array{Float64,2}(undef, tLen, bins)
    _edges = range(xMin, xMax, length=bins+1)
    _f = _edges[1:end-1] .+ (_edges[2]-_edges[1])/2
    for tInd in 1:tLen
        x_Sim_Vid = Array{Array{Float64}}(undef, nSims)
        for simId in 1:nSims
            x_vid = solEns[simId][tInd]
            x_Sim_Vid[simId] = x_vid
        end
        xGen_vid = vcat(x_Sim_Vid...)
        histData = LinearAlgebra.normalize(
            fit(Histogram, xGen_vid, _edges),
            mode=normMode
        )
        nV_t_f[tInd, :] = histData.weights
    end
    return _f, nV_t_f
end

function sampleIntFromFreq(vaf, nSamples)
    vaf<0 && (vaf=0)
    vaf>1 && (vaf=1)
    rand(Binomial(nSamples, vaf)) / nSamples
end

function variantsAboveThreshold(solEns::EnsembleSolution, f0::Real; nSamples::Union{Nothing,Real}=nothing)
    T = length(solEns[1])
    runs = length(solEns)
    nVarsAv_t = Vector{Float64}(undef, T)
    for tInd in eachindex(solEns[1])
        nVarsTot = 0
        for simId in eachindex(solEns)
            if !isnothing(nSamples)
                vafSampled_vid = (vaf->sampleIntFromFreq(vaf, nSamples)).(solEns[simId][tInd])
                detectableVariants_vid = vafSampled_vid .> f0
            else
                detectableVariants_vid = solEns[simId][tInd] .> f0
            end
            nVarsTot += sum(detectableVariants_vid)
        end
        nVarsAv_t[tInd] = nVarsTot / runs
    end
    return nVarsAv_t
end

function variantsAboveThreshold(solEns::EnsembleSolution, t::Real, f0::Real; nSamples::Union{Nothing,Real}=nothing)
    nVars_sim = Array{Float64}(undef, length(solEns))
    for simId in eachindex(solEns)
        if !isnothing(nSamples)
            vafSampled_vid = (vaf->sampleIntFromFreq(vaf, nSamples)).(solEns[simId](t))
            detectableVariants_vid = vafSampled_vid .> f0
        else
            detectableVariants_vid = solEns[simId](t) .> f0
        end
        nVars_sim[simId] = sum(detectableVariants_vid)
    end
    return mean(nVars_sim)
end


# function fitLogGrowth(_t, vaf_t, Nf)
#     @. modelLogistic(t, p) = fLogistic(t, p[1], p[2], 1/Nf, 1/2)
#     @. modelDecreasingLogistic(t, p) = fLogistic(t, p[1], p[2], p[3], 1/2)
#     corr = corspearman(_t, vaf_t)
#     fit = let 
#         if corr >= 0
#             curve_fit(modelLogistic, _t, vaf_t, [1., 0.2], lower=[-Inf,0], upper=[_t[1], +Inf])
#         else
#             curve_fit(modelDecreasingLogistic, _t, vaf_t, [1.,-0.1, 0.05], lower=[-Inf, -Inf, 0], upper=[_t[1], 0, 1])
#         end
#     end
#     return fit
# end

# function fitLogGrowthMLE(_t, vaf_t, nCov_t, Nf)

#     X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on
#     corr = corspearman(_t, vaf_t)
#     logprobX(X, t, n, b) = let 
#         if corr >= 0 # growing pop
#             logpdf(Binomial(n, fLogistic(t, b[1], b[2], 1/Nf, 0.5)), X)
#         else # shrinking pop
#             logpdf(Binomial(n, fLogistic(t, b[1], b[2], b[3], 0.5)), X)
#         end
#     end
#     logllh(b) = -([
#         logprobX(X_t[i], _t[i], nCov_t[i], b) for i in eachindex(_t)
#     ] |> sum)

#     # b[1]=t0, b[2]=γ, b[2]=x0
#     b0, lower, upper  = let 
#         if corr >= 0
#             [1., 0.05], [-1E6, 0], [_t[end], 10]
#         else
#             [1.,-0.1, vaf_t[1]], [-1E6, -Inf, 0], [_t[1], 0, 0.45]
#         end
#     end
#     #! debug
#     # println("X_t: ", X_t)
#     # for i in eachindex(_t) println("fLogistic: ", fLogistic(_t[i], b0[1], b0[2], 1/Nf, 1/2)) end
#     # for i in eachindex(_t) println("logprobX: ", logprobX(X_t[i], _t[i], nCov_t[i], b0) ) end
#     # for i in eachindex(_t) println("coverage: ", nCov_t[i] ) end
#     # println("llh: ", logllh(b0))
#     #! debug
#     alg = Fminbox(BFGS())
#     return optimize(logllh, lower, upper, b0, alg, Optim.Options(time_limit = 1.0))
#     # return optimize(logllh, lower, upper, b0)
# end

function sampleSimTimepoint(sim::RODESolution, t::Real, nSamples::Int, f0::Float64)
    _cid = findall(sim(t).>=f0) |> shuffle!
    # make sure number of samples to draw does not exceed number of available variants
    length(_cid) > nSamples && (nSamples=length(_cid))
    vaf_cid = sim(t)[_cid[1:nSamples]]
    return vaf_cid
end

function sampleSimVarTrajectories(sim::RODESolution; t::Union{Real,Nothing}=nothing, tStep=1., nTimeSteps::Int=4)

    isnothing(t) && (t = 10+rand()*80)
    _cid = findall(sim(t) .>= 0.01)
    _t = range(t; length=nTimeSteps, step=tStep)
    vaf_t_cid = Array{Float64,2}(undef, nTimeSteps, length(_cid))
    for (i,tt) in enumerate(_t)
        for (j,vaf) in enumerate(sim(tt)[_cid])
            vaf_t_cid[i,j] = vaf
        end
    end
    return _t, _cid, vaf_t_cid
end

# function fitSamplesGrowth(solEns::Union{EnsembleSolution,Vector{T}} where T, params::Dict; tMeasure::Union{Real,Tuple{Real, Real}}=(50,80), timeLimit=1.0)
#     Nf = params[:N]
#     tStep = 1.
#     nTimeSteps = 4
#     fit_tλ_cid = ElasticArray{Float64}(undef,2,0)
#     for sim in solEns
#         t = let
#             if length(tMeasure)<2
#                 tMeasure
#             else
#                 tMeasure[1]+rand()*(tMeasure[2]-tMeasure[1])
#             end
#         end # first measurement time point
#         _t, _cid, vaf_t_cid = sampleSim(sim; t, tStep, nTimeSteps)
#         for i in eachindex(_cid)
#             # fit = fitLogGrowth(_t, vaf_t_cid[:, i], Nf)
#             # append!(λ_tVaf_cid, [t, fit.param[2]])
#             vaf_t = vaf_t_cid[:, i]/2
#             vaf_t[1] > 0.45 && (vaf_t[1] = 0.45)
#             fit = DataStructuring.fitLogGrowthMLE!(_t, vaf_t, fill(1000, length(_t)), Nf; timeLimit)
#             _β = Optim.minimizer(fit)
#             append!(fit_tλ_cid, [t, _β[2]])
#         end
#     end
#     return fit_tλ_cid
# end

struct StepUniform <: Sampleable{Univariate,Continuous}
    edges::Vector{Float64}
    counts::Vector{Int64}
end

function Base.rand(rng::AbstractRNG, d::StepUniform)
    bin = StatsBase.sample(rng, range(1,length(d.counts)), Weights(d.counts))
    sample = d.edges[bin] + (d.edges[bin+1]-d.edges[bin])*rand()
    return sample
end

function fitSamplesGrowth(solEns::Union{EnsembleSolution,Vector{T}} where T, params::Dict; tMeasure::Union{Real,Tuple{Real, Real},Tuple{AbstractArray,AbstractArray}}=(50,80), timeLimit=1.0, nTimeSteps=4, tStep=1., errorQuantile=0.99)
    Nf = params[:N]
    # fit_tλ_cid = ElasticArray{Float64}(undef,2,0)
    dfVid = DataFrame(
        _t=AbstractVector{Real}[],
        vaf_t=Vector{Float64}[],
        cov_t=Vector{Int64}[]
    )

    t_sid = let 
        if length(tMeasure)<2
            nothing
        elseif eltype(tMeasure) <: Real
            tMeasure[1].+rand(length(solEns)).*(tMeasure[2]-tMeasure[1])
        else
            dist = StepUniform(tMeasure[1],tMeasure[2])
            rand(dist, length(solEns))
        end
    end

    for (j,sim) in enumerate(solEns)
        t = let
            if length(tMeasure)<2
                tMeasure
            else
                t_sid[j]
            end
        end # first measurement time point
        _t, _cid, vaf_t_cid = sampleSimVarTrajectories(sim; t, tStep, nTimeSteps)
        for i in eachindex(_cid)
            vaf_t = vaf_t_cid[:, i]/2
            cov_t = fill(1000, length(_t))
            push!(dfVid,
                Dict(
                    :_t=>_t,
                    :vaf_t=>vaf_t,
                    :cov_t=>cov_t
                )
            )
        end
    end
    # DataStructuring.fitVariantsData!(dfVid, Nf)
    DataStructuring.fitAllModels!(dfVid, Nf; errorQuantile)
    return dfVid
end

function averageTrackedVariant(solEns, simArgs)
    sizeTrackedvariant_sid_t = Array{Float64}(undef, length(solEns), length(solEns[1].t))
    for (sid, sol) in enumerate(solEns)
        trackID = simArgs[sid,:trackerID]
        for t in eachindex(sol.t)
            sizeTrackedvariant_sid_t[sid,t] = sol.u[t][trackID]
        end
    end
    return vec(mean(sizeTrackedvariant_sid_t, dims=1))
end

## ================================= ABC functions ===========================

"""
sample a variant allele frequency from a true frequency `x` with coverage `coverage`.
"""
function sampleFreqFromFreq(x, coverage)
    x<0 && (return 0.)
    x>1 && (return 1.)
    rand(Binomial(coverage, x))/coverage
end

"""
    sampleSimIDsTimepoint(sim::RODESolution, t::Real, nSamples::Int, coverage::Int; fMin::Float64=0., fMax::Float64=1.)

Sample variant trajectories from a single simulation with binomial sampling on the VAFs at times `_t`. Returns the variant id's `vid_j` of the observed variants, with a maximal length of `nSamples`.
"""
function sampleSimTimepoints(sim::RODESolution, _t::AbstractVector, nSamples::Int, coverage::Int; fMin::Float64=0., fMax::Float64=1.)
    # `x`∈[0,1]; not to be confused with `vaf`∈[0,0.5]
    x_vid_T = sim(_t)
    # The order of obsrvation, i.e. variant `vid` is observed `i-th`.
    vid_i = shuffle(1:length(x_vid_T[1]))
    j = 0   # index for accepted observed variants
    x_t_VidSampled = Vector{Vector{Float64}}(undef, nSamples)
    for vid in vid_i
        (x_vid_T[1][vid]<fMin || x_vid_T[1][vid]>fMax) && continue
        x_t = [x_vid[vid] for x_vid in x_vid_T]
        xObs_t = (x -> sampleFreqFromFreq(x, coverage)).(x_t)
        #check whether variant has enough nonzero observations; if not skip to next variant
        sum(xObs_t.==0)>=length(_t)/2 && continue
        j+=1
        x_t_VidSampled[j] = xObs_t
        j==nSamples && break
    end
    if j<nSamples
        return x_t_VidSampled[1:j]
    end
    return x_t_VidSampled
end

"""
    sampleSimTimepoint(sim::RODESolution, t::Real, nSamples::Int, coverage::Int; fMin::Float64=0., fMax::Float64=1.)

Samples variants from a single simulation with binomial sampling on the VAFs, at a single timepoint `t`.
"""
function sampleSimTimepoint(sim::RODESolution, childVid_child_Vid::Vector{Vector{T}} where T<:Integer, t::Real, nSamples::Int, coverage::Int; fMin::Float64=0., fMax::Float64=1.)
    _vidS = 1:length(sim(t)) |> shuffle
    # f_vid = sim(t) |> shuffle!
    f_vid = sim(t)
    vafObs_j = Float64[]
    j = 0
    for vidS in _vidS
        f = f_vid[vidS]
        # add size of children (if they exist)
        for childVid in childVid_child_Vid[vidS]
            f += f_vid[childVid]
        end
        (f<fMin || f>fMax) && continue
        # allele frequency of a variant is f/2
        vafObs = rand(Binomial(coverage, f/2)) / coverage
        vafObs==0 && continue
        push!(vafObs_j, vafObs)
        j+=1
        j==nSamples && break
    end
    # for f in f_vid
    #     (f<fMin || f>fMax) && continue
    #     # allele frequency of a variant is f/2
    #     vafObs = rand(Binomial(coverage, f/2)) / coverage
    #     vafObs==0 && continue
    #     push!(vafObs_j, vafObs)
    #     j+=1
    #     j==nSamples && break
    # end
    return vafObs_j
end

function buildFamilyArray(parentVid_vid::Vector{T} where T<:Integer)
    childVid_child_Vid = [Int64[] for _ in eachindex(parentVid_vid)]
    for (vid, parId) in enumerate(parentVid_vid)
        if parId==0 continue end
        push!(childVid_child_Vid[parId], vid)
    end
    return childVid_child_Vid
end

function sizeDistSims(tMeasure, solEns, parentVid_vid_Sid::Vector{Vector{T}} where T<:Integer, ctrlParams)
    # -------- sample sims --------
    vaf_vidS = Vector{Float64}(undef, ctrlParams[:simRuns]*ctrlParams[:nSamples])
    vCur = 1
    vNext = 1
    for (sid,sol) in enumerate(solEns)
        childVid_child_Vid = buildFamilyArray(parentVid_vid_Sid[sid])
        f_vidSimSample = sampleSimTimepoint(
                sol, childVid_child_Vid, tMeasure, ctrlParams[:nSamples], ctrlParams[:coverage];
                fMin=2*ctrlParams[:fMin], fMax=2*ctrlParams[:fMax]
            )
        vNext += length(f_vidSimSample)
        vaf_vidS[vCur:vNext-1] .= f_vidSimSample
        vCur = vNext
    end
    # -------- construct size distribution --------
    _f, n_f = AnalysisTools.sizeDistDens(@view vaf_vidS[1:vCur-1]; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax])
    if ctrlParams[:cumulativeDist]
        # construct cumulative size distribution
        n_f .= hcat([sum(@view n_f[1:i]) for i in eachindex(n_f)])
    else
        n_f .= hcat(n_f)
    end
    return n_f
end

function fitTrajectoriesFitness(tMeasure, solEns, ctrlParams)
    _t = range(tMeasure-(ctrlParams[:nTimeSteps]*ctrlParams[:tStep])/2, tMeasure+(ctrlParams[:nTimeSteps]*ctrlParams[:tStep])/2, length=ctrlParams[:nTimeSteps])
    if _t[end] > solEns[1].t[end]
        error("Error: fit range exceeds maximum simulation time.")
    end
    # `vaf` ∈ [0,0.5]
    vaf_t_Vid = Vector{Vector{Float64}}(undef, length(solEns)*ctrlParams[:nSamples])
    vCur = 1
    vNext = 1
    # ---- Sample trajectories ----
    for (i,sim) in enumerate(solEns)
        # sample trajectories from current sim
        vaf_t_VidSimCur = sampleSimTimepoints(
            sim, _t, ctrlParams[:nSamples], ctrlParams[:coverage];
            fMin=2*ctrlParams[:fMin], fMax=2*ctrlParams[:fMax]
            )./2
        vNext += length(vaf_t_VidSimCur)
        vaf_t_Vid[vCur:vNext-1] .= vaf_t_VidSimCur
        vCur = vNext
    end
    # ---- fit growth functions to trajectories ------ 
    models=["positive logistic", "negative logistic", "constant"]
    s_vid = DataStructuring.fitModelsLightweight((@view vaf_t_Vid[1:vCur-1]), _t, fill(ctrlParams[:coverage],length(_t)), ctrlParams[:params][:N], models; errorQuantile=ctrlParams[:sQError]
    )
    # ---- construct fitness distribution ----
    __, nVars_s = AnalysisTools.distbin(collect(skipmissing(s_vid)); bins=ctrlParams[:sBins], xMin=ctrlParams[:sFitBounds][1], xMax=ctrlParams[:sFitBounds][2], normalized=true)
    return nVars_s
end

"""
    runModelSim(paramsABC)

Perform a run of the model simulations to obtain a single particle with parameter set `paramsABC`.
"""
function runModelSim(paramsABC, ctrlParams)
    modelParams = deepcopy(ctrlParams[:params])
    # for pName in keys(paramsABC)
    for (pName, pVal) in pairs(paramsABC)
        ∈(pName, ctrlParams[:fixPar]) && continue
        modelParams[pName] = pVal
    end

    # run model sims
    model = GammaSelectionModel(modelParams[:s], modelParams[:σ], modelParams[:q])
    solEns, simArgs = evolvePopSim(modelParams, model; runs=ctrlParams[:simRuns], noDiffusion=false)
    tMeasure = (ctrlParams[:tBounds][1]+ctrlParams[:tBounds][2])/2
    parentVid_vid_Sid = simArgs[!, :parentId_vid]
    ctrlParams[:params] = modelParams

    # ============ size distribution ============
    nVars_f = sizeDistSims(tMeasure, solEns, parentVid_vid_Sid, ctrlParams)

    # ============ trajectory fitting ============
    nVars_s = fitTrajectoriesFitness(tMeasure, solEns, ctrlParams)

    # ============ mean number of variants above threshold ============
    detectedVarsAv_t = variantsAboveThreshold(solEns, 0; nSamples=500)

    return nVars_f, nVars_s, detectedVarsAv_t
end

function compareDataVSim(dataMetrics, simResults, thresholds)
    
    # unpack arguments
    nSim_f, nSim_s, nVarsThresh = simResults
    nData_f, nData_s = dataMetrics
    fThreshold, sThreshold = thresholds

    # compare size dist
    fDist = Distances.evaluate(Distances.chisq_dist, nSim_f, nData_f)
    fDistAccept = fDist <= fThreshold ? true : false
    
    # compare fitness dist
    sDist = Distances.evaluate(Distances.euclidean, mean(nSim_s), mean(nData_s))
    sDistAccept = sDist <= sThreshold ? true : false

    return (fDistAccept && sDistAccept)
end

function stdDist(n_x, _x)
    @. √( sum(_x^2*n_x) - sum(_x*n_x)^2 )
end

function meanDist(n_x, _x)
    sum(_x.*n_x)
end

function compareDataVSim2(dataMetrics, simResults, thresholds, _s)
    
    # unpack arguments
    nSim_f, nSim_s, _ = simResults
    nData_f, nData_s = dataMetrics
    fThreshold, (sAvThreshold, sStdThreshold) = thresholds

    # compare size dist
    fDist = Distances.evaluate(Distances.chisq_dist, nSim_f, nData_f)
    fDistAccept = fDist <= fThreshold ? true : false
    
    # compare fitness dist
    sAvDist = Distances.evaluate(Distances.euclidean, meanDist(nSim_s, _s), meanDist(nData_s, _s))
    sStdDist = Distances.evaluate(Distances.euclidean, stdDist(nSim_s, _s), stdDist(nData_s, _s))
    sDistAccept = (sAvDist <= sAvThreshold) && (sStdDist <= sStdThreshold)

    return (fDistAccept && sDistAccept)
end

function compareDataVSimError(dataMetrics, simResults, _s)
    
    # unpack arguments
    nSim_f, nSim_s, nVarsThresh = simResults
    nData_f, nData_s = dataMetrics

    # compare size dist
    fDist = Distances.evaluate(Distances.chisq_dist, nSim_f, nData_f)
    
    # compare fitness dist
    sAvDist = Distances.evaluate(Distances.euclidean, meanDist(nSim_s, _s), meanDist(nData_s, _s))
    sStdDist = Distances.evaluate(Distances.euclidean, stdDist(nSim_s, _s), stdDist(nData_s, _s))

    return (fDist, sAvDist, sStdDist)
end

"""
    Checks whether number of variants is is in bounds.
"""
function checkConstraintsSingleThreshold((nSim_f, nSim_s, nVarsThreshold), (lBound, uBound))
    (nVarsThreshold > lBound) && (nVarsThreshold < uBound)
end

function checkConstraintsEarlyLateThreshold((nSim_f, nSim_s, nVarsThresholdEarly, nVarsThresholdLate), (lBoundEarly, uBoundEarly), (lBoundLate, uBoundLate))
    (nVarsThresholdEarly > lBoundEarly) && (nVarsThresholdEarly < uBoundEarly) && (nVarsThresholdLate > lBoundLate) && (nVarsThresholdLate < uBoundLate)
end

function checkConstraintsDetectableVariants((nSim_f, nSim_s, nVars_t), tEarly, (lBoundEarly, uBoundEarly), tLate, (lBoundLate, uBoundLate))
    nVarsEarly = nVars_t[1+tEarly] # this is hacky af
    nVarsLate = nVars_t[1+tLate]
    (nVarsEarly>=lBoundEarly) && (nVarsEarly<=uBoundEarly) &&
    (nVarsLate>=lBoundLate) && (nVarsLate<=uBoundLate)
end

