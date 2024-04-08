

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
    nVarsDetectAv = 0
    fTotDetectAv = 0
    for simId in eachindex(solEns)
        if !isnothing(nSamples)
            vafSampled_vid = (vaf->sampleIntFromFreq(vaf, nSamples)).(solEns[simId](t))
            detectableVariants_vid = vafSampled_vid .> f0
        else
            detectableVariants_vid = solEns[simId](t) .> f0
        end
        nVarsDetectAv += sum(detectableVariants_vid) / length(solEns)
        fTotDetectAv += sum(vafSampled_vid[detectableVariants_vid]) / length(solEns)
    end
    return nVarsDetectAv, fTotDetectAv
end

"""
    Sample variant frequencies of the simulation result at timepoints `_t` with sample sizes `S_t`. Discard variants under frequency `f0`. Return a list of the total variant fractions `fTotDetectAv_t` for each timepoint and the number of variants detected `nVarsDetectAv_t`.
    `sol`     = simulation result
    `S_t`     = sample size for each timepoint
    `_t`      = list of timepoints
    `f0`      = minimum sampled frequency to accept
    returns `nVarsDetectAv_t`, `fTotDetectAv_t`
"""
function variantsAboveThreshold(solEns::EnsembleSolution, _t::AbstractVector{<:Real}, f0::Real, S_t::AbstractVector{Int})
    if length(S_t) !== length(_t)
        println("error: `S_t` should be same length as `_t`")
    end
    nVarsDetectAv_t = Vector{Float64}(undef, length(_t))
    fTotDetectAv_t = Vector{Float64}(undef, length(_t))
    for (tid, t) in enumerate(_t)
        nVarsDetectAv_t[tid], fTotDetectAv_t[tid] = variantsAboveThreshold(solEns, t, f0; nSamples=S_t[tid])
    end
    return nVarsDetectAv_t, fTotDetectAv_t
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

function drawVariantsRandomly(nVariants::Int, nVarsMax::Union{Nothing,Int})
    S = isnothing(nVarsMax) ? nVariants : minimum((nVarsMax, nVariants)) # if the number of variants is less than `nVarsMax`, use that
    vid_vidS = randperm(nVariants)[1:S]
    return vid_vidS
end

function chooseExclusiveRand(v_i, compare_j)
    vSelect = v_i[rand(1:end)]
    if !in(vSelect, compare_j)
        return vSelect
    else 
        return chooseExclusiveRand(v_i, compare_j)
    end
end
function chooseExclusiveNext(i, v_i, compare_j)
    vSelect = v_i[i]
    if !in(vSelect, compare_j)
        return vSelect
    else
        if i+1>length(v_i) return vSelect end
        return chooseExclusiveNext(i+1, v_i, compare_j)
    end
end

function drawVariantsWithPriority(
        s_vid::AbstractVector{Float64},
        nVarsMax::Union{Nothing,Int}=nothing,
        sortBias::Real=1.,
    )
    vid_vFit = (1:length(s_vid))[sortperm(s_vid, rev=true)]
    # vid_vRan = randperm(length(s_vid))
    S = 
        if isnothing(nVarsMax)
            length(s_vid)
        else
            minimum((nVarsMax, length(s_vid))) # if the number of existing variants is less than `nVarsMax`, use that
        end
    # vid_vidS = [rand()<=sortBias ? vid_vFit[i] : vid_vRan[i] for i in 1:S]
    vid_vidS = Vector{Int}(undef, S)
    for i in 1:S
        if rand()<=sortBias
            # select next fittest variant
            vid_vidS[i] = chooseExclusiveNext(i, vid_vFit, vid_vidS[1:i-1])
        else
            # select next random variant
            vid_vidS[i] = chooseExclusiveRand(vid_vFit, vid_vidS[1:i-1])
        end
    end

    return vid_vidS # variant ID's of sampled variants
end

function sampleSimTrajectories(sol::RODESolution, t::Real;
    tStep=1.,
    nTimeSteps::Int=4,
    freqCutoff=0.01,
    nVarsMax::Union{Nothing,Int}=nothing,
    s_vid::Union{Vector{Float64},Nothing}=nothing, # if passed, variant innate fitness will be used to determine sampling priority
    sortBias::Real=0, # determines the tendency of the sampler to bias towards the fittest mutants,
    verbose::Bool=false
    )
    obs_vid = sol(t) .>= freqCutoff
    if verbose
        println("total variants: ", length(obs_vid))
        println("observable variants:", sum(obs_vid))
    end
    vidO_vidS = 
        if isnothing(s_vid)
            drawVariantsRandomly(sum(obs_vid), nVarsMax)
        else
            verbose && println("drawing variants with priority bias=$(sortBias)...")
            drawVariantsWithPriority(s_vid[obs_vid], nVarsMax, sortBias)
        end
    _t = range(t; length=nTimeSteps, step=tStep)
    x_t_vidS = Array{Float64,2}(undef, nTimeSteps, length(vidO_vidS))
    for (i,tt) in enumerate(_t)
        for (j,x) in enumerate(sol(tt)[obs_vid][vidO_vidS])
            x_t_vidS[i,j] = x
        end
    end
    return _t, x_t_vidS
end

struct StepUniform <: Sampleable{Univariate,Continuous}
    edges::Vector{Float64}
    counts::Vector{Int64}
end

function Base.rand(rng::AbstractRNG, d::StepUniform)
    bin = StatsBase.sample(rng, range(1,length(d.counts)), Weights(d.counts))
    sample = d.edges[bin] + (d.edges[bin+1]-d.edges[bin])*rand()
    return sample
end

function fitSamplesGrowth(solEns::Union{EnsembleSolution,Vector{T}} where T, params::Dict;
        tMeasure::Union{Real,Tuple{Real, Real},Tuple{AbstractArray,AbstractArray}}=(50,80), # time measurement range
        timeLimit=1.0,
        nTimeSteps=4,
        tStep=1.,
        errorQuantile=0.99,
        freqCutoff=0.01,
        nVarsMax=4,
        s_vid_Sim::Union{Nothing,Vector{Vector{Float64}}}=nothing, # if passed, variant innate fitness will be used to determine sampling priority
        sortBias::Real=0., # determines the tendency of the sampler to bias towards the fittest mutants
        verbose::Bool=false,
    )
    Nf = params[:N]
    dfVid = DataFrame(
        _t=AbstractVector{Real}[],
        vaf_t=Vector{Float64}[],
        cov_t=Vector{Int64}[]
    )

    # first measurement time of every patient
    t_pid = let
        if length(tMeasure)<2
            tMeasure
        elseif eltype(tMeasure) <: Real #in this case draw from uniform with boundaries given by `tMeasure`
            rand(length(solEns)).*(tMeasure[2]-tMeasure[1]).+tMeasure[1]
        else
            dist = StepUniform(tMeasure[1],tMeasure[2])
            rand(dist, length(solEns))
        end
    end
    if verbose
        println("t_pid: ")
        display(t_pid)
    end
    # sample simulation and fit variants
    for (j,sol) in enumerate(solEns)
        t = t_pid[j] # first measurement time point
        _t, x_t_cid = sampleSimTrajectories(
            sol, t;
            tStep,
            nTimeSteps,
            freqCutoff,
            nVarsMax,
            s_vid = !isnothing(s_vid_Sim) ? s_vid_Sim[j] : nothing,
            sortBias,
            verbose,
        )
        verbose && println("detected variants: ", size(x_t_cid,2))
        for cid in axes(x_t_cid, 2)
            vaf_t = x_t_cid[:, cid]/2
            # make sure vaf>0.5 and vaf<0 cases are taken care of
            for (tid,vaf) in enumerate(vaf_t)
                if vaf>=0.5
                    vaf_t[tid] = 0.49999
                elseif vaf<=0
                    vaf_t[tid] = 0.000001
                end
            end
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
function sampleSimTimepoints(sim::RODESolution, _t::AbstractVector, nSamples::Int, coverage::Int; 
        fMin::Float64=0.,
        fMax::Float64=1.,
        s_vid::Union{Nothing, AbstractVector{Float64}}=nothing,
        sortBias::Real=0.,
    )
    # `x`∈[0,1]; not to be confused with `vaf`∈[0,0.5]
    x_vid_T = sim(_t)
    # select only observable variants
    obs_vid = fMin .<= x_vid_T[1] .<= fMax
    vidO_vidS = 
        if !isnothing(s_vid)
            # we take 4*nSamples to allow for the possibility of variants being lost later during binomial sampling
            drawVariantsWithPriority(s_vid[obs_vid], 4*nSamples, sortBias)
        else
            drawVariantsRandomly(length(s_vid[obs_vid]), 4*nSamples)
        end
    j = 0   # index for accepted observed variants
    x_t_VidSampled = Vector{Vector{Float64}}(undef, nSamples)
    for vidO in vidO_vidS
        x_t = [x_vid[obs_vid][vidO] for x_vid in x_vid_T]
        xObs_t = (x -> sampleFreqFromFreq(x, coverage)).(x_t)
        #check whether variant has enough nonzero observations; if not skip to next variant (this is why _vidS has length 4*nSamples)
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
function sampleSimTimepoint(
        sim::RODESolution,
        vid_child_Vid::Vector{Vector{T}} where T<:Integer,
        t::Real,
        nSamples::Int,
        coverage::Int;
        fMin::Float64=0.,
        fMax::Float64=1.,
        s_vid::Union{Nothing, AbstractVector{Float64}}=nothing,
        sortBias::Real=0.,
    )
    obs_vid = fMin .<= sim(t) .<= fMax
    vidO_vidS = 
        if isnothing(s_vid)
            drawVariantsRandomly(length(s_vid[obs_vid]), 4*nSamples)
        else
            drawVariantsWithPriority(s_vid[obs_vid], 4*nSamples, sortBias)
        end
    vafObs_j = Vector{Float64}(undef, nSamples)
    j = 0
    for vidO in vidO_vidS
        x = @view(sim(t)[obs_vid])[vidO]
        # add size of children (if they exist)
        for vid in vid_child_Vid[obs_vid][vidO]
            x += sim(t)[vid]
        end
        (x<fMin || x>fMax) && continue
        # allele frequency of a variant is f=x/2
        vafObs = rand(Binomial(coverage, x/2)) / coverage
        vafObs==0 && continue
        j+=1
        vafObs_j[j]=vafObs
        j==nSamples && break
    end
    if j<nSamples
        return vafObs_j[1:j]
    end
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

function sizeDistSims(
    tMeasure,
    solEns,
    parentVid_vid_Sid::Vector{Vector{T}} where T<:Integer,
    s_vid_Sid::Vector{Vector{Float64}},
    ctrlParams;
    sortBias=0.,
    )
    # -------- sample sims --------
    vaf_vidS = Vector{Float64}(undef, ctrlParams[:simRuns]*ctrlParams[:nSamples])
    vCur = 1
    vNext = 1
    for (sid,sol) in enumerate(solEns)
        childVid_child_Vid = buildFamilyArray(parentVid_vid_Sid[sid])
        f_vidSimSample = sampleSimTimepoint(
                sol,
                childVid_child_Vid,
                tMeasure,
                ctrlParams[:nSamples],
                ctrlParams[:coverage];
                fMin = 2*ctrlParams[:fMin],
                fMax = 2*ctrlParams[:fMax],
                s_vid = s_vid_Sid[sid],
                sortBias = sortBias,
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

function fitTrajectoriesFitness(
        tMeasure,
        solEns,
        s_vid_Sid::Vector{Vector{Float64}},
        ctrlParams;
        sortBias=0.
    )
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
            fMin=2*ctrlParams[:fMin],
            fMax=2*ctrlParams[:fMax],
            # s_vid=simArgs[i,:s_vid],
            s_vid=s_vid_Sid[i],
            sortBias,
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
    Sample variant frequencies of the simulation result at timepoints `_t` with sample sizes `S_t`. Discard variants under frequency `f0`. Return a list `s_vid_T` with for each timepoint a list of the sampled variant frequencies above f0.
    `sol`     = simulation result
    `S_t`     = sample size for each timepoint
    `_t`      = list of timepoints
    `f0`      = minimum sampled frequency to accept
    returns `s_vid_T`
"""
function samplePatientSim(sol, S_t::Vector{Int}, _t; f0=0)
    nVars = length(sol[1])
    s_vid_T = Vector{Vector{Float64}}(undef, length(_t))
    for (tInd,t) in enumerate(_t)
        s_vid_T[tInd] = Float64[]
        for vid in 1:nVars
            # draw S times with success prob fVid to get frequency of vid sVid in sample
            binomDist = Binomial(S_t[tInd], sol(t)[vid])
            sVid = rand(binomDist) / S_t[tInd]
            if sVid>f0 push!(s_vid_T[tInd], sVid) end
        end
    end
    return s_vid_T
end

function samplePatientSim(sol, S::Int, _t; f0=0)
    S_t = fill(S, length(_t))
    samplePatientSim(sol, S_t, _t; f0)
end

function runModelSimFixedFitness(paramsABC, ctrlParams)
    # fix `s` based on value of τ
    ctrlParams[:params][:s] = ctrlParams[:sFixed] * paramsABC[:τ]
    runModelSim(paramsABC, ctrlParams)
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

    if ctrlParams[:normalizedFitnessDist]
        modelParams[:s] = modelParams[:s]*modelParams[:τ]
        modelParams[:σ] = modelParams[:σ]*modelParams[:τ]
    end

    # run model sims
    model = GammaSelectionModel(modelParams[:s], modelParams[:σ], modelParams[:q])
    solEns, simArgs = evolvePopSim(modelParams, model; runs=ctrlParams[:simRuns], noDiffusion=false)
    tMeasure = (ctrlParams[:tBounds][1]+ctrlParams[:tBounds][2])/2
    parentVid_vid_Sid = simArgs[!, :parentId_vid]
    ctrlParams[:params] = modelParams

    # ============ size distribution ============
    sortBias =
        if haskey(ctrlParams, :sortBias)
            ctrlParams[:sortBias]
        else
            modelParams[:sortBias]
        end
    nVars_f = sizeDistSims(
        tMeasure,
        solEns,
        parentVid_vid_Sid,
        simArgs[!,:s_vid],
        ctrlParams;
        sortBias=sortBias,
    )

    # ============ logistic fitness fitting ============
    nVars_s = fitTrajectoriesFitness(
        tMeasure,
        solEns,
        simArgs[!,:s_vid],
        ctrlParams;
        sortBias=sortBias
    )

    # ============ mean number of variants and variant sizes above threshold ============
    f0 = ctrlParams[:fSampleThreshold]
    _tS = ctrlParams[:_tS]
    S_tS = ctrlParams[:S_tS]
    nVarsDetectAv_t, fTotDetectAv_t = variantsAboveThreshold(solEns, _tS, f0, S_tS)


    return nVars_f, nVars_s, nVarsDetectAv_t, fTotDetectAv_t
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
    nSim_f, nSim_s, _, _ = simResults
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

function distDataVSim(dataMetrics, simResults, _s)
    
    # unpack arguments
    nSim_f, nSim_s, _, _ = simResults
    nData_f, nData_s = dataMetrics  

    # compare size dist
    fDist = Distances.evaluate(Distances.chisq_dist, nSim_f, nData_f)
    
    # compare fitness dist
    sAvDist = Distances.evaluate(Distances.euclidean, meanDist(nSim_s, _s), meanDist(nData_s, _s))
    sStdDist = Distances.evaluate(Distances.euclidean, stdDist(nSim_s, _s), stdDist(nData_s, _s))

    return fDist, sAvDist, sStdDist
end

function compareVariantNumbers(nData_t, simResults, threshold)
    # unpack arguments
    _, _, nSim_t, _ = simResults

    # compare data and sim number of variants
    rsqerror = .√((nData_t .- nSim_t ./ ((n->n==0 ? 1 : n).(nData_t))).^2) |> sum
    return rsqerror <= threshold
end

function compareVariantSizes(xData_t, simResults, threshold)
    # unpack arguments
    _, _, _, xSim_t = simResults

    # compare data and sim variant sizes
    error = .√((xData_t .- xSim_t ./ ((n->n==0 ? minimum(xData_t[xData_t.>0]) : n).(xData_t))).^2) |> sum
    return error <= threshold
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


## ----------------------------------------
#region - Plotting AnalysisTools

function simResultArray(sol::RODESolution)
    _t = sol.t
    size_t_vid = Array{Float64,2}(undef, length(_t), length(sol[1]))
    for vid in eachindex(sol[1])
        for tid in eachindex(sol.u)
            size_t_vid[tid,vid] = sol.u[tid][vid]
        end
    end
    return _t, size_t_vid
end

function mullerPlotArgs(sol::RODESolution, detectThreshold)
    _t, size_t_vid = simResultArray(sol)
    visible_vid = [any(size_t_vid[:,vid] .> detectThreshold) for vid in eachindex(sol[1])]
    visible_vid0 = prepend!(copy(visible_vid), [true,])
    n_t_vid0 = Array{Float64,2}(undef, length(_t), length(sol[1])+1)
    for i in 1:length(sol[1])
        n_t_vid0[:,1+i] = [solT[i] for solT in sol.u]
    end
    n_t_vid0[:,1] = [1-sum(n_t_vid0[t,2:end]) for t in 1:size(n_t_vid0,1)]
    n_t_vid0[findall(x->x<0, n_t_vid0)] .= 0
    return n_t_vid0, visible_vid0
end

#endregion