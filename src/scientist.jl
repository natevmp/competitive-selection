

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

# """
# """
# function variantsPopulationCoverage(solEns::EnsembleSolution; threshold=0)
#     _t = solEns[1].t
#     for (tInd, t) in enumerate(_t)
#         solEns[i]

# end

function variantsAboveThreshold(solEns::EnsembleSolution, f0=0.02)
    nVars_sim_t = Array{Float64}(undef, runs, length(solEns[1]))
    for simId in 1:runs
        for tInd in 1:paramsPop[:T]
            detectableVariants_vid = solEns[simId][tInd] .> f0
            nVars_sim_t[simId,tInd] = sum(detectableVariants_vid)
        end
    end
    return vec(mean(nVars_sim_t,dims=1))
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

function fitSamplesGrowth(solEns::Union{EnsembleSolution,Vector{T}} where T, params::Dict; tMeasure::Union{Real,Tuple{Real, Real}}=(50,80), timeLimit=1.0)
    Nf = params[:N]
    tStep = 1.
    nTimeSteps = 4
    # fit_tλ_cid = ElasticArray{Float64}(undef,2,0)
    dfVid = DataFrame(
        _t=AbstractVector{Real}[],
        vaf_t=Vector{Float64}[],
        cov_t=Vector{Int64}[]
    )
    for sim in solEns
        t = let
            if length(tMeasure)<2
                tMeasure
            else
                tMeasure[1]+rand()*(tMeasure[2]-tMeasure[1])
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
    DataStructuring.fitVariantsData!(dfVid, Nf)
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


