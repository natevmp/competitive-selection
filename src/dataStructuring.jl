module DataStructuring

using CSV, DataFrames, Statistics, LsqFit, StatsBase, Distributions, Optim, Optimization, Ipopt, ProgressMeter
using OptimizationBBO, OptimizationMOI
include("./analysisTools.jl")
using .AnalysisTools

export analyseVariants, meanVariantParamTime, dropdupes, binVariantParamTime

include("fitModels.jl")

"""
    fixArrTypeToType(arr, T)

Create an array `newArr` with eltype `T` from parsing values in array `arr`.
"""
function fixArrTypeToType(arr, T::DataType)
    eltype(arr) == T && return arr
    newArr = Vector{T}(undef, length(arr))
    for (i,x) in enumerate(arr)
        try 
            newArr[i] = parse(T,x)
        catch
            newArr[i] = 0.
        end
    end
    return newArr
end

"""
    structureData!(df)

Alter dataframe `df` with variant call information in the following ways:
- Remove rows where key `:Gene` has value "NA".
- Change column datatypes where numerical data is specified in strings to relevant numerical type (Int64 or Float64).
"""
function structureData!(df::DataFrame)
    filter!(:Gene => g -> g .!= "NA", df)
    df[!,:VAF] .= fixArrTypeToType(df[!,:VAF], Float64)
    df[!,:TOTALcount] .= fixArrTypeToType(df[!,:TOTALcount], Int64)
    df[!,:WTcount] .= fixArrTypeToType(df[!,:WTcount], Int64)
    df[!,:MUTcount_Xadj] .= fixArrTypeToType(df[!,:MUTcount_Xadj], Float64)
    return nothing
end

"""
    distbin(x_vid)

Bin elements of `x_vid` withing range `xMin` to `xMax`.
"""
function distbin(x_vid::Vector{T} where T<:Real; bins::Int=25, xMin::Real=0, xMax::Real=1., normalized::Bool=true)
    _xEdges = range(xMin, xMax, length=bins+1)
    dx = _xEdges[2]-_xEdges[1]
    _x = _xEdges[1:end-1] .+ dx/2
    withinRangeBool_vid = (x -> x>=xMin && x<xMax).(x_vid)   #apply thresholds
    nV_x = zeros(Float64, bins)
    for x in x_vid[withinRangeBool_vid]
        k = ((x-xMin)/dx |> floor |> Int) + 1 # get index of bin containing `x`
        nV_x[k] += 1
    end
    normalized && (nV_x ./= (length(x_vid)*dx))
    return _x, nV_x
end

"""
    dropdupes(arr)

Create an array with only the unique elements in `arr`.
"""
function dropdupes(arr::Vector{<:Real})
    newArr = eltype(arr)[]
    for el in arr
        el in newArr ? continue : push!(newArr, el)
    end
    return newArr
end
function dropdupes(arr::AbstractVector)
    newArr = eltype(arr)[]
    for el in arr
        el in newArr ? continue : push!(newArr, el)
    end
    return newArr
end

"""
    predict(_t, shape)

Compute the predicted curve for a model with form and parameters defined by `shape`.
"""
function predict(_t, shape::Union{LogisticModel, MaxAdjustedLogisticModel})
    [fLogistic(t, shape.t0, shape.γ, shape.x0, shape.xF) for t in _t]
end

function predict(_t, shape::ConstantModel)
    [fConstant(t, shape.x) for t in _t]
end

"""
    fitLogGrowthMLE(_t, vaf_t, nCov_t, Nf)

Compute the best fit parameter values of `t0`, `γ`, and `x0` in the logistic growth model
    `fLogistic(t, t0, r, x0) = 0.5 / ( 1 + (0.5-x0)/x0 * exp(-γ*(t-t0)) )`
using a Maximum likelihood estimation whereby sampling is assumed to have been performed binomially.
"""
function fitLogGrowthMLE!(_t, vaf_t, nCov_t, Nf, modelShape::ModelShape; timeLimit=1.0)

    X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on
    expansionType =
        if corspearman(_t, vaf_t)>=0
            PositiveExpansion()
        else
            NegativeExpansion()
        end
    model = ModelFit(modelShape, expansionType)
    initial_β = setInitialParamsMLE!(model, _t, vaf_t, Nf)

    function logprobX(X, t, n, β)
        p = fLogModel(t,β,model)
        if p <= 0 || p > 1 || isnan(p)
            #! debug
            #region
            # display(model)
            # println("\nX: ", X)
            # println("β: ", β)
            # println("p: ", p)
            # println("t: ", t)
            #endregion
            p = 1E-20
        end
        logpdf(Binomial(n, p), X)
    end
    logllh(β, p) = 
        -([logprobX(X_t[i], _t[i], nCov_t[i], β) for i in eachindex(_t)] |> sum)
    #! debug
    #region
    # function logllh(β)
    #     try
    #         -([
    #             logprobX(X_t[i], _t[i], nCov_t[i], β) for i in eachindex(_t)
    #         ] |> sum)
    #     catch
    #         println("X_t: ", X_t)
    #     end
    # end
    #endregion

    sol = 
        if typeof(model)==ModelFit{MaxAdjustedLogisticModel, NegativeExpansion}
            # maxAdjLog, NE
            constraintF, lcons, ucons = constraintsMLE(model, _t)
            optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff(), cons=constraintF)
            prob = Optimization.OptimizationProblem(optprob, initial_β, lcons=lcons, ucons=ucons)
            solve(prob, Ipopt.Optimizer(); print_level=0)
        else
            # other
            lcons, ucons = boundsMLE(model, _t, vaf_t)
            optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff())
            prob = Optimization.OptimizationProblem(optprob, initial_β, lb=lcons, ub=ucons)
            solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
            # solve(prob, LBFGS())
        end
    
    addFitParams!(model, sol.u)
    return true

    # alg = Fminbox(BFGS())
    # lower_β, upper_β = boundsMLE(model, _t, vaf_t)
    # fit = optimize(logllh, lower_β, upper_β, initial_β, alg, Optim.Options(time_limit = timeLimit))
    # #! debug
    # #region
    # # fit = 
    # #     try
    # #         optimize(logllh, lower_β, upper_β, initial_β, alg, Optim.Options(time_limit = timeLimit))
    # #     catch e
    # #         println("lower_β: ", lower_β)
    # #         println("upper_β: ", upper_β)
    # #         throw(e)
    # #     end
    # #endregion
    # _β = Optim.minimizer(fit)
    # addFitParams!(model, _β)
    # return Optim.converged(fit)
end

function fLogCompModel(t, Nf, b)
    # b[1]=t0, b[2]=s, b[3]=a, b[4]=μ
    AnalysisTools.fIndFitnessWithAvCompetition(t, b[1], b[2], b[3], b[4], Nf)
end

"""
    fitLogCompGrowthMLE(_t, vaf_t, nCov_t, Nf)

Compute the best fit parameter values of `t0`, `s`, `a`, and `μ` in the competitive logistic growth model using a Maximum likelihood estimation whereby sampling is assumed to have been performed binomially.
"""
function fitLogCompGrowthMLE(_t, vaf_t, nCov_t, Nf; timeLimit=1.0)

    X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on

    function logprobX(X, t, n, b)
        # length(b)>3 && b[4]<0 && (b[4]=0.0001)
        p = fLogCompModel(t,Nf,b)
        #! debug
        # if p < 0 || p > 1 || isnan(p)
        #     println("\nb: ", b)
        #     println("p: ", p)
        #     println("t: ", t)
        # end
        #! debug
        p > 1 && (p=0.99)
        p < 0 && (p=0.001)
        logpdf(Binomial(n, p), X)
    end

    logllh(b) = -([
        logprobX(X_t[i], _t[i], nCov_t[i], b) for i in eachindex(_t)
    ] |> sum)

    #! debug
    # println("X_t: ", X_t)
    # for i in eachindex(_t) println("fLogistic: ", fLogModel(_t[i], Nf, logitModel.b0, logitModel)) end
    # for i in eachindex(_t) println("logprobX: ", logprobX(X_t[i], _t[i], nCov_t[i], logitModel.b0) ) end
    # for i in eachindex(_t) println("coverage: ", nCov_t[i] ) end
    # println("llh: ", logllh(logitModel.b0))
    #! debug

    # b[1]=t0, b[2]=s, b[3]=a, b[4]=μ
    lower = [0, 1E-4, 1E-4, 1E-2]
    upper = [_t[1], 10., 10., 10.]
    b0 = [5., 0.2, 0.2, 1.]
    # lower = [0, 1E-4, 0.5]
    # upper = [_t[1], 5., 5.]
    # b0 = [5., 0.1, 0.6]

    alg = Fminbox(BFGS())
    return optimize(logllh, lower, upper, b0, alg, Optim.Options(time_limit = timeLimit))
end

"""
    sortDataPerVariant(df)

Create a new DataFrame `dfVid` from the original Fabre data, whereby the rows correpond to the distinct variants observed. I.e. each row describes a single variant trajectory.
"""
function sortDataPerVariant(df::DataFrame)
    _pid = dropdupes(df[!,:SardID])
    dfVid = DataFrame(
        _t=Vector{Float64}[],
        vaf_t=Vector{Float64}[],
        pid=Int64[],
        gene=String7[],
        cov_t=Vector{Int64}[], # coverage
        n_t=Vector{Float64}[], # mutation counts
        vid=Int16[], # variant id number
    )
    vidCur::Int16 = 0
    for pid in _pid
        _vid = 
        dropdupes(
            subset(df, :SardID => s->s.==pid)[!,"AAChange.refGene"]
        )
        for (i,vid) in enumerate(_vid)
            subDf = subset(
                df,
                "AAChange.refGene" => ref->ref.==vid,
                "SardID" => s->s.==pid,
            )
            
            
            length( findall((v->v>0).(subDf[:,:VAF]) )) < 2 && continue #skip variant if only 1 non-zero VAF measured
            t0Ind = findfirst(subDf[:,:VAF].>0) # get first timepoint where the variant is detected
            length(subDf[t0Ind:end,:Age]) < 3 && continue # skip variant if less than 3 datapoints
            vidCur += 1

            # _t = subDf[1:end,:Age]
            # vaf_t = subDf[1:end,:VAF]
            # cov_t = subDf[1:end, :TOTALcount] # coverage
            # n_t = subDf[1:end, :MUTcount_Xadj] # mutation counts
            _t = subDf[t0Ind:end,:Age]
            vaf_t = subDf[t0Ind:end,:VAF]
            cov_t = subDf[t0Ind:end, :TOTALcount] # coverage
            n_t = subDf[t0Ind:end, :MUTcount_Xadj] # mutation counts
            gene = subDf[1,:Gene]
            vid = vidCur
            push!(dfVid, (_t, vaf_t, pid, gene, cov_t, n_t, vid))
        end
    end
    return dfVid
end

"""
    fitVariantsData!(dfVid)

Compute the best fit parameters of the logistic model for each variant in the DataFrame `dfVid`, and add the results as new columns.
"""
function fitVariantsData!(dfVid::DataFrame, Nf::Real; competitionSaturate::Bool=false)
    nVars = size(dfVid,1)

    modelShape = begin
        !competitionSaturate ? LogisticModel() : MaxAdjustedLogisticModel()
    end
    fitVars = fieldnames(typeof(modelShape))
    for varName in fitVars
        dfVid[:, varName] = Vector{Float64}(undef, nVars)
    end
    dfVid.ρ = Vector{Float64}(undef, nVars)
    dfVid.fitConverge = Vector{Bool}(undef, nVars)
    dfVid._tt = Vector{
        StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    }(undef, nVars)
    dfVid.vafTh_tt = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafTh_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafThQl_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafThQu_t = Vector{Vector{Float64}}(undef, nVars)
    @showprogress for dfVCur in eachrow(dfVid)
        converged = fitLogGrowthMLE!(dfVCur[:_t], dfVCur[:vaf_t], dfVCur[:cov_t], Nf, modelShape)
        dfVCur[:ρ] = corspearman(dfVCur[:_t], dfVCur[:vaf_t])
        dfVCur[:fitConverge] = converged
        for param in fieldnames(typeof(modelShape))
            dfVCur[param] = getfield(modelShape, param)
        end
        dfVCur[:_tt] = range(dfVCur[:_t][1],dfVCur[:_t][end], length=25)
        dfVCur[:vafTh_tt] = predict(dfVCur[:_tt], modelShape)
        dfVCur[:vafTh_t] = predict(dfVCur[:_t], modelShape)
        dfVCur[:vafThQl_t] = [
            quantile(Binomial(dfVCur[:cov_t][tInd], dfVCur[:vafTh_t][tInd]), 0.01) / dfVCur[:cov_t][tInd]
            for tInd in eachindex(dfVCur[:_t])
        ]
        dfVCur[:vafThQu_t] = [
            quantile(Binomial(dfVCur[:cov_t][tInd], dfVCur[:vafTh_t][tInd]), 0.99) / dfVCur[:cov_t][tInd]
            for tInd in eachindex(dfVCur[:_t])
        ]
    end
end

"""
    fitVariantsCompData!(dfVid)

Compute the best fit parameters of the logistic competition model for each variant in the DataFrame `dfVid`, and add the results as new columns.
"""
function fitVariantsCompData!(dfVid::DataFrame, Nf::Real)
    # b[1]=t0, b[2]=s, b[3]=a, b[4]=μ
    nVars = size(dfVid,1)
    dfVid.t0 = Vector{Float64}(undef, nVars)
    dfVid.s = Vector{Float64}(undef, nVars)
    dfVid.a = Vector{Float64}(undef, nVars)
    dfVid.μ = Vector{Float64}(undef, nVars)
    dfVid._tt = Vector{
        StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    }(undef, nVars)
    dfVid.vafTh_tt = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafTh_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.fitConverge = Vector{Bool}(undef, nVars)
    dfVid.vafThQl_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafThQu_t = Vector{Vector{Float64}}(undef, nVars)
    @showprogress for dfVCur in eachrow(dfVid)
        fit = fitLogCompGrowthMLE(dfVCur[:_t], dfVCur[:vaf_t], dfVCur[:cov_t], Nf)
        dfVCur[:fitConverge] = Optim.converged(fit)
        _β = Optim.minimizer(fit)
        dfVCur[:t0] = _β[1]
        dfVCur[:s] = _β[2]
        dfVCur[:a] = _β[3]
        dfVCur[:μ] = _β[4]
        dfVCur[:_tt] = range(dfVCur[:_t][1],dfVCur[:_t][end], length=25)
        dfVCur[:vafTh_tt] = [
            fIndFitnessWithAvCompetition(t, dfVCur[:t0], dfVCur[:s], dfVCur[:a], dfVCur[:μ], Nf)
            for t in dfVCur[:_tt]
        ]
        dfVCur[:vafTh_t] = [
            fIndFitnessWithAvCompetition(t, dfVCur[:t0], dfVCur[:s], dfVCur[:a], dfVCur[:μ], Nf)
            for t in dfVCur[:_t]
        ]
        dfVCur[:vafThQl_t] = [
            quantile(Binomial(dfVCur[:cov_t][tInd], dfVCur[:vafTh_t][tInd]), 0.01) / dfVCur[:cov_t][tInd]
            for tInd in eachindex(dfVCur[:_t])
        ]
        dfVCur[:vafThQu_t] = [
            quantile(Binomial(dfVCur[:cov_t][tInd], dfVCur[:vafTh_t][tInd]), 0.99) / dfVCur[:cov_t][tInd]
            for tInd in eachindex(dfVCur[:_t])
        ]
    end
end

"""
    analyseVariants(df, Nf)

Create a new DataFrame `dfVid` with per variant rows from the Fabre data `df`, and compute and add logistic fit information for each variant.
"""
function analyseVariants(df::DataFrame, Nf::Real; competitionSaturate::Bool=false, competition::Bool=false)
    structureData!(df)
    dfVid = sortDataPerVariant(df)
    if !competition
        fitVariantsData!(dfVid, Nf; competitionSaturate)
    else
        fitVariantsCompData!(dfVid, Nf)
    end
    return dfVid
end


"""
    binVariantParamTime(age_vid, p_vid; bins::Int=20) -> _tBinEdge, p_vid_TBin

Bin variant estimated logistic growth rates `λ_vid` measured at times `age_vid` in time-spaced bins. 
"""
function binVariantParamTime(age_vid, p_vid; bins::Int=20)
    _tBinEdge = range(minimum(age_vid),maximum(age_vid)+1, length=bins+1)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    p_vid_TBin = [eltype(p_vid)[] for i in eachindex(_tBin)]
    for (vid,t) in enumerate(age_vid)
        tBin = 1 + floor((t-_tBinEdge[1])/Float64(_tBinEdge.step)) |> Int
        push!(p_vid_TBin[tBin], p_vid[vid])
    end
    return _tBinEdge, p_vid_TBin
end

"""
    binVariantGRTime(dfVid::DataFrame)

Bin variant estimated logistic growth rates `λ_vid` measured at times `age_vid` in time-spaced bins. 
"""
function binVariantGRTime(dfVid::DataFrame; bins::Int=20)
    age_vid = [dfVidCur[:_t][1] for dfVidCur in eachrow(dfVid[!,:])]
    λ_vid = dfVid[!, :λ]
    binVariantParamTime(age_vid, λ_vid; bins)
end

"""
    meanVariantParamTime(age_vid, λ_vid)

Compute the mean of a fitted parameter at binned timespans from variant measure times `age_vid`.
Return the bin midpoints `_tBin` and mean growth rates per bin `vafAv_tBin`.
"""
function meanVariantParamTime(age_vid, p_vid; bins::Int=20)
    _tBinEdge, p_vid_TBin = binVariantParamTime(age_vid, p_vid; bins=bins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    pAv_tBin = Vector{Float64}(undef, length(_tBin))
    for i in eachindex(p_vid_TBin)
        pAv_tBin[i] = mean(p_vid_TBin[i])
    end
    return _tBin, pAv_tBin
end



"""
    sizeDistBinTime(age_vid, v_vid; bins::Int=21)

Compute the size distribution of `v_vid` at binned timespans from measurement times `age_vid`.
"""
function sizeDistBinTime(age_vid, f_vid; tBins::Int=20, fBins::Int=25)
    _tBinEdge, f_vid_TBin = binVariantParamTime(age_vid, f_vid; bins=tBins)
    _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
    _fBinEdge = range(0, 0.5, length=fBins+1)
    _fBin = _fBinEdge[1:end-1] .+ Float64(_fBinEdge.step)/2
    densV_f_TBin = Vector{Vector{Float64}}(undef, tBins)
    for (tInd,f_vid) in enumerate(f_vid_TBin)
        # create distribution for every timepoint
        nVars_f = zeros(Int, fBins)
        for f in f_vid
            binInd = 1 + ((f-_fBinEdge[1])/Float64(_fBinEdge.step) |> floor) |> Int
            nVars_f[binInd] += 1
        end
        nTot = sum(nVars_f)
        if nTot==0
            densV_f_TBin[tInd] = nVars_f
            continue
        end
        densV_f_TBin[tInd] = (n -> n*_fBinEdge.step/nTot).(nVars_f)
    end
    return _tBin, _fBin, densV_f_TBin
end

"""
Compute the size distribution of `f_vid` in absolute number of elements per bin.
"""
function sizeDistAbs(f_vid; fBins::Int=25, fMin=0., fMax=0.5)
    _fBinEdge = range(fMin, fMax, length=fBins+1)
    n_f = zeros(Int64, fBins)
    for f in f_vid
        binInd = 1 + ((f-fMin)/Float64(_fBinEdge.step) |> floor) |> Int
        n_f[binInd] += 1
    end
    return _fBinEdge, n_f
end

"""
Compute the size distribution of `f_vid` as a density per bin.
"""
function sizeDistDens(f_vid; fBins, fMin, fMax)
    _fBinEdge, n_f = sizeDistAbs(f_vid; fBins, fMin, fMax)
    nDens_f = (n -> n*Float64(_fBinEdge.step)/sum(n_f)).(n_f)
    # _fBinCentre = _fBinEdge[1:end-1] .+ Float64(_fBinEdge.step)/2
    return _fBinEdge, nDens_f
end

function sizeDistDens(f_vid, t_vid, tBounds; fBins::Int=25, fMin=0., fMax=0.5)
    tMask_vid = ( t -> t>tBounds[1] && t<tBounds[2] ).(t_vid)
    _fBinEdge, nDens_f = sizeDistDens(@view f_vid[tMask_vid]; fBins, fMin, fMax)
    return _fBinEdge, nDens_f
end

# """
# mapTimeBinned(age_vid, v_vid, fmap; bins::Int=21)

# Perform the function `fmap` on `v_vid` at binned timespans from variant measure times `age_vid`.
# """
# function mapTimeBinned(age_vid, v_vid, fmap; bins::Int=21)
#     _tBinEdge, λ_vid_TBin = binVariantsTime(age_vid, v_vid; bins=bins)
#     _tBin = _tBinEdge[1:end-1] .+ Float64(_tBinEdge.step)/2
#     λAv_tBin = Vector{Float64}(undef, length(_tBin))
#     for i in eachindex(λ_vid_TBin)
#         f_tBin[i] = mean(λ_vid_TBin[i])
#     end
#     return _tBin, λAv_tBin
# end


## =================================================================

function fitAllModels!(dfVid::DataFrame, Nf::Real, models::Vector{String}=["positive logistic", "negative logistic", "constant"]; errorQuantile=0.99)
    nVars = size(dfVid,1)

    dfVid.fitType = Vector{String}(undef, nVars)
    dfVid.ρ = Vector{Float64}(undef, nVars)
    dfVid._tt = Vector{
        StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    }(undef, nVars)
    dfVid.vafTh_tt = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafTh_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafThQl_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.vafThQu_t = Vector{Vector{Float64}}(undef, nVars)
    dfVid.fit_t = Vector{Vector{Bool}}(undef, nVars)
    dfVid.goodFit = Vector{Bool}(undef, nVars)

    modelVars = Symbol[]
    for name in models
        append!(modelVars, fieldnames(modelShapeFromName(name)))
    end
    # varNames = [:γ, :t0, :x0, :xF, :x]
    varNames = dropdupes(modelVars)
    for varName in varNames
        dfVid[:, varName] = Vector{Union{Missing,Float64}}(missing, nVars)
    end

    @showprogress for dfVCur in eachrow(dfVid)
        _models = [createModelFromName(model) for model in models]
        llh_model = [fitGrowthModel!(dfVCur[:_t], dfVCur[:vaf_t], dfVCur[:cov_t], Nf, model)
            for model in _models]
        llhMax, indMax = findmax(llh_model)
        modelFit = _models[indMax]
        dfVCur[:fitType] = modelName(modelFit)

        dfVCur[:ρ] = corspearman(dfVCur[:_t], dfVCur[:vaf_t])
        for param in fieldnames(typeof(modelFit.shape))
            dfVCur[param] = getfield(modelFit.shape, param)
        end

        dfVCur[:_tt] = range(dfVCur[:_t][1],dfVCur[:_t][end], length=25)
        dfVCur[:vafTh_tt] = predict(dfVCur[:_tt], modelFit.shape)
        dfVCur[:vafTh_t] = predict(dfVCur[:_t], modelFit.shape)
        dfVCur[:vafThQl_t], dfVCur[:vafThQu_t], dfVCur[:fit_t] = testFit(dfVCur[:vaf_t], dfVCur[:vafTh_t], dfVCur[:cov_t], errorQuantile)
        dfVCur[:goodFit] = all(dfVCur[:fit_t])
    end
end

"""
    Fits the input variant sizes `vaf_t_Vid` at times `_t`, measured with coverage `nCov_t` to the models specified in `models`. All logistic models transition between `0` and `0.5`.
    Returns vector of fitted selective advantages `s_vid`. Fits that do not succeed are given a `missing` value.
"""
function fitModelsLightweight(vaf_t_Vid, _t, cov_t, Nf::Real, models::Vector{String}=["positive logistic", "negative logistic", "constant"]; errorQuantile=0.99 )
    s_vid = Vector{Union{Missing,Float64}}(undef, length(vaf_t_Vid))
    for (i,vaf_t) in enumerate(vaf_t_Vid)
        _models = [createModelFromName(model) for model in models]
        llh_model = [fitGrowthModel!(_t, vaf_t, cov_t, Nf, model)
            for model in _models]
        llhMax, indMax = findmax(llh_model)
        modelFit = _models[indMax]
        s_vid[i] = modelFit.shape.γ
        vafThQl_t, vafThQu_t, fit_t = testFit(vaf_t, predict(_t, modelFit.shape), cov_t, errorQuantile)
        !all(fit_t) && (s_vid[i] = missing)
        # println(s_vid[i])
    end
    return s_vid
end

"""
    Fits the input variant sizes `vaf_t` at times `_t`, measured with coverage `nCov_t` to the model specified in `model`. 
"""
function fitGrowthModel!(_t, vaf_t, nCov_t, Nf, model::ModelFit; timeLimit=1.0)
    
    # round because X-chromosome mutations can be at half integer abundance
    X_t = round.(vaf_t.*nCov_t) .|> Int #! we need a better correction for the X chromosome later on

    initial_β = setInitialParamsMLE!(model, _t, vaf_t, Nf)

    function logprobX(X, t, n, β)
        p = fLogModel(t,β,model)
        if p <= 0 || p > 1 || isnan(p)
            p = 1E-20
        end
        logpdf(Binomial(n, p), X)
    end
    logllh(β, p) = 
        -([logprobX(X_t[i], _t[i], nCov_t[i], β) for i in eachindex(_t)] |> sum)

    lcons, ucons = boundsMLE(model, _t, vaf_t)
    optprob = OptimizationFunction(logllh, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optprob, initial_β, lb=lcons, ub=ucons)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
    addFitParams!(model, sol.u)
    logllhFit =  - logllh(sol.u, 0)
    return logllhFit
end

function testFit!(dfVid::DataFrame, errorQuantile)
    nVars = size(dfVid,1)
    if !hasproperty(dfVid, :vafThQl_t)
        dfVid.vafThQl_t = Vector{Vector{Float64}}(undef, nVars)
    end
    if !hasproperty(dfVid, :vafThQu_t)
        dfVid.vafThQu_t = Vector{Vector{Float64}}(undef, nVars)
    end
    if !hasproperty(dfVid, :goodFit)
        dfVid.goodFit = Vector{Bool}(undef, nVars)
    end
    if !hasproperty(dfVid, :_goodFit)
        dfVid._goodFit = Vector{Vector{Bool}}(undef, nVars)
    end

    for dfVCur in eachrow(dfVid)
        dfVCur[:vafThQl_t], dfVCur[:vafThQu_t], dfVCur[:_goodFit] = testFit(dfVCur[:vaf_t], dfVCur[:vafTh_t], dfVCur[:cov_t], errorQuantile)
        dfVCur[:goodFit] = all(dfVCur[:_goodFit])
    end
end

function testFit(vafS_t, vafTh_t, cov_t, errorQuantile)
    vafThQL_t = [
            quantile(
                Binomial(cov_t[tInd], vafTh_t[tInd]), 1-errorQuantile
                ) / cov_t[tInd]
            for tInd in eachindex(cov_t)
        ]
    vafThQU_t = [
            quantile(
                Binomial(cov_t[tInd], vafTh_t[tInd]), errorQuantile
                ) / cov_t[tInd]
            for tInd in eachindex(cov_t)
        ]
    fit_t = ((vafS_t.>vafThQL_t) .&& (vafS_t.<vafThQU_t))
    return vafThQL_t, vafThQU_t, fit_t
end



end