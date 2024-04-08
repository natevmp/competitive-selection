
function prepareSims(params, model, _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}, runs) where {U,S<:Real}
    T = params[:T]
    μ = params[:μ]
    simArgs = DataFrame()
    simArgs.k = rand(Poisson(μ*T), runs)    # total number of variants arising
    simArgs.t₀_vid = (k->(T*rand(k) |> sort!)).(simArgs[!,:k])
    simArgs._trackerID = [Vector{Int}(undef,length(_trackerVariant)) for j in 1:runs]
    for row in eachrow(simArgs)
        for ts in _trackerVariant # `ts` is tuple (t0, s)
            row.k += 1
            push!(row.t₀_vid, ts[1])
        end
        sort!(row.t₀_vid)
        for (i,ts) in enumerate(_trackerVariant)
            row._trackerID[i] = findfirst(row.t₀_vid .== ts[1])
        end
    end
    simArgs.s_vid = (k->selectionParamModel(model, k)).(simArgs[!,:k])
    if haskey(params, :sMax)
        for s_vid in simArgs.s_vid
            maxS_vid = findall(s_vid.>params[:sMax])
            s_vid[maxS_vid] .= params[:sMax]
        end
    end
    simArgs.init_vid = falses.(simArgs[!,:k])
    simArgs.x₀_vid = zeros.(simArgs[!,:k])
    simArgs.parentId_vid = zeros.(Int64, simArgs[!,:k])

    length(_trackerVariant)<1 && (return simArgs)
    eltype(_trackerVariant)==Tuple{Real} && (return simArgs)
    for rowArgs in eachrow(simArgs)
        rowArgs.s_vid[rowArgs._trackerID] .= [ts[2] for ts in _trackerVariant]
    end
    return simArgs
end

function selectionParamModel(model::FixedSelectionModel, k::Int)
    fill(model.s, k)
end

function selectionParamModel(model::FreeFixedModel, k::Int)
    fill(model.s, k)
end

function selectionParamModel(model::GaussianSelectionModel, k::Int)
    rand(Normal(model.s, model.σ), k)
end

function selectionParamModel(model::ExponentialSelectionModel, k::Int)
    rand(Exponential(model.s), k)
end

function selectionParamModel(model::GammaSelectionModel, k::Int)
    rand(Gamma(model.s^2/model.σ^2 , model.σ^2/model.s), k)
end

function fitnessDoubleHit(sChild, sParent)
    sParent + sChild + sChild*sParent
end

"""
    Sets the parent clone of a new variant, and optionally performs a double hit if `q>0`.
"""
function initiateVariant!(childId::Integer, init_vid, s_vid, x_vid, parentId_vid::Vector{Int64}, q; sMax::Real=Inf)
    pRand = rand()  # random variable ∈(0,1) to select the parent
    xTot = 0.
    for (parentId,x) in enumerate(@view x_vid[init_vid])
        xTot += x
        if xTot<pRand continue end # continue if current ID is not parent
        parentId_vid[childId] = parentId
        if (q==0 || rand()>q) break end # break if double hit fails
        s_vid[childId] = min(fitnessDoubleHit(s_vid[childId], s_vid[parentId]), sMax)
        break
    end
    init_vid[childId] = true
end

function drift(model::Model, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t; sMax::Real=Inf)
    sX = 0   # non-allocating sum
    for i in eachindex(dx_vid)
        t < t0_vid[i] && break
        sX += x_vid[i]*s_vid[i]
    end
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break # skip variant if not yet initiated
        !init_vid[i] && initiateVariant!(i, init_vid, s_vid, x_vid, parentId_vid, model.q; sMax=sMax) # initiate variant if new
        if x<0
            x_vid[i] = 0
            dx_vid[i] = 0
            continue
        end
        if x>1
            x_vid[i] = 1
            dx_vid[i] = 0
            continue
        end
        dx_vid[i] = α*x*(s_vid[i] - sX) # evolve variant
    end
end

function diffusion(model::Model, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t)
    sX = 0  # non-allocating sum
    for i in eachindex(dx_vid)
        t < t0_vid[i] && break
        sX += x_vid[i]*s_vid[i]
    end
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        dx_vid[i] = (x>0 && x<1) ? √( α*(2 + s_vid[i] + sX)*x/N ) : 0
    end
end

function drift(model::FreeFixedModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t; sMax::Real=Inf)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        !init_vid[id] && initiateVariant!(i, init_vid, s_vid, x_vid, parentId_vid, model.q)
        if x>0
            dx_vid[i] = α*s_vid[i]*x*(1 - x)
        else
            x_vid[i] = 0
            dx_vid[i] = 0
        end
    end
end

function diffusion(model::FreeFixedModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        dx_vid[i] = (x>0) ? √( α*(2+s_vid[i]*(1+x))*x/N ) : 0
    end
end

"""
    evolvePopSim()

`_trackerVariant` is a vector containing clones to be tracked. The elements are either the birth times, or a tuple containing both birth time and fitness of the form `(t0, s)`.
"""
function evolvePopSim(params::Dict; runs::Int=1, _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}=Vector{Float64}[], noDiffusion::Bool=false) where {U,S<:Real}
    model =
        if params[:sType]=="fixed"
            FixedSelectionModel(params[:s], params[:q])
        elseif params[:sType]=="exponential"
            ExponentialSelectionModel(params[:s], params[:q])
        elseif params[:sType]=="gaussian"
            GaussianSelectionModel(params[:s], params[:σ], params[:q])
        elseif params[:sType]=="gamma"
            GammaSelectionModel(params[:s], params[:σ], params[:q])
        elseif params[:sType]=="free"
            FreeFixedModel(params[:s], params[:q])
        else
            println("error: selection model undefined.")
        end
    evolvePopSim(params, model; runs, _trackerVariant, noDiffusion)
end

function evolvePopSim(params::Dict, model::Model; runs::Int=1, _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}=Vector{Float64}[], noDiffusion::Bool=false) where {U,S<:Real}
    
    N = params[:N]
    α =
        if haskey(params, :α)
            params[:α]
        else
            1/params[:τ]
        end
    s = params[:s]
    T = params[:T]
    μ = params[:μ]

    f!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = 
        drift(model, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t; sMax=(haskey(params, :sMax) ? params[:sMax] : Inf))
    g!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = begin
        noDiffusion ? 0 : diffusion(model, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t)
    end

    _t = range(0,T, length=T+1)

    simArgs = prepareSims(params, model, _trackerVariant, runs)
    # t0Min = minimum([t_vid[1] for t_vid in simArgs[!,:t₀_vid]])
    # get allocation-free minimum
    t0Min = 1.
    for t_vid in simArgs[!,:t₀_vid]
        length(t_vid)==0 && continue # ingnore empty arrays
        t_vid[1]<t0Min && (t0Min = t_vid[1])
    end

    condt0(u,t,integrator) = t==t0Min
    function addStops!(integrator)
        for t0 in integrator.p[1]
            add_tstop!(integrator, t0)
        end
    end
    callBackAddStops = DiscreteCallback(condt0, addStops!; save_positions=(false,false))
    function affect!(integrator)
        i = findfirst(integrator.p[1].==integrator.t)
        integrator.u[i] = 1/integrator.p[5]
    end
    condition(u,t,integrator) = t ∈ integrator.p[1]
    callBackAddClone = DiscreteCallback(condition, affect!; save_positions=(false,false))
    callbacks = CallbackSet(callBackAddStops, callBackAddClone)

    prob = SDEProblem(f!, g!, simArgs[1,:x₀_vid], (0.0, T), [simArgs[1,:t₀_vid], simArgs[1,:init_vid], α, simArgs[1,:s_vid], N, simArgs[1,:parentId_vid]])
    # solver=SOSRI2
    # solver = SOSRA2
    # solver = SRA3
    # solver = ImplicitEulerHeun
    solver = LambaEM    #fastest
    # solver = SKenCarp #fails
    function probFunc(prob, i, repeat)
        remake(prob, u0=simArgs[i,:x₀_vid], p=[simArgs[i,:t₀_vid], simArgs[i,:init_vid], α, simArgs[i, :s_vid], N, simArgs[i,:parentId_vid]])
    end
    ensembleProb = EnsembleProblem(prob, prob_func=probFunc)
    # return ensembleProb
    solEns = solve(ensembleProb, solver(), EnsembleThreads(); callback=callbacks, tstops=[t0Min,], saveat=_t, trajectories=runs)
    return solEns, simArgs
end

# function evolvePopRand(params::Dict)

#     N = params[:N]
#     α = params[:α]
#     s = params[:s]
#     T = params[:T]
#     μ = params[:μ]
#     σ = params[:σ]
#     # k = params[:k]

#     _t = range(0,T, length=101)
#     k = rand(Poisson(μ*T))
#     println(k)
#     t₀_vid = T*rand(k) |> sort!
#     x₀_vid = ones(length(t₀_vid)) * 1/N
#     s_vid = rand(Normal(μ, σ), k)
#     prob = SDEProblem(drift, diffusion, x₀_vid, (0.0, T), (t₀_vid, α, s_vid, N))
#     sol = solve(prob, SOSRI2(), saveat=_t)
#     # ensembleProb = EnsembleProblem(prob)
#     # solEns = solve(ensembleProb, SOSRI2(), EnsembleThreads(), trajectories=runs, saveat=_t)
#     return sol
# end


# function evolveClonesDiff(model::BiCloneModel, n0::Real, _t::AbstractVector{T} where T, runs::Int; onlySurvivors::Bool=false)
#     f(du, u, (N,ρ,s), t) = drift(du, u, (N,ρ,s), t, model)
#     g(du, u, (N,ρ,s), t) = diffusion(du, u, (N,ρ,s), t, model)
#     x₀=[n0, model.N-n0]
#     # x₀=[model.N*f0, model.N-model.N*f0]
#     prob = SDEProblem(f, g, x₀, (0.0, _t[end]), (model.N,model.ρ,model.s))
#     ensembleProb =
#         if onlySurvivors
#             EnsembleProblem(prob, output_func=conditionSurvival)
#         else
#             EnsembleProblem(prob)
#         end
#     solEns = solve(ensembleProb, SOSRI2(), EnsembleThreads(), trajectories=runs, saveat=_t)
#     # solEns = solve(ensembleProb, SOSRI(), EnsembleThreads(), trajectories=runs)
# end
