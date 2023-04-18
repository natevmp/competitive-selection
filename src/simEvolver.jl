

function doubleHitVariant!(vid::Integer, init_vid, s_vid, x_vid, q)
    if q==0  # if q==0 we skip the double hit process
        init_vid[vid] = true
        return
    end
    if rand()>q # if double hit doesn't occur we break
        init_vid[vid] = true
        return
    end
    pRand = rand()
    xTot = 0.
    for (i,x) in enumerate(x_vid[init_vid])
        xTot += x
        if pRand<xTot
            s_vid[vid] += s_vid[i] + s_vid[vid]*s_vid[i]
            break
        end
    end
    init_vid[vid] = true
end


function drift(model::Model, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N), t)
    sX = 0   # non-allocating sum
    for i in eachindex(dx_vid)
        t < t0_vid[i] && break
        sX += x_vid[i]*s_vid[i]
    end
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        !init_vid[i] && doubleHitVariant!(i, init_vid, s_vid, x_vid, model.q)
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
        dx_vid[i] = α*x*(s_vid[i] - sX)
    end
end

function diffusion(model::Model, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N), t)
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

function drift(model::FreeFixedModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N), t)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        !init_vid[id] && doubleHitVariant!(i, init_vid, s_vid, x_vid, model.q)
        if x>0
            dx_vid[i] = α*s_vid[i]*x*(1 - x)
        else
            x_vid[i] = 0
            dx_vid[i] = 0
        end
    end
end

function diffusion(model::FreeFixedModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N), t)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        dx_vid[i] = (x>0) ? √( α*(2+s_vid[i]*(1+x))*x/N ) : 0
    end
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

function prepareSims(params, model, trackerVariant, runs)
    T = params[:T]
    μ = params[:μ]
    simArgs = DataFrame()
    simArgs.k = rand(Poisson(μ*T), runs)
    simArgs.t₀_vid = (k->(T*rand(k) |> sort!)).(simArgs[!,:k])
    if !isnothing(trackerVariant)
        simArgs.trackerID = Vector{Int}(undef, runs)
        for row in eachrow(simArgs)
            row.k += 1
            push!(row.t₀_vid, trackerVariant[1])
            sort!(row.t₀_vid)
            row.trackerID = findfirst(row.t₀_vid .== trackerVariant[1])
        end
    end
    simArgs.s_vid = (k->selectionParamModel(model, k)).(simArgs[!,:k])
    simArgs.init_vid = falses.(simArgs[!,:k])
    simArgs.x₀_vid = zeros.(simArgs[!,:k])

    if isnothing(trackerVariant)
        return simArgs
    end
    
    length(trackerVariant)==1 && (return simArgs)

    for rowArgs in eachrow(simArgs)
        rowArgs.s_vid[rowArgs.trackerID] = trackerVariant[2]
    end
    println(simArgs[1, :s_vid])
    println(simArgs[1,:s_vid][simArgs[1,:trackerID]])
    println()
    return simArgs
end

function evolvePopSim(params::Dict; runs::Int=1, trackerVariant::Union{Nothing, Real,Tuple{Real,Real}}=nothing, noDiffusion::Bool=false)
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
    evolvePopSim(params, model; runs, trackerVariant, noDiffusion)
end

function evolvePopSim(params::Dict, model::Model; runs::Int=1, trackerVariant::Union{Nothing, Real,Tuple{Real,Real}}=nothing, noDiffusion::Bool=false)
    
    N = params[:N]
    α = params[:α]
    s = params[:s]
    T = params[:T]
    μ = params[:μ]

    f!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N), t) = 
        drift(model, dx_vid, x_vid, (t0_vid, init_vid, α, s, N), t)
    g!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N), t) = begin
        noDiffusion ? 0 : diffusion(model, dx_vid, x_vid, (t0_vid, init_vid, α, s, N), t)
    end

    _t = range(0,T, length=T+1)

    simArgs = prepareSims(params, model, trackerVariant, runs)
    # t0Min = minimum([t_vid[1] for t_vid in simArgs[!,:t₀_vid]])
    # get allocation-free minimum, ingnoring empty arrays
    t0Min = 1.
    for t_vid in simArgs[!,:t₀_vid]
        length(t_vid)==0 && continue
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

    prob = SDEProblem(f!, g!, simArgs[1,:x₀_vid], (0.0, T), [simArgs[1,:t₀_vid], simArgs[1,:init_vid], α, simArgs[1,:s_vid], N])
    # solver=SOSRI2
    # solver = SOSRA2
    # solver = SRA3
    solver = LambaEM    #fastest
    # solver = SKenCarp #fails
    function probFunc(prob, i, repeat)
        remake(prob, u0=simArgs[i,:x₀_vid], p=[simArgs[i,:t₀_vid], simArgs[i,:init_vid], α, simArgs[i, :s_vid], N])
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
