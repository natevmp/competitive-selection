
struct VFreqspace{F<:AbstractFloat}
    _f::Vector{F}
    n_f::Vector{F}
end

"Create a discretized frequency space with variable grid spacing.
N   == size of the population
l   == number of gridpoints in the discretization (including f=0 and f=1)
K   == offset of f=1/N gridpoint; K=0 implies 1/N is at i=2, K=1 means 1/N is at i=3, etc.
"
function VFreqspace(N::Integer, l::Integer, K::Integer=4)
    _f = zeros(Float64, l)
    _f[1] = 0
    for (i,k) in enumerate(K:-1:1)
        _f[1+i] = 1/N * (1/2)^k
    end

    M = l - K
    a = 2 * ( N-(M-1) )/( N*(M-1)*(M-2) )
    df_m = [1/N + (i-1)*a for i in 1:M-1]
    for m in 1:M-1
        _f[K+1+m] = sum(df_m[1:m])
    end

    n_f = zeros(Float64, l)
    VFreqspace(_f, n_f)
end

Base.length(vfs::VFreqspace) = length(vfs._f)

"""Evolve growing population VAF spectrum with Moran and pure birth dynamics"""
function evolveExistingVariants(vfs::VFreqspace, params::Dict, t::Real; alg=nothing, order=2, reltol=1e-6)
    N = params[:N]
    α = params[:α]
    n₀ = params[:n₀]

    fluxInd = findall(f->f==1/N, vfs._f)[1] - 1

    inL = length(vfs)-2
    dx_i = (vfs._f[2:end] .- vfs._f[1:end-1])
    in_x = vfs._f[2:end-1]
    
    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates

    Δ = CenteredDifference(2, order, dx_i, inL)
    B = spzeros(inL, inL)
    fB(x) = α/N * x*(1-x)
    Bu_x = similar(vfs.n_f[2:end-1])
    ΔBu_x = similar(vfs.n_f[2:end-1])

    function step!(du_x, u_x, (B, Bu_x, ΔBu_x), t)
        B[diagind(B)] .= (x -> fB(x)).(in_x)
        mul!(Bu_x, B, u_x)
        mul!(ΔBu_x, Δ*BC, Bu_x)
        du_x .= ΔBu_x
    end

    t0 = 0.0
    t1 = t
    n₀_x = copy(vfs.n_f[2:end-1])
    n₀_x[fluxInd] = n₀

    prob = ODEProblem(step!, n₀_x, (t0, t1), (B, Bu_x, ΔBu_x))
    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f[2:end-1] .= sol.u[end]
end


"""Evolve fixed sized population VAF spectrum with Moran dynamics as a diffusion PDE"""
function evolveSizeDist(vfs::VFreqspace, params::Dict, t::Real;  alg=nothing, order=2, reltol=1e-6)

    N = params[:N]
    α = params[:α]
    s = params[:s]
    μ = params[:μ]
    
    inL = length(vfs._f)-2
    df_f = (vfs._f[2:end] .- vfs._f[1:end-1])
    _fIn = vfs._f[2:end-1]

    n0_f = copy(vfs.n_f[2:end-1])
    fluxInd = findfirst(_fIn .== 1/N)
    dϵ = (df_f[fluxInd] + df_f[fluxInd+1])/2 # scales intensity of flux delta function from frequency step
    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates

    ∇ = CenteredDifference(1, order, df_f, inL)
    Δ = CenteredDifference(2, order, df_f, inL)
    c_f = spzeros(inL)
    A = spzeros(inL, inL)
    B = spzeros(inL, inL)
    fA(x, xTot) = - α*s*x*( 1 - xTot )
    fB(x, xTot) = α/(2N)*x*( 2 + s + s*xTot )

    An_f = similar(n0_f)
    Bn_f = similar(n0_f)

    ∇An_f = similar(n0_f)
    ΔBn_f = similar(n0_f)

    function step!(dn_f, n_f, (A, B, An_f, Bn_f, ∇An_f, ΔBn_f, c_f), t)
        n_f[n_f.<0] .= 0
        nSpl = Spline1D(_fIn, n_f; k=2, bc="zero")
        # xTot, e = quadgk((x -> x*nSpl(x)),0,1; order=10)
        xTot, e = quadgk((x -> x*nSpl(x)),0,1)
        
        A[diagind(A)] .= (x -> fA(x,xTot)).(_fIn)
        B[diagind(B)] .= (x -> fB(x,xTot)).(_fIn)
        c_f[fluxInd] = μ / dϵ
        # c_f[fluxInd] = μ

        mul!(An_f, A, n_f)
        mul!(Bn_f, B, n_f)
        mul!(∇An_f, ∇*BC, An_f)
        mul!(ΔBn_f, Δ*BC, Bn_f)
        dn_f .= ∇An_f .+ ΔBn_f .+ c_f
    end

    t0 = 0.0
    t1 = t

    if μ==0
        n0_f[fluxInd] = 1 / dϵ
        println(1/ dϵ)
    end

    prob = ODEProblem(step!, n0_f, (t0, t1), (A, B, An_f, Bn_f, ∇An_f, ΔBn_f, c_f))

    sol = solve(prob, save_everystep=false)

    # vfs.n_f[2:end-1] .= sol.u[end] .* vfs._f[2:end-1]
    vfs.n_f[2:end-1] .= sol.u[end]

end

function evolveDriftPDE(vfs::VFreqspace, params::Dict, t::Real;  alg=nothing, order=2, reltol=1e-6)
    N = params[:N]
    α = params[:α]
    s = params[:s]
    μ = params[:μ]

    inL = length(vfs._f)-2
    df_f = (vfs._f[2:end] .- vfs._f[1:end-1])
    _f = vfs._f[2:end-1]
    n0_f = copy(vfs.n_f[2:end-1])
    fluxInd = findfirst(_f .== 1/N)
    dϵ = (df_f[fluxInd] + df_f[fluxInd+1])/2 # scales intensity of flux delta function from frequency step
    BC = Dirichlet0BC(Float64)
    ∇ = CenteredDifference(1, order, df_f, inL)
    c_f = spzeros(inL)
    ∇n_f = similar(n0_f)
    if μ==0
        n0_f[fluxInd] = 1 / dϵ
    end

    # ∂_t n(x,t) = α(1+s)( ∱y⋅n(y,t)dy - 1 )( n(x,t) + x ∂_x n(x,t) )
    # ∂ₜn(x,xTot,n,∂ₓn,c) = α*(1+s)*(xTot-1)*(n + x*∂ₓn) + c
    ∂ₜn(x,xTot,n,∂ₓn,c) = α*(1+s)*(xTot-1)*(n + x*∂ₓn) + c
    function step!(dn_f, n_f, (c_f, ∇n_f), t)
        # n_f[n_f.<0] .= 0    # jank enforce positive image
    
        # nSpl = Spline1D(_f, n_f; k=2)
        # xTot, e = quadgk((x -> x*nSpl(x)),0,1; order=7)
        c_f[fluxInd] = μ / dϵ
        mul!(∇n_f, ∇*BC, n_f)

        dn_f .= ∂ₜn.(_f, 0.9, n_f, ∇n_f, c_f)
    end

    t0 = 0.0
    t1 = t
    prob = ODEProblem(step!, n0_f, (t0, t1), (c_f, ∇n_f))

    sol = solve(prob, save_everystep=false)
    vfs.n_f[2:end-1] .= sol.u[end]

end

function evolveExpectedGrowth()

    # create discretized grid of u (variant birth time)

    # set up evolver for x_u(t) and y(t)

    # set initial values for x_u(t) and y(t)




end