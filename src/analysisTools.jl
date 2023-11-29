module AnalysisTools

using CSV
using DataFrames
using Statistics
using LsqFit, StatsBase

# export fitLogGrowth
export fLogistic, fIndFitnessWithAvCompetition, fConstant

# function fitLogGrowth(_t, vaf_t, Nf)
#     fLogistic(t, t0, r) = 1 / ( 1 + (1 - 1/Nf)/(1/Nf) * exp(-r*(t-t0)) )
#     fAntiLog(t, t0, r, b0) = 1 / ( 1 + (1-b0)/b0 * exp(-r*(t-t0)) )
#     @. modelLogistic(t, p) = fLogistic(t, p[1], p[2])

#     # b0 /(b0 + (1 - b0) E^(-t \[Lambda]))

#     @. modelAntiLog(t, p) = fAntiLog(t, p[1], p[2], p[3])
#     corr = corspearman(_t, vaf_t)
#     fit = let 
#         if corr >= 0
#             curve_fit(modelLogistic, _t, vaf_t, [1., 0.2], lower=[-Inf,0], upper=[_t[1], +Inf])
#         else
#             curve_fit(modelAntiLog, _t, vaf_t, [1.,-0.1, 0.05], lower=[-Inf, -Inf, 0], upper=[_t[1], 0, 1])
#         end
#     end
#     return fit
# end

fLogistic(t, t0, r, x0) = 1 / ( 1 + (1-x0)/x0 * exp(-r*(t-t0)) )

fLogistic(t, t0, r, x0, xF) = xF / ( 1 + (xF-x0)/x0 * exp(-r*(t-t0)) )

# fLogistic(t, t0, r, x0, xF, z) = (xF-z/r) / ( 1 + (xF-z/r-x0)/x0 * exp(-(r-z/xF)*(t-t0)) )

function fLogistic(t, t0, r, x0, xF, z)
    if z < r*xF - 0.001
        return (xF-z/r) / ( 1 + (xF-z/r-x0)/x0 * exp(-(r-z/xF)*(t-t0)) )
    elseif z > r*xF + 0.001
        return (xF-z/r) / ( 1 + (xF-z/r-x0)/x0 * exp(-(r-z/xF)*(t-t0)) )
    else
        return xF*x0 / (xF + r*x0*t)
    end
end

fConstant(t,x) = x

fIndFitnessWithAvCompetition(t, t0, s, a, μ, N) = exp( -(a-s)*(t-t0) )*( a*exp(-(a+μ/N)*t0) +  μ/N ) /
    ( N * (a*exp(-(a+μ/N)*t) + μ/N) )

"""
get size binned histogram from list of variant sizes.
"""
function distbin(x_vid; bins::Int=25, xMin::Real=0, xMax::Union{Nothing,Real}=nothing, normalized::Bool=true)
    isnothing(xMax) && (xMax = maximum(x_vid)+0.01)
    _xEdge = range(xMin, xMax, bins+1)
    dx = Float64(_xEdge.step)
    # _x = _xEdges[1:end-1] .+ dx/2
    inRange_vid = (x -> x>=xMin && x<xMax).(x_vid)   #apply thresholds
    n_x = zeros(Float64, bins)
    # for x in @view x_vid[inRange_vid]
    for x in x_vid[inRange_vid]
        k = ((x-xMin)/dx |> floor |> Int) + 1
        n_x[k] += 1
    end
    normalized && (n_x .*= dx/sum(inRange_vid))
    return _xEdge, n_x
end

"""
Compute the size distribution of `f_vid` in absolute number of elements per bin.
"""
function sizeDistAbs(x_vid; bins::Int=25, xMin=0., xMax=0.5)
    _xEdge = range(xMin, xMax, bins+1)
    dx = Float64(_xEdge.step)
    n_x = zeros(Int64, bins)
    for x in x_vid
        k = ((x-xMin)/dx |> floor |> Int) + 1
        n_x[k] += 1
    end
    return _xEdge, n_x
end

"""
Compute the size distribution of `f_vid` as a density per bin.
"""
function sizeDistDens(x_vid; bins::Int=25, xMin=0., xMax=0.5)
    _xEdge = range(xMin, xMax, bins+1)
    dx = Float64(_xEdge.step)
    n_x = zeros(Float64, bins)
    for x in x_vid
        k = ((x-xMin)/dx |> floor |> Int) + 1
        if k > bins
            n_x[end] += 1.
            continue
        end
        n_x[k] += 1.
    end
    nDens_x = (n -> n*dx/sum(n_x)).(n_x)
    # _fBinCentre = _fBinEdge[1:end-1] .+ Float64(_fBinEdge.step)/2
    return _xEdge, nDens_x
end

function sizeDistDens(x_vid, t_vid, tBounds; bins::Int=25, xMin=0., xMax=0.5)
    tMask_vid = ( t -> t>tBounds[1] && t<tBounds[2] ).(t_vid)
    _xEdge, nDens_x = sizeDistDens(@view x_vid[tMask_vid]; bins, xMin, xMax)
    return _xEdge, nDens_x
end

    
end