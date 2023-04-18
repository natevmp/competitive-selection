module CompetitiveSelection

# using StaticArrays
using ElasticArrays
using DataFrames
using Statistics, StatsBase, LinearAlgebra, SparseArrays
using DifferentialEquations
using DiffEqOperators, Dierckx, QuadGK
using LsqFit, StatsBase
using Distributions, Optim, LineSearches
using Random
include("dataStructuring.jl")
using .DataStructuring

# ---- API -----
# theory function exports
export fLogistic
# simulation and PDE exports
export VFreqspace, evolvePopSim
# analysis exports
export fitLogGrowth


# struct Variant{T<:Real}
#     x_t::Vector{T}
#     _t::Vector{T}
# end
# Variant(N::Real) = Variant([n0,], [0.,])

abstract type Model end

abstract type StochasticSelectionModel <: Model end

struct GaussianSelectionModel <: StochasticSelectionModel
    s::Float64
    σ::Float64
    q::Float64
end

struct ExponentialSelectionModel <: StochasticSelectionModel
    s::Float64
    q::Float64
end

struct GammaSelectionModel <: StochasticSelectionModel
    s::Float64
    σ::Float64
    q::Float64
end

struct FixedSelectionModel <: Model
    s::Float64
    q::Float64
end

struct FreeFixedModel <: Model
    s::Float64
    q::Float64
end

include("simEvolver.jl")
include("pdeEvolver.jl")
include("scientist.jl")


end

