
using Pkg

# Pkg.add("JLD2")
# Pkg.add("ProgressMeter")
# Pkg.add("Random")
# Pkg.add("Distributions")
# Pkg.add("DifferentialEquations")
# Pkg.add("DiffEqOperators")
# Pkg.add("OrdinaryDiffEq")
# Pkg.add("LinearAlgebra")
# Pkg.add("StatsBase")
# Pkg.add("Statistics")
# Pkg.add("SparseArrays")
# Pkg.add("StaticArrays")
# Pkg.add("Dierckx")
# Pkg.add("QuadGK")
# Pkg.add("Roots")
# Pkg.add("Interpolations")

Pkg.add(["CSV", "DataFrames", "Statistics", "LsqFit", "StatsBase", "Distributions", "Optim", "Optimization", "OptimizationBBO", "OptimizationMOI", "Ipopt", "ProgressMeter"])

Pkg.add(["ElasticArrays", "LinearAlgebra", "SparseArrays", "DifferentialEquations", "DiffEqOperators", "Dierckx", "QuadGK", "LineSearches", "Random"])

Pkg.add(["GpABC","Distances"])

Pkg.add("JLD2")