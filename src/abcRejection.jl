module ABCRejection

using Random
using Distributions
include("competitiveSelection.jl")
using .CompetitiveSelection
include("analysisTools.jl")
using .AnalysisTools
include("dataStructuring.jl")
using .DataStructuring
using DifferentialEquations
using Distances

export runABCParticles, runParticle, testABCParticles, testParticle

## ====================== ABC framework ======================
struct Particle{T}
    paramSet::NamedTuple
    simResults::T
end

function runABCParticles(runModelSim::Function, priorDists_pid::NamedTuple, nParticles::Integer, ctrlParams::Dict; verbose::Bool=false)
    particle_tid = Vector{Particle}(undef, nParticles)
    for tid in 1:nParticles
        particle_tid[tid] = runParticle(runModelSim, priorDists_pid, ctrlParams; verbose)
    end
    return particle_tid
end

function runParticle(runModelSim::Function, pDist_pid::NamedTuple, ctrlParams::Dict; verbose::Bool=false)
    pVal_pid = [rand(pDist) for pDist in pDist_pid]
    if verbose
        println("parameter values of particle:")
        for (i,pid) in enumerate(propertynames(pDist_pid))
            println(string(pid)*": ", string(pVal_pid[i]))
        end
    end
    paramSet = (; zip(keys(pDist_pid), pVal_pid)...)
    simResults = runModelSim(paramSet, ctrlParams)
    return Particle(
        paramSet,
        simResults
    )
end

function testABCParticles(compareDataVSim::Function, particle_tid::Vector, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
    accepted_tid = falses(length(particle_tid))
    for (tid,particle) in enumerate(particle_tid)
        accepted_tid[tid] = testParticle(compareDataVSim, particle.simResults, dataMetrics, errorThresholds, constraintsF)
    end
    return accepted_tid
end

function testParticle(compareDataVSim::Function, simResults, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
    acceptedMetrics = compareDataVSim(dataMetrics, simResults, errorThresholds)
    if isnothing(constraintsF)
        return all(acceptedMetrics)
    else
        acceptedConstraints = constraintsF(simResults)
        return all(acceptedMetrics) && all(acceptedConstraints)
    end
end

function measureABCParticles(compareDataVSimError::Function, particle_tid::Vector, dataMetrics)
    dists1 = compareDataVSimError(dataMetrics, particle_tid[1].simResults)
    dists_tid = Vector{typeof(dists1)}(undef, length(particle_tid))
    dists_tid[1] = dists1
    for (tid,particle) in enumerate(@view particle_tid[2:end])
        dists_tid[tid] = compareDataVSimError(dataMetrics, particle.simResults)
    end
    return dists_tid
end

function runABC(runModelSim::Function, compareDataVSim::Function, priorDists_pid::NamedTuple, dataMetrics, ctrlParams::Dict, nParticles)

    # priorDists_pid = (s=Uniform(1,2),σ=Uniform(1,2),μ=Uniform(1,2))
    particle_tid = runABCParticles(runModelSim, priorDists_pid, ctrlParams, nParticles)

    accepted_tid = testABCParticles(compareDataVSim, particles_tid, dataMetrics, ctrlParams[:thresholds])

    return particle_tid, accepted_tid
end

function acceptedParams(particle_tid, accepted_tid)
    [particle.paramSet for particle in particle_tid[accepted_tid]]
end

end